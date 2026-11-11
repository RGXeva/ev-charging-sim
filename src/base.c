/* SPDX-License-Identifier: MIT */
/*
 * Base station: a single MPI process; one thread performs all MPI calls.
 *
 * Versions of the receive loop (all poll periodically so the base can also
 * decide when to terminate the simulation):
 *   B1  MPI_Iprobe, then blocking MPI_Recv; replies with blocking MPI_Send
 *   B2  one pre-posted MPI_Irecv + MPI_Test; replies with MPI_Isend
 *   B3  R persistent receives (MPI_Recv_init/MPI_Startall) + MPI_Testsome,
 *       a whole batch of reports is handled per poll; replies with MPI_Isend
 */
#define _POSIX_C_SOURCE 200809L
#include "common.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <unistd.h>

#define REPORT_BYTES (REPORT_DOUBLES * (int)sizeof(double) + REPORT_INTS * (int)sizeof(int))
#define B3_SLOTS 16

typedef struct { double *v; int n, cap; } dvec_t;

typedef struct {
    const config_t *cfg;
    int N;
    MPI_Datatype rtype;
    double t0;
    double *last_alert, *busy_until;    /* indexed by cart rank, base clock (s since t0) */
    FILE *log, *acsv, *scsv;
    long alerts, r_cand, r_none_near, r_none_grid, unexpected;
    long msg_recv, msg_sent, bytes_recv, bytes_sent;
    double svc_sum, svc_max;
    int finals;
    double *nst;                        /* N x ST_COUNT final node statistics */
    int term_sent;
    double t_term, t_last_final;
    int (*rbuf)[RD_MAX_INTS];           /* one reply buffer per node (<= 1 alert in flight) */
    MPI_Request *rreq, *treq;
    int term_val;
    dvec_t s_nbr, s_alert, s_net, s_sync, s_svc;
} base_t;

static void dv_push(dvec_t *d, double x)
{
    if (d->n == d->cap) {
        int cap = d->cap ? 2 * d->cap : 256;
        double *p = realloc(d->v, (size_t)cap * sizeof(double));
        if (!p) return;
        d->v = p;
        d->cap = cap;
    }
    d->v[d->n++] = x;
}

static int cmpd(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

/* nearest-rank percentile; vector is sorted in place */
static double pct(dvec_t *d, double p)
{
    if (d->n == 0) return NAN;
    qsort(d->v, d->n, sizeof(double), cmpd);
    int i = (int)ceil(p * d->n) - 1;
    if (i < 0) i = 0;
    if (i >= d->n) i = d->n - 1;
    return d->v[i];
}

static double mean(const dvec_t *d)
{
    if (d->n == 0) return NAN;
    double s = 0;
    for (int i = 0; i < d->n; i++) s += d->v[i];
    return s / d->n;
}

static void now_str(char *b, size_t len)
{
    int ts[7];
    stamp_now(ts, NULL);
    stamp_fmt(ts, b, len);
}

/* ------------------------------------------------------- alert handling */

static void handle_alert(base_t *b, const report_t *rep, int src, double t_recv, int blocking)
{
    const config_t *cfg = b->cfg;
    double now = t_recv - b->t0;
    int r = rep->cart_rank;
    b->alerts++;
    b->msg_recv++;
    b->bytes_recv += REPORT_BYTES;

    if (r < 0 || r >= b->N || rep->n_nbrs < 0 || rep->n_nbrs > MAX_NBR) {
        b->unexpected++;
        fprintf(b->log, "invalid report from world rank %d ignored\n", src);
        return;
    }

    /* Did every neighbour report congestion within the window W? */
    int all_reported = 1;
    for (int j = 0; j < rep->n_nbrs; j++) {
        int nb = rep->nbr_rank[j];
        if (nb < 0 || nb >= b->N || !(now - b->last_alert[nb] <= cfg->window_s)) all_reported = 0;
    }

    /* the report itself proves the sender and its neighbours are busy now */
    b->last_alert[r] = now;
    b->busy_until[r] = now + cfg->window_s;
    for (int j = 0; j < rep->n_nbrs; j++) {
        int nb = rep->nbr_rank[j];
        if (nb >= 0 && nb < b->N && b->busy_until[nb] < now + cfg->window_s)
            b->busy_until[nb] = now + cfg->window_s;
    }

    int cand[MAX_CAND], nc = 0, dist = 0, status;
    if (all_reported) {
        status = RD_NONE_NEARBY;
    } else {
        nc = find_candidates(cfg->m, cfg->n, rep->row, rep->col, b->busy_until, now,
                             cand, MAX_CAND, &dist);
        status = nc > 0 ? RD_CANDIDATES : RD_NONE_IN_GRID;
    }

    int *out = b->rbuf[src];
    if (!blocking) MPI_Wait(&b->rreq[src], MPI_STATUS_IGNORE);   /* previous reply to this node */
    out[0] = status; out[1] = dist; out[2] = nc; out[3] = rep->alert_seq;
    for (int i = 0; i < nc; i++) {
        out[RD_HDR + 3 * i]     = cand[i];
        out[RD_HDR + 3 * i + 1] = cand[i] / cfg->n;
        out[RD_HDR + 3 * i + 2] = cand[i] % cfg->n;
    }
    int cnt = RD_HDR + 3 * nc;
    double svc = MPI_Wtime() - t_recv;
    out[4] = (int)(svc * 1e6 + 0.5);
    if (blocking) MPI_Send(out, cnt, MPI_INT, src, TAG_REDIRECT, MPI_COMM_WORLD);
    else MPI_Isend(out, cnt, MPI_INT, src, TAG_REDIRECT, MPI_COMM_WORLD, &b->rreq[src]);
    b->msg_sent++;
    b->bytes_sent += cnt * (long)sizeof(int);
    b->svc_sum += svc;
    if (svc > b->svc_max) b->svc_max = svc;
    dv_push(&b->s_svc, svc);
    if (status == RD_CANDIDATES) b->r_cand++;
    else if (status == RD_NONE_NEARBY) b->r_none_near++;
    else b->r_none_grid++;

    /* logging happens after the reply is posted: file I/O stays off the critical path */
    char ts[32], nbrs[256] = "", cands[1024] = "";
    stamp_fmt(rep->ts, ts, sizeof ts);
    int off = 0;
    for (int j = 0; j < rep->n_nbrs; j++)
        off += snprintf(nbrs + off, sizeof nbrs - off, "%s%d(%d,%d)=%d/%d", j ? " " : "",
                        rep->nbr_rank[j], rep->nbr_row[j], rep->nbr_col[j],
                        rep->nbr_free[j], rep->nbr_total[j]);
    off = 0;
    for (int i = 0; i < nc && off < (int)sizeof cands - 32; i++)
        off += snprintf(cands + off, sizeof cands - off, "%s%d(%d,%d)", i ? " " : "",
                        cand[i], cand[i] / cfg->n, cand[i] % cfg->n);
    char dstr[32] = "";
    if (nc) snprintf(dstr, sizeof dstr, " (distance %d)", dist);
    const char *verdict = status == RD_CANDIDATES ? "suggest nearest available"
                        : status == RD_NONE_NEARBY ? "no available node nearby"
                        : "no available node in grid";

    fprintf(b->log,
            "[t=%9.3fs] ALERT #%d from node %d (%d,%d) world %d | detected %s | free %d/%d (threshold %d)\n"
            "             neighbours: %s\n"
            "             neighbour exchange: %d msgs, %d bytes, %.3f ms | all neighbours reported within W=%gs: %s\n"
            "             -> %s%s%s%s | base service %.1f us\n",
            now, rep->alert_seq, r, rep->row, rep->col, src, ts, rep->free_ports, rep->total_ports,
            rep->threshold, rep->n_nbrs ? nbrs : "(none)", rep->msgs_exchanged, rep->bytes_exchanged,
            rep->nbr_rtt * 1e3, cfg->window_s, all_reported ? "yes" : "no", verdict,
            nc ? ": " : "", cands, dstr, svc * 1e6);
    fprintf(b->acsv, "%.6f,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,\"%s\",%d,%d,%.6f,%d,%d,%d,%d,\"%s\",%.3f\n",
            now, ts, rep->alert_seq, r, src, rep->row, rep->col, rep->free_ports, rep->total_ports,
            rep->threshold, rep->n_nbrs, nbrs, rep->msgs_exchanged, rep->bytes_exchanged,
            rep->nbr_rtt * 1e3, all_reported, status, dist, nc, cands, svc * 1e6);
    if (!cfg->quiet)
        printf("[base t=%7.3fs] alert #%d node %d (%d,%d) free %d/%d -> %s%s%s\n", now, rep->alert_seq,
               r, rep->row, rep->col, rep->free_ports, rep->total_ports, verdict, nc ? ": " : "", cands);
}

static void handle_final(base_t *b, int src)
{
    double st[ST_COUNT];
    MPI_Status s;
    int cnt;
    MPI_Recv(st, ST_COUNT, MPI_DOUBLE, src, TAG_FINAL, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    memcpy(b->nst + (size_t)src * ST_COUNT, st, sizeof st);
    MPI_Probe(src, TAG_SAMPLES, MPI_COMM_WORLD, &s);
    MPI_Get_count(&s, MPI_DOUBLE, &cnt);
    double *buf = malloc((size_t)(cnt > 0 ? cnt : 1) * sizeof(double));
    if (!buf) MPI_Abort(MPI_COMM_WORLD, 2);
    MPI_Recv(buf, cnt, MPI_DOUBLE, src, TAG_SAMPLES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    for (int i = 0; i + 2 < cnt; i += 3) {
        int type = (int)buf[i];
        fprintf(b->scsv, "%d,%d,%.9f,%.9f\n", src, type, buf[i + 1], buf[i + 2]);
        if (type == SMP_NBR) dv_push(&b->s_nbr, buf[i + 1]);
        else if (type == SMP_ALERT) { dv_push(&b->s_alert, buf[i + 1]); dv_push(&b->s_net, buf[i + 2]); }
        else if (type == SMP_SYNC) dv_push(&b->s_sync, buf[i + 1]);
    }
    free(buf);
    b->finals++;
    b->t_last_final = MPI_Wtime();
}

static void check_term(base_t *b, int blocking)
{
    if (b->term_sent || MPI_Wtime() - b->t0 < b->cfg->sim_s) return;
    b->t_term = MPI_Wtime();
    for (int i = 0; i < b->N; i++) {
        if (blocking) MPI_Send(&b->term_val, 1, MPI_INT, i, TAG_TERM, MPI_COMM_WORLD);
        else MPI_Isend(&b->term_val, 1, MPI_INT, i, TAG_TERM, MPI_COMM_WORLD, &b->treq[i]);
    }
    b->term_sent = 1;
    b->msg_sent += b->N;
    b->bytes_sent += b->N * (long)sizeof(int);
    fprintf(b->log, "[t=%9.3fs] simulation time reached: TERMINATE sent to %d nodes\n",
            b->t_term - b->t0, b->N);
}

/* -------------------------------------------------------- receive loops */

static void discard_unexpected(base_t *b, MPI_Status *s)
{
    int cnt;
    MPI_Get_count(s, MPI_BYTE, &cnt);
    char *tmp = malloc((size_t)(cnt > 0 ? cnt : 1));
    if (!tmp) MPI_Abort(MPI_COMM_WORLD, 2);
    MPI_Recv(tmp, cnt, MPI_BYTE, s->MPI_SOURCE, s->MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    free(tmp);
    b->unexpected++;
    fprintf(b->log, "unexpected tag %d from %d discarded\n", s->MPI_TAG, s->MPI_SOURCE);
}

static void loop_b1(base_t *b)
{
    report_t rep;
    int idle = 0;
    while (b->finals < b->N) {
        check_term(b, 1);
        int flag;
        MPI_Status s;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &s);
        idle_wait(b->cfg, flag, &idle);
        if (!flag) continue;
        if (s.MPI_TAG == TAG_ALERT) {
            MPI_Recv(&rep, 1, b->rtype, s.MPI_SOURCE, TAG_ALERT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            handle_alert(b, &rep, s.MPI_SOURCE, MPI_Wtime(), 1);
        } else if (s.MPI_TAG == TAG_FINAL) {
            handle_final(b, s.MPI_SOURCE);
        } else {
            discard_unexpected(b, &s);
        }
    }
}

/* FINAL messages arrive only once per node, at shutdown: probe for them */
static int poll_finals(base_t *b)
{
    int flag, got = 0;
    MPI_Status s;
    for (;;) {
        MPI_Iprobe(MPI_ANY_SOURCE, TAG_FINAL, MPI_COMM_WORLD, &flag, &s);
        if (!flag) return got;
        handle_final(b, s.MPI_SOURCE);
        got++;
    }
}

static void loop_b2(base_t *b)
{
    report_t buf, rep;
    MPI_Request req;
    MPI_Status s;
    int idle = 0;
    MPI_Irecv(&buf, 1, b->rtype, MPI_ANY_SOURCE, TAG_ALERT, MPI_COMM_WORLD, &req);
    while (b->finals < b->N) {
        check_term(b, 0);
        int flag, progress = 0;
        MPI_Test(&req, &flag, &s);
        if (flag) {
            double t = MPI_Wtime();
            rep = buf;   /* copy out, then re-post at once so the next report can land */
            MPI_Irecv(&buf, 1, b->rtype, MPI_ANY_SOURCE, TAG_ALERT, MPI_COMM_WORLD, &req);
            handle_alert(b, &rep, s.MPI_SOURCE, t, 0);
            progress = 1;
        }
        if (b->term_sent) progress += poll_finals(b);   /* FINALs only follow TERMINATE */
        idle_wait(b->cfg, progress, &idle);
    }
    MPI_Cancel(&req);
    MPI_Wait(&req, MPI_STATUS_IGNORE);
}

static void loop_b3(base_t *b)
{
    int R = b->N < B3_SLOTS ? b->N : B3_SLOTS;
    report_t slots[B3_SLOTS];
    MPI_Request req[B3_SLOTS];
    MPI_Status sts[B3_SLOTS];
    int idx[B3_SLOTS], idle = 0;
    for (int i = 0; i < R; i++)
        MPI_Recv_init(&slots[i], 1, b->rtype, MPI_ANY_SOURCE, TAG_ALERT, MPI_COMM_WORLD, &req[i]);
    MPI_Startall(R, req);
    while (b->finals < b->N) {
        check_term(b, 0);
        int outc, progress = 0;
        MPI_Testsome(R, req, &outc, idx, sts);
        if (outc != MPI_UNDEFINED && outc > 0) {
            for (int i = 0; i < outc; i++) {
                int j = idx[i];
                report_t rep = slots[j];
                MPI_Start(&req[j]);        /* re-arm the slot before the work */
                handle_alert(b, &rep, sts[i].MPI_SOURCE, MPI_Wtime(), 0);
            }
            progress = 1;
        }
        if (b->term_sent) progress += poll_finals(b);   /* FINALs only follow TERMINATE */
        idle_wait(b->cfg, progress, &idle);
    }
    for (int i = 0; i < R; i++) {
        MPI_Cancel(&req[i]);
        MPI_Wait(&req[i], MPI_STATUS_IGNORE);
        MPI_Request_free(&req[i]);
    }
}

/* ---------------------------------------------------------------- driver */

static FILE *open_out(const config_t *cfg, const char *name, const char *mode)
{
    char path[600];
    snprintf(path, sizeof path, "%s/%s", cfg->outdir, name);
    FILE *f = fopen(path, mode);
    if (!f) {
        fprintf(stderr, "base: cannot open %s, logging to /dev/null\n", path);
        f = fopen("/dev/null", "w");
    }
    if (f) setvbuf(f, NULL, _IOFBF, 1 << 16);
    return f;
}

static void write_summary(base_t *b, double sim_time, double cpu_s)
{
    const config_t *cfg = b->cfg;
    double tot[ST_COUNT] = {0}, mx_nbr = 0, mx_alert = 0, min_free = cfg->k;
    for (int i = 0; i < b->N; i++) {
        const double *st = b->nst + (size_t)i * ST_COUNT;
        for (int j = 0; j < ST_COUNT; j++) tot[j] += st[j];
        if (st[ST_NBR_RTT_MAX] > mx_nbr) mx_nbr = st[ST_NBR_RTT_MAX];
        if (st[ST_ALERT_RTT_MAX] > mx_alert) mx_alert = st[ST_ALERT_RTT_MAX];
        if (st[ST_MIN_FREE] < min_free) min_free = st[ST_MIN_FREE];
    }
    long total_msgs = (long)tot[ST_MSG_SENT] + b->msg_sent;
    long total_bytes = (long)tot[ST_BYTES_SENT] + b->bytes_sent;
    double shutdown = b->term_sent ? b->t_last_final - b->t_term : NAN;
    double nbr_mean = mean(&b->s_nbr), alert_mean = mean(&b->s_alert), net_mean = mean(&b->s_net);
    double svc_mean = mean(&b->s_svc);
    double nbr_med = pct(&b->s_nbr, 0.5), nbr_p95 = pct(&b->s_nbr, 0.95);
    double al_med = pct(&b->s_alert, 0.5), al_p95 = pct(&b->s_alert, 0.95);
    double net_med = pct(&b->s_net, 0.5), net_p95 = pct(&b->s_net, 0.95);
    double sync_med = pct(&b->s_sync, 0.5);
    char host[256] = "?";
    if (gethostname(host, sizeof host) != 0) strcpy(host, "?");
    host[sizeof host - 1] = '\0';

    FILE *L = b->log;
    fprintf(L, "\n================ SIMULATION SUMMARY ================\n");
    fprintf(L, "grid %d x %d (%d nodes) + 1 base | k=%d ports | cycle %d ms | threshold %d | W %.2f s\n",
            cfg->m, cfg->n, b->N, cfg->k, cfg->cycle_ms, cfg->threshold, cfg->window_s);
    fprintf(L, "versions N%d / B%d | poll %d us | recheck %d | seed %u\n", cfg->node_ver, cfg->base_ver,
            cfg->poll_us, cfg->recheck, cfg->seed);
    fprintf(L, "simulation time               %.3f s\n", sim_time);
    fprintf(L, "alerts received               %ld\n", b->alerts);
    fprintf(L, "  replies with candidates     %ld\n", b->r_cand);
    fprintf(L, "  replies 'none nearby'       %ld\n", b->r_none_near);
    fprintf(L, "  replies 'none in grid'      %ld\n", b->r_none_grid);
    fprintf(L, "local redirects (neighbours)  %.0f\n", tot[ST_LOCAL_REDIRECT]);
    fprintf(L, "neighbour exchanges           %.0f\n", tot[ST_EXCHANGES]);
    fprintf(L, "messages: nodes sent %.0f (req %.0f, resp %.0f, alert %.0f), base sent %ld (redirect+terminate)\n",
            tot[ST_MSG_SENT], tot[ST_REQ_SENT], tot[ST_RESP_SENT], tot[ST_ALERTS], b->msg_sent);
    fprintf(L, "total messages / bytes        %ld / %ld\n", total_msgs, total_bytes);
    fprintf(L, "neighbour exchange time       n=%d mean %.3f ms median %.3f ms p95 %.3f ms max %.3f ms\n",
            b->s_nbr.n, nbr_mean * 1e3, nbr_med * 1e3, nbr_p95 * 1e3, mx_nbr * 1e3);
    fprintf(L, "alert->redirect round trip    n=%d mean %.3f ms median %.3f ms p95 %.3f ms max %.3f ms\n",
            b->s_alert.n, alert_mean * 1e3, al_med * 1e3, al_p95 * 1e3, mx_alert * 1e3);
    fprintf(L, "  minus base service time     mean %.3f ms median %.3f ms p95 %.3f ms\n",
            net_mean * 1e3, net_med * 1e3, net_p95 * 1e3);
    fprintf(L, "base service time             mean %.1f us max %.1f us\n", svc_mean * 1e6, b->svc_max * 1e6);
    if (b->s_sync.n)
        fprintf(L, "termination Allreduce (N1/N3) n=%d median %.3f ms\n", b->s_sync.n, sync_med * 1e3);
    fprintf(L, "shutdown (TERMINATE -> last FINAL) %.3f ms\n", shutdown * 1e3);
    fprintf(L, "base CPU time                 %.3f s\n", cpu_s);
    fprintf(L, "unexpected/invalid messages   %ld\n", b->unexpected);
    fprintf(L, "lowest availability observed  %.0f / %d\n", min_free, cfg->k);
    fprintf(L, "\nper-node statistics\n%-5s", "node");
    for (int j = 0; j < ST_COUNT; j++) fprintf(L, " %s", STAT_NAMES[j]);
    fputc('\n', L);
    for (int i = 0; i < b->N; i++) {
        fprintf(L, "%-5d", i);
        for (int j = 0; j < ST_COUNT; j++) fprintf(L, " %.6g", b->nst[(size_t)i * ST_COUNT + j]);
        fputc('\n', L);
    }

    FILE *S = open_out(cfg, "summary.csv", "w");
    if (S) {
        fprintf(S, "m,n,nodes,k,cycle_ms,sim_s,threshold,window_s,node_ver,base_ver,poll_us,recheck,seed,p_arrive,"
                   "sim_time,alerts,cand,none_nearby,none_grid,local_redirects,exchanges,total_msgs,total_bytes,"
                   "nbr_n,nbr_mean_ms,nbr_med_ms,nbr_p95_ms,alert_n,alert_mean_ms,alert_med_ms,alert_p95_ms,"
                   "net_mean_ms,net_med_ms,net_p95_ms,svc_mean_us,svc_max_us,sync_med_ms,shutdown_ms,base_cpu_s,"
                   "unexpected,host\n");
        fprintf(S, "%d,%d,%d,%d,%d,%.3f,%d,%.3f,%d,%d,%d,%d,%u,%.3f,%.6f,%ld,%ld,%ld,%ld,%.0f,%.0f,%ld,%ld,"
                   "%d,%.6f,%.6f,%.6f,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.3f,%.3f,%.6f,%.6f,%.6f,%ld,%s\n",
                cfg->m, cfg->n, b->N, cfg->k, cfg->cycle_ms, cfg->sim_s, cfg->threshold, cfg->window_s,
                cfg->node_ver, cfg->base_ver, cfg->poll_us, cfg->recheck, cfg->seed, cfg->p_arrive, sim_time, b->alerts,
                b->r_cand, b->r_none_near, b->r_none_grid, tot[ST_LOCAL_REDIRECT], tot[ST_EXCHANGES],
                total_msgs, total_bytes, b->s_nbr.n, nbr_mean * 1e3, nbr_med * 1e3, nbr_p95 * 1e3,
                b->s_alert.n, alert_mean * 1e3, al_med * 1e3, al_p95 * 1e3, net_mean * 1e3, net_med * 1e3,
                net_p95 * 1e3, svc_mean * 1e6, b->svc_max * 1e6, sync_med * 1e3, shutdown * 1e3, cpu_s,
                b->unexpected, host);
        fclose(S);
    }

    if (!cfg->quiet) {
        printf("\n[base] %dx%d N%d/B%d: %.2f s, %ld alerts (%ld candidates, %ld none nearby, %ld none in grid), "
               "%.0f local redirects, %ld messages\n",
               cfg->m, cfg->n, cfg->node_ver, cfg->base_ver, sim_time, b->alerts, b->r_cand, b->r_none_near,
               b->r_none_grid, tot[ST_LOCAL_REDIRECT], total_msgs);
        printf("[base] neighbour exchange median %.3f ms | alert round trip median %.3f ms | shutdown %.3f ms\n",
               nbr_med * 1e3, al_med * 1e3, shutdown * 1e3);
        printf("[base] logs written to %s/\n", cfg->outdir);
    }
}

int run_base(const config_t *cfg, int n_nodes, MPI_Datatype report_type)
{
    base_t *b = calloc(1, sizeof *b);
    if (!b) MPI_Abort(MPI_COMM_WORLD, 2);
    b->cfg = cfg;
    b->N = n_nodes;
    b->rtype = report_type;
    b->last_alert = malloc((size_t)n_nodes * sizeof(double));
    b->busy_until = malloc((size_t)n_nodes * sizeof(double));
    b->nst = calloc((size_t)n_nodes * ST_COUNT, sizeof(double));
    b->rbuf = malloc((size_t)n_nodes * sizeof(*b->rbuf));
    b->rreq = malloc((size_t)n_nodes * sizeof(MPI_Request));
    b->treq = malloc((size_t)n_nodes * sizeof(MPI_Request));
    if (!b->last_alert || !b->busy_until || !b->nst || !b->rbuf || !b->rreq || !b->treq)
        MPI_Abort(MPI_COMM_WORLD, 2);
    for (int i = 0; i < n_nodes; i++) {
        b->last_alert[i] = -1e30;
        b->busy_until[i] = -1e30;
        b->rreq[i] = b->treq[i] = MPI_REQUEST_NULL;
    }

    if (mkdir_p(cfg->outdir) != 0) fprintf(stderr, "base: cannot create %s\n", cfg->outdir);
    b->log = open_out(cfg, "base_log.txt", "w");
    b->acsv = open_out(cfg, "alerts.csv", "w");
    b->scsv = open_out(cfg, "samples.csv", "w");
    if (!b->log || !b->acsv || !b->scsv) MPI_Abort(MPI_COMM_WORLD, 3);
    fprintf(b->acsv, "t_s,detected,seq,cart_rank,world_rank,row,col,free,total,threshold,n_nbrs,neighbours,"
                     "nbr_msgs,nbr_bytes,nbr_rtt_ms,all_reported,status,distance,n_cand,candidates,service_us\n");
    fprintf(b->scsv, "node,type,v1,v2\n");

    char start[32];
    now_str(start, sizeof start);
    fprintf(b->log, "EV charging network simulation - base station log\n");
    fprintf(b->log, "start %s | grid %d x %d (%d nodes) | k=%d | cycle %d ms | sim %.1f s | threshold %d | "
                    "W %.2f s | N%d/B%d | poll %d us | recheck %d | seed %u | p_arrive %.2f p_depart %.2f | hot [%s] p %.2f/%.2f | "
                    "cooldown %d\n\n",
            start, cfg->m, cfg->n, n_nodes, cfg->k, cfg->cycle_ms, cfg->sim_s, cfg->threshold, cfg->window_s,
            cfg->node_ver, cfg->base_ver, cfg->poll_us, cfg->recheck, cfg->seed, cfg->p_arrive, cfg->p_depart, cfg->hot,
            cfg->p_arrive_hot, cfg->p_depart_hot, cfg->cooldown);

    MPI_Barrier(MPI_COMM_WORLD);      /* common start */
    b->t0 = MPI_Wtime();

    if (cfg->base_ver == 1) loop_b1(b);
    else if (cfg->base_ver == 2) loop_b2(b);
    else loop_b3(b);

    MPI_Waitall(n_nodes, b->rreq, MPI_STATUSES_IGNORE);
    MPI_Waitall(n_nodes, b->treq, MPI_STATUSES_IGNORE);
    double sim_time = MPI_Wtime() - b->t0;

    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    double cpu = ru.ru_utime.tv_sec + ru.ru_utime.tv_usec * 1e-6 + ru.ru_stime.tv_sec + ru.ru_stime.tv_usec * 1e-6;
    write_summary(b, sim_time, cpu);

    char end[32];
    now_str(end, sizeof end);
    fprintf(b->log, "\nend %s - all %d nodes confirmed shutdown\n", end, b->finals);
    fclose(b->log);
    fclose(b->acsv);
    fclose(b->scsv);
    free(b->last_alert); free(b->busy_until); free(b->nst);
    free(b->rbuf); free(b->rreq); free(b->treq);
    free(b->s_nbr.v); free(b->s_alert.v); free(b->s_net.v); free(b->s_sync.v); free(b->s_svc.v);
    free(b);
    return 0;
}
