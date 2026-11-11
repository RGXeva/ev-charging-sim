/* SPDX-License-Identifier: MIT */
/*
 * EV charging node: one MPI process per node inside a 2-D Cartesian communicator.
 *
 * Threads (OpenMP parallel region with k+1 threads):
 *   thread 0      controller - the only thread that calls MPI (MPI_THREAD_FUNNELED)
 *   threads 1..k  charging ports - each periodically updates its own slot of the
 *                 shared array ports.in_use[] (atomic read/write, no locks)
 *
 * Versions of the controller:
 *   N1  blocking point-to-point : every cycle MPI_Sendrecv with the 4 neighbours,
 *       blocking MPI_Send/MPI_Recv to the base, MPI_Allreduce to agree on termination
 *   N2  non-blocking on demand  : neighbours are only contacted when the node is
 *       congested; persistent receives, MPI_Isend/MPI_Irecv/MPI_Test* event loop
 *   N3  neighbourhood collective: every cycle MPI_Neighbor_allgather on comm2D
 */
#define _POSIX_C_SOURCE 200809L
#include "common.h"

#include <omp.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

enum { S_IDLE, S_WAIT_NBR, S_WAIT_BASE };

#define REQ_BYTES    (2 * (int)sizeof(int))
#define INFO_BYTES   (NBR_INFO_INTS * (int)sizeof(int))
#define REPORT_BYTES (REPORT_DOUBLES * (int)sizeof(double) + REPORT_INTS * (int)sizeof(int))

typedef struct {
    int ts[7];
    int cycle;
    int free_ports;
} hist_t;

typedef struct {
    int in_use[MAX_PORTS];          /* shared array: 1 = port occupied        */
    double last_change[MAX_PORTS];  /* written only by the owning port thread */
    long changes[MAX_PORTS];        /* written only by the owning port thread */
    int stop;                       /* set by the controller at shutdown      */
} ports_t;

typedef struct {
    const config_t *cfg;
    MPI_Comm comm;
    int base;                       /* world rank of the base station */
    MPI_Datatype rtype;
    int rank, world_rank, row, col, k, thr;
    int dir[MAX_NBR];               /* up, down, left, right (may be MPI_PROC_NULL) */
    int nv, nbr[MAX_NBR];           /* valid neighbours, same order */
    ports_t ports;
    double t0;                      /* common start (monotonic clock) */
    int cycle, alert_seq, cooldown_until;
    double st[ST_COUNT];
    double *smp; int nsmp, smp_cap;
    hist_t *hist; int nhist, hist_cap;
    FILE *log;
    /* N2 state */
    MPI_Request req_in[MAX_NBR], resp_out[MAX_NBR], req_out[MAX_NBR], resp_in[MAX_NBR];
    int req_in_buf[MAX_NBR][2], req_out_buf[2];
    nbr_info_t resp_out_buf[MAX_NBR], resp_in_buf[MAX_NBR];
    MPI_Request alert_req, base_req;
    report_t rep;
    int base_buf[RD_MAX_INTS];
} node_t;

static void nlog(node_t *nd, const char *fmt, ...)
{
    if (!nd->log) return;
    int ts[7];
    char b[32];
    stamp_now(ts, NULL);
    stamp_fmt(ts, b, sizeof b);
    fprintf(nd->log, "[%s] cycle %-4d ", b, nd->cycle);
    va_list ap;
    va_start(ap, fmt);
    vfprintf(nd->log, fmt, ap);
    va_end(ap);
    fputc('\n', nd->log);
}

static void add_sample(node_t *nd, int type, double v1, double v2)
{
    if (nd->nsmp == nd->smp_cap) {
        int cap = nd->smp_cap ? 2 * nd->smp_cap : 256;
        double *p = realloc(nd->smp, (size_t)cap * 3 * sizeof(double));
        if (!p) return;
        nd->smp = p;
        nd->smp_cap = cap;
    }
    double *s = nd->smp + 3 * nd->nsmp++;
    s[0] = type; s[1] = v1; s[2] = v2;
}

/* ------------------------------------------------------------------ ports */

static void port_loop(node_t *nd, int port)
{
    const config_t *cfg = nd->cfg;
    unsigned s = cfg->seed * 2654435761u + (unsigned)nd->rank * 7919u + (unsigned)port * 104729u + 1u;
    int hot = is_hot(cfg, nd->rank);
    double pa = hot ? cfg->p_arrive_hot : cfg->p_arrive;
    double pd = hot ? cfg->p_depart_hot : cfg->p_depart;
    double cyc = cfg->cycle_ms / 1000.0;

    /* start from the steady-state occupancy pa / (pa + pd) */
    double r0 = rand_r(&s) / ((double)RAND_MAX + 1.0);
    #pragma omp atomic write
    nd->ports.in_use[port] = (pa + pd > 0) && r0 < pa / (pa + pd);

    /* ports update mid-cycle so they never race the controller's sample */
    double next = nd->t0 + 0.5 * cyc;
    for (;;) {
        int stop;
        #pragma omp atomic read
        stop = nd->ports.stop;
        if (stop) break;
        double now = mono_now();
        if (now < next) {
            double w = next - now;
            sleep_us(w > 0.02 ? 20000 : (long)(w * 1e6) + 1);
            continue;
        }
        next += cyc;
        double r = rand_r(&s) / ((double)RAND_MAX + 1.0);
        int cur;
        #pragma omp atomic read
        cur = nd->ports.in_use[port];
        int newv = cur;
        if (!cur && r < pa) newv = 1;
        else if (cur && r < pd) newv = 0;
        if (newv != cur) {
            #pragma omp atomic write
            nd->ports.in_use[port] = newv;
            nd->ports.last_change[port] = now;
            nd->ports.changes[port]++;
        }
    }
}

static int cur_free(node_t *nd)
{
    int used = 0;
    for (int i = 0; i < nd->k; i++) {
        int v;
        #pragma omp atomic read
        v = nd->ports.in_use[i];
        used += v;
    }
    return nd->k - used;
}

/* start of a cycle: read the shared array and keep a timestamped history */
static int sample_cycle(node_t *nd)
{
    nd->cycle++;
    int f = cur_free(nd);
    if (nd->nhist == nd->hist_cap) {
        int cap = nd->hist_cap ? 2 * nd->hist_cap : 128;
        hist_t *p = realloc(nd->hist, (size_t)cap * sizeof(hist_t));
        if (p) { nd->hist = p; nd->hist_cap = cap; }
    }
    if (nd->nhist < nd->hist_cap) {
        hist_t *h = &nd->hist[nd->nhist++];
        stamp_now(h->ts, NULL);
        h->cycle = nd->cycle;
        h->free_ports = f;
    }
    nd->st[ST_CYCLES]++;
    nd->st[ST_SUM_FREE] += f;
    if (f < nd->st[ST_MIN_FREE]) nd->st[ST_MIN_FREE] = f;
    return f;
}

static void fill_info(node_t *nd, nbr_info_t *in, int free_ports)
{
    in->rank = nd->rank; in->row = nd->row; in->col = nd->col;
    in->free_ports = free_ports; in->total_ports = nd->k; in->cycle = nd->cycle;
}

/* ------------------------------------------------- shared decision logic */

/* Decide after a neighbour exchange. Returns 1 if the base must be alerted;
 * the report is then ready in nd->rep. */
static int evaluate(node_t *nd, int free_now, const nbr_info_t *info, int ninfo,
                    int msgs, int bytes, double rtt)
{
    int best = -1;
    for (int j = 0; j < ninfo; j++)
        if (info[j].free_ports >= nd->thr &&
            (best < 0 || info[j].free_ports > info[best].free_ports))
            best = j;

    if (best >= 0) {
        nd->st[ST_LOCAL_REDIRECT]++;
        nlog(nd, "neighbour exchange %.3f ms: redirect incoming EVs to node %d (%d,%d), %d/%d free",
             rtt * 1e3, info[best].rank, info[best].row, info[best].col,
             info[best].free_ports, info[best].total_ports);
        nd->cooldown_until = nd->cycle + nd->cfg->cooldown;
        return 0;
    }

    report_t *r = &nd->rep;
    memset(r, 0, sizeof *r);
    stamp_now(r->ts, &r->t_wall);
    r->nbr_rtt = rtt;
    r->cycle = nd->cycle;
    r->cart_rank = nd->rank; r->world_rank = nd->world_rank;
    r->row = nd->row; r->col = nd->col;
    r->total_ports = nd->k; r->free_ports = free_now; r->threshold = nd->thr;
    r->n_nbrs = ninfo;
    for (int j = 0; j < MAX_NBR; j++) {
        int ok = j < ninfo;
        r->nbr_rank[j]  = ok ? info[j].rank : -1;
        r->nbr_row[j]   = ok ? info[j].row : -1;
        r->nbr_col[j]   = ok ? info[j].col : -1;
        r->nbr_free[j]  = ok ? info[j].free_ports : -1;
        r->nbr_total[j] = ok ? info[j].total_ports : -1;
    }
    r->msgs_exchanged = msgs;
    r->bytes_exchanged = bytes;
    r->alert_seq = ++nd->alert_seq;
    nlog(nd, "all %d neighbours below threshold (exchange %.3f ms) -> ALERT #%d to base",
         ninfo, rtt * 1e3, r->alert_seq);
    return 1;
}

static void handle_redirect(node_t *nd, const int *buf, int count, double rtt)
{
    int status = buf[0], dist = buf[1], nc = buf[2], seq = buf[3], svc = buf[4];
    if (count < RD_HDR || nc < 0 || count != RD_HDR + 3 * nc) {
        nlog(nd, "malformed redirect (count %d)", count);
        nc = 0;
    }
    double net = rtt - svc * 1e-6;
    add_sample(nd, SMP_ALERT, rtt, net);
    nd->st[ST_ALERT_RTT_SUM] += rtt;
    nd->st[ST_ALERT_NET_SUM] += net;
    if (rtt > nd->st[ST_ALERT_RTT_MAX]) nd->st[ST_ALERT_RTT_MAX] = rtt;
    nd->st[ST_MSG_RECV]++;
    nd->st[ST_BYTES_RECV] += count * (int)sizeof(int);

    if (status == RD_CANDIDATES && nc > 0) {
        nd->st[ST_BASE_CAND]++;
        if (nd->log) {
            char line[1024];
            int off = 0;
            for (int i = 0; i < nc && off < (int)sizeof line - 32; i++)
                off += snprintf(line + off, sizeof line - off, " %d(%d,%d)",
                                buf[RD_HDR + 3 * i], buf[RD_HDR + 3 * i + 1], buf[RD_HDR + 3 * i + 2]);
            nlog(nd, "base reply #%d (rtt %.3f ms, base service %d us): %d nearest node(s) at distance %d:%s -> forwarded to EVs",
                 seq, rtt * 1e3, svc, nc, dist, line);
        }
    } else {
        nd->st[ST_BASE_NONE]++;
        nlog(nd, "base reply #%d (rtt %.3f ms): %s", seq, rtt * 1e3,
             status == RD_NONE_NEARBY ? "no available node nearby" : "no available node in the grid");
    }
    nd->cooldown_until = nd->cycle + nd->cfg->cooldown;
}

/* ------------------------------------------------ N2: non-blocking, on demand */

/* answer every neighbour request that has arrived; returns how many */
static int serve_requests(node_t *nd)
{
    if (nd->nv == 0) return 0;
    int outc, idx[MAX_NBR];
    MPI_Testsome(nd->nv, nd->req_in, &outc, idx, MPI_STATUSES_IGNORE);
    if (outc == MPI_UNDEFINED) return 0;
    for (int i = 0; i < outc; i++) {
        int j = idx[i];
        MPI_Wait(&nd->resp_out[j], MPI_STATUS_IGNORE);   /* previous answer must be done */
        fill_info(nd, &nd->resp_out_buf[j], cur_free(nd));
        MPI_Isend(&nd->resp_out_buf[j], NBR_INFO_INTS, MPI_INT, nd->nbr[j], TAG_NBR_RESP,
                  nd->comm, &nd->resp_out[j]);
        MPI_Start(&nd->req_in[j]);
        nd->st[ST_REQ_RECV]++;  nd->st[ST_MSG_RECV]++; nd->st[ST_BYTES_RECV] += REQ_BYTES;
        nd->st[ST_RESP_SENT]++; nd->st[ST_MSG_SENT]++; nd->st[ST_BYTES_SENT] += INFO_BYTES;
    }
    return outc;
}

static void controller_n2(node_t *nd)
{
    const config_t *cfg = nd->cfg;
    double cyc = cfg->cycle_ms / 1000.0, next = nd->t0 + cyc;
    int state = S_IDLE, terminating = 0, free_at_ex = 0, idle = 0;
    double t_ex = 0, t_alert = 0;

    for (int j = 0; j < nd->nv; j++) {
        MPI_Recv_init(nd->req_in_buf[j], 2, MPI_INT, nd->nbr[j], TAG_NBR_REQ, nd->comm, &nd->req_in[j]);
        nd->resp_out[j] = nd->req_out[j] = nd->resp_in[j] = MPI_REQUEST_NULL;
    }
    if (nd->nv) MPI_Startall(nd->nv, nd->req_in);
    nd->alert_req = MPI_REQUEST_NULL;
    MPI_Irecv(nd->base_buf, RD_MAX_INTS, MPI_INT, nd->base, MPI_ANY_TAG, MPI_COMM_WORLD, &nd->base_req);

    for (;;) {
        int progress = serve_requests(nd);

        /* messages from the base station: REDIRECT or TERMINATE */
        int flag;
        MPI_Status s;
        MPI_Test(&nd->base_req, &flag, &s);
        if (flag) {
            progress = 1;
            if (s.MPI_TAG == TAG_TERM) {
                terminating = 1;
                nd->st[ST_MSG_RECV]++; nd->st[ST_BYTES_RECV] += sizeof(int);
                nlog(nd, "TERMINATE received from base station");
            } else if (s.MPI_TAG == TAG_REDIRECT && state == S_WAIT_BASE) {
                int count;
                MPI_Get_count(&s, MPI_INT, &count);
                MPI_Wait(&nd->alert_req, MPI_STATUS_IGNORE);
                handle_redirect(nd, nd->base_buf, count, MPI_Wtime() - t_alert);
                state = S_IDLE;
            } else {
                nlog(nd, "unexpected message tag %d from base", s.MPI_TAG);
            }
            MPI_Irecv(nd->base_buf, RD_MAX_INTS, MPI_INT, nd->base, MPI_ANY_TAG, MPI_COMM_WORLD, &nd->base_req);
        }

        /* new cycle: sample availability, start an exchange if congested */
        double now = mono_now();
        if (now >= next) {
            next += cyc;
            if (next <= now) next = now + cyc;     /* fell behind: skip, do not burst */
            int f = sample_cycle(nd);
            if (!terminating && state == S_IDLE && f < nd->thr && nd->cycle >= nd->cooldown_until) {
                nd->st[ST_EXCHANGES]++;
                nlog(nd, "free %d/%d < threshold %d: requesting data from %d neighbour(s)",
                     f, nd->k, nd->thr, nd->nv);
                free_at_ex = f;
                t_ex = MPI_Wtime();
                MPI_Waitall(nd->nv, nd->req_out, MPI_STATUSES_IGNORE);
                nd->req_out_buf[0] = nd->cycle;
                nd->req_out_buf[1] = f;
                for (int j = 0; j < nd->nv; j++) {
                    MPI_Irecv(&nd->resp_in_buf[j], NBR_INFO_INTS, MPI_INT, nd->nbr[j], TAG_NBR_RESP,
                              nd->comm, &nd->resp_in[j]);
                    MPI_Isend(nd->req_out_buf, 2, MPI_INT, nd->nbr[j], TAG_NBR_REQ, nd->comm, &nd->req_out[j]);
                }
                nd->st[ST_REQ_SENT] += nd->nv; nd->st[ST_MSG_SENT] += nd->nv;
                nd->st[ST_BYTES_SENT] += nd->nv * REQ_BYTES;
                state = S_WAIT_NBR;
                progress = 1;
            }
        }

        if (state == S_WAIT_NBR) {
            int done = 1;
            if (nd->nv) MPI_Testall(nd->nv, nd->resp_in, &done, MPI_STATUSES_IGNORE);
            if (done) {
                double rtt = MPI_Wtime() - t_ex;
                add_sample(nd, SMP_NBR, rtt, nd->nv);
                nd->st[ST_NBR_RTT_SUM] += rtt;
                if (rtt > nd->st[ST_NBR_RTT_MAX]) nd->st[ST_NBR_RTT_MAX] = rtt;
                nd->st[ST_RESP_RECV] += nd->nv; nd->st[ST_MSG_RECV] += nd->nv;
                nd->st[ST_BYTES_RECV] += nd->nv * INFO_BYTES;
                if (terminating) {
                    nlog(nd, "neighbour exchange finished after TERMINATE: no new alert");
                    state = S_IDLE;
                } else if (evaluate(nd, free_at_ex, nd->resp_in_buf, nd->nv, 2 * nd->nv,
                                    nd->nv * (REQ_BYTES + INFO_BYTES), rtt)) {
                    t_alert = MPI_Wtime();
                    MPI_Isend(&nd->rep, 1, nd->rtype, nd->base, TAG_ALERT, MPI_COMM_WORLD, &nd->alert_req);
                    nd->st[ST_ALERTS]++; nd->st[ST_MSG_SENT]++; nd->st[ST_BYTES_SENT] += REPORT_BYTES;
                    state = S_WAIT_BASE;
                } else {
                    state = S_IDLE;
                }
                progress = 1;
            }
        }

        if (terminating && state == S_IDLE) break;
        idle_wait(cfg, progress, &idle);
    }
}

/* After TERMINATE: keep answering neighbours until every node has stopped
 * initiating requests (MPI_Ibarrier), then retire all pending requests. */
static void shutdown_n2(node_t *nd)
{
    MPI_Request bar;
    int idle = 0;
    MPI_Ibarrier(nd->comm, &bar);
    for (;;) {
        int p = serve_requests(nd), done;
        MPI_Test(&bar, &done, MPI_STATUS_IGNORE);
        if (done) break;
        idle_wait(nd->cfg, p, &idle);
    }
    MPI_Status s;
    int cancelled;
    for (int j = 0; j < nd->nv; j++) {
        MPI_Cancel(&nd->req_in[j]);
        MPI_Wait(&nd->req_in[j], &s);
        MPI_Test_cancelled(&s, &cancelled);
        if (!cancelled) nlog(nd, "warning: neighbour request arrived during shutdown");
        MPI_Request_free(&nd->req_in[j]);
    }
    MPI_Cancel(&nd->base_req);
    MPI_Wait(&nd->base_req, &s);
    MPI_Test_cancelled(&s, &cancelled);
    if (!cancelled) nlog(nd, "warning: base message arrived during shutdown");
    MPI_Waitall(nd->nv, nd->resp_out, MPI_STATUSES_IGNORE);
    MPI_Waitall(nd->nv, nd->req_out, MPI_STATUSES_IGNORE);
    MPI_Wait(&nd->alert_req, MPI_STATUS_IGNORE);
}

/* ----------------------------------- N1 / N3: periodic, blocking exchange */

static void controller_periodic(node_t *nd, int collective)
{
    const config_t *cfg = nd->cfg;
    double cyc = cfg->cycle_ms / 1000.0, next = nd->t0 + cyc;
    int term_local = 0, x;

    for (;;) {
        double now = mono_now();
        if (now < next) sleep_us((long)((next - now) * 1e6));
        next += cyc;                       /* lock-step: one exchange per iteration */
        int f = sample_cycle(nd);

        if (!term_local) {
            int flag;
            MPI_Iprobe(nd->base, TAG_TERM, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
            if (flag) {
                MPI_Recv(&x, 1, MPI_INT, nd->base, TAG_TERM, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                term_local = 1;
                nd->st[ST_MSG_RECV]++; nd->st[ST_BYTES_RECV] += sizeof(int);
                nlog(nd, "TERMINATE received from base station");
            }
        }
        /* all nodes must leave the lock-step loop in the same iteration */
        int term_all;
        double ts = MPI_Wtime();
        MPI_Allreduce(&term_local, &term_all, 1, MPI_INT, MPI_MAX, nd->comm);
        double tsync = MPI_Wtime() - ts;
        add_sample(nd, SMP_SYNC, tsync, 0);
        nd->st[ST_SYNC_SUM] += tsync;
        if (term_all) break;

        nbr_info_t me, got[MAX_NBR];
        fill_info(nd, &me, f);
        for (int i = 0; i < MAX_NBR; i++) got[i].rank = -1;
        double te = MPI_Wtime();
        if (collective) {
            /* Cartesian neighbour order: up, down, left, right */
            MPI_Neighbor_allgather(&me, NBR_INFO_INTS, MPI_INT, got, NBR_INFO_INTS, MPI_INT, nd->comm);
        } else {
            for (int d = 0; d < 2; d++) {
                int lo = nd->dir[2 * d], hi = nd->dir[2 * d + 1];
                MPI_Sendrecv(&me, NBR_INFO_INTS, MPI_INT, hi, TAG_NBR_RESP,
                             &got[2 * d], NBR_INFO_INTS, MPI_INT, lo, TAG_NBR_RESP, nd->comm, MPI_STATUS_IGNORE);
                MPI_Sendrecv(&me, NBR_INFO_INTS, MPI_INT, lo, TAG_NBR_RESP,
                             &got[2 * d + 1], NBR_INFO_INTS, MPI_INT, hi, TAG_NBR_RESP, nd->comm, MPI_STATUS_IGNORE);
            }
        }
        double rtt = MPI_Wtime() - te;
        nbr_info_t info[MAX_NBR];
        int ninfo = 0;
        for (int i = 0; i < MAX_NBR; i++)
            if (nd->dir[i] != MPI_PROC_NULL) info[ninfo++] = got[i];
        add_sample(nd, SMP_NBR, rtt, ninfo);
        nd->st[ST_EXCHANGES]++;
        nd->st[ST_NBR_RTT_SUM] += rtt;
        if (rtt > nd->st[ST_NBR_RTT_MAX]) nd->st[ST_NBR_RTT_MAX] = rtt;
        nd->st[ST_RESP_SENT] += ninfo; nd->st[ST_RESP_RECV] += ninfo;
        nd->st[ST_MSG_SENT] += ninfo;  nd->st[ST_MSG_RECV] += ninfo;
        nd->st[ST_BYTES_SENT] += ninfo * INFO_BYTES; nd->st[ST_BYTES_RECV] += ninfo * INFO_BYTES;

        if (f < nd->thr && nd->cycle >= nd->cooldown_until) {
            nlog(nd, "free %d/%d < threshold %d: checking %d neighbour record(s)", f, nd->k, nd->thr, ninfo);
            if (evaluate(nd, f, info, ninfo, 2 * ninfo, 2 * ninfo * INFO_BYTES, rtt)) {
                double t = MPI_Wtime();
                MPI_Send(&nd->rep, 1, nd->rtype, nd->base, TAG_ALERT, MPI_COMM_WORLD);
                nd->st[ST_ALERTS]++; nd->st[ST_MSG_SENT]++; nd->st[ST_BYTES_SENT] += REPORT_BYTES;
                MPI_Status s;
                int count;
                MPI_Recv(nd->base_buf, RD_MAX_INTS, MPI_INT, nd->base, TAG_REDIRECT, MPI_COMM_WORLD, &s);
                MPI_Get_count(&s, MPI_INT, &count);
                handle_redirect(nd, nd->base_buf, count, MPI_Wtime() - t);
            }
        }
    }
    if (!term_local) {   /* another node saw TERMINATE first; consume ours */
        MPI_Recv(&x, 1, MPI_INT, nd->base, TAG_TERM, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        nd->st[ST_MSG_RECV]++; nd->st[ST_BYTES_RECV] += sizeof(int);
        nlog(nd, "TERMINATE received from base station");
    }
}

/* ------------------------------------------------------------------ main */

static void write_node_log_tail(node_t *nd)
{
    if (!nd->log) return;
    fprintf(nd->log, "\nAvailability history (shared port array sampled every %d ms)\n", nd->cfg->cycle_ms);
    fprintf(nd->log, "%-6s %-5s %-3s %-4s %-6s %-6s %-6s %s\n",
            "Year", "Month", "Day", "Hour", "Minute", "Second", "Cycle", "Availability");
    for (int i = 0; i < nd->nhist; i++) {
        const hist_t *h = &nd->hist[i];
        fprintf(nd->log, "%-6d %02d    %02d  %02d   %02d     %02d.%03d %-6d %d\n",
                h->ts[0], h->ts[1], h->ts[2], h->ts[3], h->ts[4], h->ts[5], h->ts[6], h->cycle, h->free_ports);
    }
    fprintf(nd->log, "\nPort state changes:");
    for (int i = 0; i < nd->k; i++) fprintf(nd->log, " p%d=%ld", i, nd->ports.changes[i]);
    fprintf(nd->log, "\n\nStatistics\n");
    for (int i = 0; i < ST_COUNT; i++) fprintf(nd->log, "  %-16s %.6g\n", STAT_NAMES[i], nd->st[i]);
}

int run_node(const config_t *cfg, MPI_Comm comm2D, int base_world, MPI_Datatype report_type)
{
    node_t *nd = calloc(1, sizeof *nd);
    if (!nd) MPI_Abort(MPI_COMM_WORLD, 2);
    nd->cfg = cfg;
    nd->comm = comm2D;
    nd->base = base_world;
    nd->rtype = report_type;
    nd->k = cfg->k;
    nd->thr = cfg->threshold;
    int coords[2];
    MPI_Comm_rank(comm2D, &nd->rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &nd->world_rank);
    MPI_Cart_coords(comm2D, nd->rank, 2, coords);
    nd->row = coords[0];
    nd->col = coords[1];
    MPI_Cart_shift(comm2D, 0, 1, &nd->dir[0], &nd->dir[1]);   /* up, down    */
    MPI_Cart_shift(comm2D, 1, 1, &nd->dir[2], &nd->dir[3]);   /* left, right */
    for (int i = 0; i < MAX_NBR; i++)
        if (nd->dir[i] != MPI_PROC_NULL) nd->nbr[nd->nv++] = nd->dir[i];
    nd->st[ST_MIN_FREE] = nd->k;

    MPI_Barrier(MPI_COMM_WORLD);          /* base has created the output directory */
    nd->t0 = mono_now();

    if (cfg->node_logs) {
        char path[600];
        snprintf(path, sizeof path, "%s/node_%03d.log", cfg->outdir, nd->rank);
        nd->log = fopen(path, "w");
        if (nd->log) {
            setvbuf(nd->log, NULL, _IOFBF, 1 << 16);
            fprintf(nd->log, "EV charging node %d  coord (%d,%d)  world rank %d  version N%d\n",
                    nd->rank, nd->row, nd->col, nd->world_rank, cfg->node_ver);
            fprintf(nd->log, "neighbours up=%d down=%d left=%d right=%d (-ve = none)  ports k=%d  threshold %d%s\n\n",
                    nd->dir[0], nd->dir[1], nd->dir[2], nd->dir[3], nd->k, nd->thr,
                    is_hot(cfg, nd->rank) ? "  [hot node]" : "");
        }
    }

    omp_set_dynamic(0);
    #pragma omp parallel num_threads(cfg->k + 1)
    {
        int tid = omp_get_thread_num();
        if (tid == 0) {
            if (omp_get_num_threads() != cfg->k + 1)
                nlog(nd, "warning: only %d OpenMP threads", omp_get_num_threads());
            if (cfg->node_ver == 2) controller_n2(nd);
            else controller_periodic(nd, cfg->node_ver == 3);
            #pragma omp atomic write
            nd->ports.stop = 1;
        } else if (tid - 1 < cfg->k) {
            port_loop(nd, tid - 1);
        }
    }   /* implicit barrier: every port thread has stopped here */

    if (cfg->node_ver == 2) shutdown_n2(nd);

    for (int i = 0; i < nd->k; i++) nd->st[ST_PORT_CHANGES] += nd->ports.changes[i];
    MPI_Send(nd->st, ST_COUNT, MPI_DOUBLE, nd->base, TAG_FINAL, MPI_COMM_WORLD);
    MPI_Send(nd->smp, 3 * nd->nsmp, MPI_DOUBLE, nd->base, TAG_SAMPLES, MPI_COMM_WORLD);

    nlog(nd, "shutdown complete");
    write_node_log_tail(nd);
    if (nd->log) fclose(nd->log);
    free(nd->smp);
    free(nd->hist);
    free(nd);
    return 0;
}
