/* SPDX-License-Identifier: MIT */
/* Utilities: argument parsing, timestamps, derived datatype, small helpers. */
#define _POSIX_C_SOURCE 200809L
#include "common.h"

#include <errno.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

const char *STAT_NAMES[ST_COUNT] = {
    "cycles", "alerts", "local_redirects", "base_candidates", "base_none",
    "exchanges", "req_sent", "req_recv", "resp_sent", "resp_recv",
    "msg_sent", "msg_recv", "bytes_sent", "bytes_recv",
    "nbr_rtt_sum", "nbr_rtt_max", "alert_rtt_sum", "alert_net_sum",
    "alert_rtt_max", "sync_sum", "min_free", "sum_free", "port_changes"
};

void print_usage(FILE *f)
{
    fprintf(f,
        "usage: mpirun -np <m*n+1> wsn [m n] [options]\n"
        "  m n          grid rows and columns (default: MPI_Dims_create on size-1)\n"
        "  -k ports     charging ports per node            (default 6)\n"
        "  -c ms        port update cycle in milliseconds  (default 1000)\n"
        "  -t seconds   simulation length                  (default 30)\n"
        "  -T thr       congested when free ports < thr    (default 2)\n"
        "  -W seconds   base station window W              (default 5)\n"
        "  -N 1|2|3     node version: 1 blocking p2p, 2 non-blocking on demand,\n"
        "               3 neighbourhood collective         (default 2)\n"
        "  -B 1|2|3     base version: 1 probe+blocking, 2 Irecv/Test,\n"
        "               3 persistent batch                 (default 2)\n"
        "  -p us        polling interval in microseconds   (default 200)\n"
        "  -R n         immediate re-polls before sleeping (default 1, 0 = off)\n"
        "  -s seed      random seed                        (default 42)\n"
        "  -a p         arrival probability per cycle      (default 0.5)\n"
        "  -d p         departure probability per cycle    (default 0.3)\n"
        "  -H list      hot nodes, e.g. 0,1,3              (default none)\n"
        "  -A p         arrival probability of hot nodes   (default 0.95)\n"
        "  -D p         departure probability of hot nodes (default 0.05)\n"
        "  -C cycles    cool-down after an episode         (default 2)\n"
        "  -o dir       output directory                   (default out)\n"
        "  -l 0|1       per-node logs                      (default 1)\n"
        "  -q           quiet console output\n");
}

static int to_int(const char *s, int *out)
{
    char *end;
    errno = 0;
    long v = strtol(s, &end, 10);
    if (errno || *end || end == s) return -1;
    *out = (int)v;
    return 0;
}

static int to_double(const char *s, double *out)
{
    char *end;
    errno = 0;
    double v = strtod(s, &end);
    if (errno || *end || end == s) return -1;
    *out = v;
    return 0;
}

int parse_args(int argc, char **argv, config_t *cfg, char *err, size_t errlen)
{
    memset(cfg, 0, sizeof(*cfg));
    cfg->k = 6; cfg->cycle_ms = 1000; cfg->sim_s = 30; cfg->threshold = 2;
    cfg->window_s = 5; cfg->node_ver = 2; cfg->base_ver = 2; cfg->poll_us = 200; cfg->recheck = 1;
    cfg->seed = 42; cfg->p_arrive = 0.5; cfg->p_depart = 0.3; cfg->p_arrive_hot = 0.95; cfg->p_depart_hot = 0.05;
    cfg->cooldown = 2; cfg->node_logs = 1;
    strcpy(cfg->outdir, "out");

    int npos = 0, pos[2] = {0, 0};
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (a[0] != '-' || a[1] == '\0' || a[2] != '\0') {
            if (npos >= 2 || to_int(a, &pos[npos]) != 0) {
                snprintf(err, errlen, "unexpected argument '%s'", a);
                return -1;
            }
            npos++;
            continue;
        }
        char opt = a[1];
        if (opt == 'q') { cfg->quiet = 1; continue; }
        if (opt == 'h') { snprintf(err, errlen, "help"); return -1; }
        if (i + 1 >= argc) { snprintf(err, errlen, "option -%c needs a value", opt); return -1; }
        const char *v = argv[++i];
        int bad = 0;
        switch (opt) {
        case 'k': bad = to_int(v, &cfg->k); break;
        case 'c': bad = to_int(v, &cfg->cycle_ms); break;
        case 't': bad = to_double(v, &cfg->sim_s); break;
        case 'T': bad = to_int(v, &cfg->threshold); break;
        case 'W': bad = to_double(v, &cfg->window_s); break;
        case 'N': bad = to_int(v, &cfg->node_ver); break;
        case 'B': bad = to_int(v, &cfg->base_ver); break;
        case 'p': bad = to_int(v, &cfg->poll_us); break;
        case 'R': bad = to_int(v, &cfg->recheck); break;
        case 's': { int s = 0; bad = to_int(v, &s); cfg->seed = (unsigned)s; break; }
        case 'a': bad = to_double(v, &cfg->p_arrive); break;
        case 'd': bad = to_double(v, &cfg->p_depart); break;
        case 'A': bad = to_double(v, &cfg->p_arrive_hot); break;
        case 'D': bad = to_double(v, &cfg->p_depart_hot); break;
        case 'H': snprintf(cfg->hot, sizeof cfg->hot, "%s", v); break;
        case 'C': bad = to_int(v, &cfg->cooldown); break;
        case 'o': snprintf(cfg->outdir, sizeof cfg->outdir, "%s", v); break;
        case 'l': bad = to_int(v, &cfg->node_logs); break;
        default: snprintf(err, errlen, "unknown option -%c", opt); return -1;
        }
        if (bad) { snprintf(err, errlen, "bad value '%s' for -%c", v, opt); return -1; }
    }
    if (npos == 1) { snprintf(err, errlen, "give both m and n"); return -1; }
    if (npos == 2) { cfg->m = pos[0]; cfg->n = pos[1]; }

    if (npos == 2 && (cfg->m < 1 || cfg->n < 1)) { snprintf(err, errlen, "m and n must be >= 1"); return -1; }
    if (cfg->k < 1 || cfg->k > MAX_PORTS) { snprintf(err, errlen, "k must be 1..%d", MAX_PORTS); return -1; }
    if (cfg->cycle_ms < 10) { snprintf(err, errlen, "cycle must be >= 10 ms"); return -1; }
    if (cfg->sim_s <= 0 || cfg->window_s < 0) { snprintf(err, errlen, "bad -t or -W"); return -1; }
    if (cfg->threshold < 1 || cfg->threshold > cfg->k) { snprintf(err, errlen, "threshold must be 1..k"); return -1; }
    if (cfg->node_ver < 1 || cfg->node_ver > 3 || cfg->base_ver < 1 || cfg->base_ver > 3) {
        snprintf(err, errlen, "versions must be 1..3"); return -1;
    }
    if (cfg->poll_us < 0 || cfg->cooldown < 0 || cfg->recheck < 0) { snprintf(err, errlen, "bad -p or -C"); return -1; }
    return 0;
}

void stamp_now(int ts[7], double *epoch)
{
    struct timespec t;
    struct tm tmv;
    clock_gettime(CLOCK_REALTIME, &t);
    localtime_r(&t.tv_sec, &tmv);
    ts[0] = tmv.tm_year + 1900; ts[1] = tmv.tm_mon + 1; ts[2] = tmv.tm_mday;
    ts[3] = tmv.tm_hour; ts[4] = tmv.tm_min; ts[5] = tmv.tm_sec;
    ts[6] = (int)(t.tv_nsec / 1000000);
    if (epoch) *epoch = (double)t.tv_sec + t.tv_nsec * 1e-9;
}

void stamp_fmt(const int ts[7], char *buf, size_t len)
{
    snprintf(buf, len, "%04d-%02d-%02d %02d:%02d:%02d.%03d",
             ts[0], ts[1], ts[2], ts[3], ts[4], ts[5], ts[6]);
}

int mkdir_p(const char *path)
{
    char tmp[512];
    snprintf(tmp, sizeof tmp, "%s", path);
    for (char *p = tmp + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            if (mkdir(tmp, 0755) != 0 && errno != EEXIST) return -1;
            *p = '/';
        }
    }
    if (mkdir(tmp, 0755) != 0 && errno != EEXIST) return -1;
    return 0;
}

void sleep_us(long us)
{
    if (us <= 0) return;
    struct timespec t = { us / 1000000, (us % 1000000) * 1000 };
    while (nanosleep(&t, &t) != 0 && errno == EINTR) { }
}

double mono_now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + t.tv_nsec * 1e-9;
}

/* Sleep only after (recheck + 1) consecutive idle polls. Open MPI's MPI_Iprobe and
 * MPI_Testsome drive the progress engine *after* checking, so a message that lands
 * during the call is only reported by the next call; polling once more before
 * sleeping removes that extra sleep period from the detection latency. */
void idle_wait(const config_t *cfg, int progress, int *idle)
{
    if (progress) { *idle = 0; return; }
    if (++*idle > cfg->recheck) {
        sleep_us(cfg->poll_us);
        *idle = 0;
    }
}

int is_hot(const config_t *cfg, int cart_rank)
{
    const char *p = cfg->hot;
    while (*p) {
        char *end;
        long v = strtol(p, &end, 10);
        if (end == p) break;
        if (v == cart_rank) return 1;
        p = (*end == ',') ? end + 1 : end;
    }
    return 0;
}

/* One derived datatype for the whole report: a block of doubles followed by a
 * block of ints, so a congestion report is a single MPI message. */
MPI_Datatype make_report_type(void)
{
    _Static_assert(offsetof(report_t, alert_seq) + sizeof(int)
                   == offsetof(report_t, ts) + REPORT_INTS * sizeof(int),
                   "REPORT_INTS does not match report_t");
    int blocklen[2] = { REPORT_DOUBLES, REPORT_INTS };
    MPI_Aint disp[2] = { offsetof(report_t, t_wall), offsetof(report_t, ts) };
    MPI_Datatype types[2] = { MPI_DOUBLE, MPI_INT };
    MPI_Datatype tmp, rt;
    MPI_Type_create_struct(2, blocklen, disp, types, &tmp);
    MPI_Type_create_resized(tmp, 0, sizeof(report_t), &rt);
    MPI_Type_free(&tmp);
    MPI_Type_commit(&rt);
    return rt;
}
