/* SPDX-License-Identifier: MIT */
/*
 * Distributed wireless sensor network simulation of an EV charging grid
 * Shared definitions: configuration, message tags, message layouts, statistics.
 *
 * Process layout (MPI_COMM_WORLD, size = m*n + 1):
 *   world ranks 0 .. m*n-1  -> EV charging nodes (2-D Cartesian communicator)
 *   world rank  m*n         -> base station (single process, single MPI thread)
 */
#ifndef WSN_COMMON_H
#define WSN_COMMON_H

#include <mpi.h>
#include <stdio.h>

#define MAX_NBR     4      /* up, down, left, right */
#define MAX_PORTS   64     /* fixed upper bound for the shared port array */
#define MAX_CAND    128    /* upper bound on candidates in one redirect */

/* ---- message tags (numbering kept from the Applied 1 design) ---- */
#define TAG_NBR_REQ   101  /* node -> neighbour : request availability      */
#define TAG_NBR_RESP  102  /* neighbour -> node : availability record       */
#define TAG_ALERT     202  /* node -> base      : congestion report         */
#define TAG_REDIRECT  203  /* base -> node      : nearest available nodes   */
#define TAG_TERM      204  /* base -> node      : terminate simulation      */
#define TAG_FINAL     205  /* node -> base      : final per-node statistics */
#define TAG_SAMPLES   206  /* node -> base      : timing samples (shutdown) */

/* redirect status codes */
#define RD_CANDIDATES    1   /* candidate list attached                       */
#define RD_NONE_NEARBY   0   /* all neighbours reported within window W       */
#define RD_NONE_IN_GRID -1   /* no node outside the known-busy set exists     */

typedef struct {
    int m, n;              /* grid rows, columns                              */
    int k;                 /* charging ports per node (OpenMP threads)        */
    int cycle_ms;          /* port update / availability sampling cycle       */
    double sim_s;          /* simulation length, enforced by the base station */
    int threshold;         /* congested when free ports < threshold           */
    double window_s;       /* base station "predefined period" W              */
    int node_ver;          /* 1 blocking p2p, 2 non-blocking on demand, 3 neighbourhood collective */
    int base_ver;          /* 1 probe+blocking, 2 Irecv/Test, 3 persistent batch */
    int poll_us;           /* polling interval of event loops                 */
    int recheck;           /* extra immediate re-polls before sleeping (0 = off) */
    unsigned seed;
    double p_arrive;       /* per cycle probability a free port becomes used  */
    double p_depart;       /* per cycle probability a used port becomes free  */
    double p_arrive_hot;   /* p_arrive used by "hot" nodes                    */
    double p_depart_hot;   /* p_depart used by "hot" nodes                    */
    char hot[256];         /* comma separated cart ranks of hot nodes         */
    int cooldown;          /* cycles to wait after an episode before the next */
    char outdir[256];
    int node_logs;         /* 1: each node writes node_<rank>.log at shutdown */
    int quiet;
} config_t;

/* availability record exchanged between neighbours (6 ints) */
typedef struct {
    int rank, row, col, free_ports, total_ports, cycle;
} nbr_info_t;
#define NBR_INFO_INTS 6

/* congestion report node -> base, sent as one derived-datatype message */
typedef struct {
    double t_wall;         /* epoch seconds when congestion was confirmed     */
    double nbr_rtt;        /* time of the neighbour exchange (s)              */
    double t_reserved;
    int ts[7];             /* Y M D h m s ms (local time)                     */
    int cycle;
    int cart_rank, world_rank, row, col;
    int total_ports, free_ports, threshold;
    int n_nbrs;
    int nbr_rank[MAX_NBR], nbr_row[MAX_NBR], nbr_col[MAX_NBR];
    int nbr_free[MAX_NBR], nbr_total[MAX_NBR];
    int msgs_exchanged;    /* neighbour messages for this episode (sent+recv) */
    int bytes_exchanged;   /* payload bytes for this episode                  */
    int alert_seq;         /* sequence number of this node's alert            */
} report_t;
#define REPORT_DOUBLES 3
#define REPORT_INTS    (7 + 1 + 4 + 3 + 1 + 5 * MAX_NBR + 3)

/* redirect base -> node: int array, variable length 5 + 3*n_cand */
#define RD_HDR 5           /* status, distance, n_cand, alert_seq, service_us */
#define RD_MAX_INTS (RD_HDR + 3 * MAX_CAND)

/* per-node statistics sent to the base in TAG_FINAL (array of doubles) */
enum {
    ST_CYCLES, ST_ALERTS, ST_LOCAL_REDIRECT, ST_BASE_CAND, ST_BASE_NONE,
    ST_EXCHANGES, ST_REQ_SENT, ST_REQ_RECV, ST_RESP_SENT, ST_RESP_RECV,
    ST_MSG_SENT, ST_MSG_RECV, ST_BYTES_SENT, ST_BYTES_RECV,
    ST_NBR_RTT_SUM, ST_NBR_RTT_MAX, ST_ALERT_RTT_SUM, ST_ALERT_NET_SUM,
    ST_ALERT_RTT_MAX, ST_SYNC_SUM, ST_MIN_FREE, ST_SUM_FREE, ST_PORT_CHANGES,
    ST_COUNT
};
extern const char *STAT_NAMES[ST_COUNT];

/* timing sample types (TAG_SAMPLES carries triplets: type, v1, v2) */
#define SMP_NBR    1       /* v1 = neighbour exchange time, v2 = #neighbours  */
#define SMP_ALERT  2       /* v1 = alert->redirect RTT, v2 = RTT minus base service */
#define SMP_SYNC   3       /* v1 = termination-agreement Allreduce time (N1/N3) */

/* ---- util.c ---- */
int  parse_args(int argc, char **argv, config_t *cfg, char *err, size_t errlen);
void print_usage(FILE *f);
void stamp_now(int ts[7], double *epoch);
void stamp_fmt(const int ts[7], char *buf, size_t len);
int  mkdir_p(const char *path);
void sleep_us(long us);
double mono_now(void);
void idle_wait(const config_t *cfg, int progress, int *idle);
int  is_hot(const config_t *cfg, int cart_rank);
MPI_Datatype make_report_type(void);

/* ---- node.c / base.c ---- */
int run_node(const config_t *cfg, MPI_Comm comm2D, int base_world, MPI_Datatype report_type);
int run_base(const config_t *cfg, int n_nodes, MPI_Datatype report_type);

/* ---- base_logic.c (pure, unit-testable) ---- */
int find_candidates(int m, int n, int row, int col, const double *busy_until,
                    double now, int *cand, int max_cand, int *distance);

#endif
