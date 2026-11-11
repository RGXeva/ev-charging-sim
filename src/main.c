/* SPDX-License-Identifier: MIT */
/*
 * Distributed wireless sensor network simulation of an EV charging grid
 * Open MPI (processes, virtual topology) + OpenMP (charging ports).
 *
 *   mpirun -np <m*n+1> ./wsn m n [options]      (see -h)
 */
#include "common.h"

#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv)
{
    int provided, size, wrank;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &wrank);

    config_t cfg;
    char err[256];
    int bad = 0;
    if (provided < MPI_THREAD_FUNNELED) {
        if (wrank == 0) fprintf(stderr, "error: MPI library does not provide MPI_THREAD_FUNNELED\n");
        bad = 1;
    } else if (parse_args(argc, argv, &cfg, err, sizeof err) != 0) {
        if (wrank == 0) {
            if (strcmp(err, "help") != 0) fprintf(stderr, "error: %s\n", err);
            print_usage(stderr);
        }
        bad = 1;
    } else if (size < 2) {
        if (wrank == 0) fprintf(stderr, "error: need at least 2 processes (1 node + 1 base station)\n");
        bad = 1;
    } else if (cfg.m == 0) {
        int dims[2] = {0, 0};
        MPI_Dims_create(size - 1, 2, dims);
        cfg.m = dims[0];
        cfg.n = dims[1];
    } else if (cfg.m * cfg.n != size - 1) {
        if (wrank == 0)
            fprintf(stderr, "error: m*n = %d x %d = %d nodes needs %d processes, got %d\n",
                    cfg.m, cfg.n, cfg.m * cfg.n, cfg.m * cfg.n + 1, size);
        bad = 1;
    }
    if (bad) {                        /* every rank reached the same decision */
        MPI_Finalize();
        return 1;
    }

    int n_nodes = size - 1, base_world = size - 1;
    int is_base = (wrank == base_world);

    /* split WORLD into {charging nodes} and {base station} */
    MPI_Comm node_comm;
    MPI_Comm_split(MPI_COMM_WORLD, is_base, wrank, &node_comm);
    MPI_Datatype report_type = make_report_type();

    int rc;
    if (is_base) {
        rc = run_base(&cfg, n_nodes, report_type);
    } else {
        int dims[2] = {cfg.m, cfg.n}, periods[2] = {0, 0};
        MPI_Comm comm2D;
        MPI_Cart_create(node_comm, 2, dims, periods, 1, &comm2D);
        rc = run_node(&cfg, comm2D, base_world, report_type);
        MPI_Comm_free(&comm2D);
    }

    MPI_Type_free(&report_type);
    MPI_Comm_free(&node_comm);
    MPI_Finalize();
    return rc;
}
