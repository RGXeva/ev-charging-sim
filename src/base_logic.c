/* SPDX-License-Identifier: MIT */
/* Base station decision logic, kept free of MPI so it can be unit tested. */
#include "common.h"

#include <stdlib.h>

/*
 * Search outwards from (row, col) in Manhattan rings d = 2, 3, ... and return
 * ALL nodes on the first ring that contains a node not known to be busy
 * (busy_until[rank] <= now). Distance 1 is skipped: the reporting node's
 * neighbours have just reported themselves busy in the same report.
 * Candidates come out in ascending cart-rank order (row-major traversal).
 * Returns the number of candidates (0 = no free node anywhere).
 */
int find_candidates(int m, int n, int row, int col, const double *busy_until,
                    double now, int *cand, int max_cand, int *distance)
{
    int maxd = (m - 1) + (n - 1);
    for (int d = 2; d <= maxd; d++) {
        int nc = 0;
        for (int r = 0; r < m; r++) {
            int dr = abs(r - row);
            if (dr > d) continue;
            int dc = d - dr;
            int cols[2] = { col - dc, col + dc };
            int ncols = (dc == 0) ? 1 : 2;
            for (int i = 0; i < ncols; i++) {
                int c = cols[i];
                if (c < 0 || c >= n) continue;
                int rank = r * n + c;
                if (busy_until[rank] <= now && nc < max_cand) cand[nc++] = rank;
            }
        }
        if (nc > 0) { *distance = d; return nc; }
    }
    *distance = 0;
    return 0;
}
