/* SPDX-License-Identifier: MIT */
/* Unit tests for the base-station candidate search (no MPI calls). */
#include "common.h"

#include <stdio.h>
#include <string.h>

static int failures = 0;

static void expect(const char *name, int m, int n, int row, int col, const double *busy,
                   const int *want, int nwant, int want_dist)
{
    int cand[MAX_CAND], dist;
    int nc = find_candidates(m, n, row, col, busy, 10.0, cand, MAX_CAND, &dist);
    int ok = nc == nwant && (nwant == 0 || dist == want_dist);
    for (int i = 0; ok && i < nc; i++) ok = cand[i] == want[i];
    printf("%-44s %s (got %d candidate(s) at distance %d:", name, ok ? "PASS" : "FAIL", nc, dist);
    for (int i = 0; i < nc; i++) printf(" %d", cand[i]);
    printf(")\n");
    if (!ok) failures++;
}

int main(void)
{
    double busy[64];
    for (int i = 0; i < 64; i++) busy[i] = -1e30;   /* nothing known busy */

    /* Reference example, 3x3: node 0 reports, neighbours 1 and 3 busy -> {2, 4, 6} */
    busy[0] = busy[1] = busy[3] = 15.0;
    { int w[] = {2, 4, 6}; expect("3x3 corner node 0, reference example", 3, 3, 0, 0, busy, w, 3, 2); }

    /* node 1 already reported earlier (its report showed 2 and 4 busy) -> only 6 left */
    busy[2] = busy[4] = 15.0;
    { int w[] = {6}; expect("3x3 node 0, 2 and 4 known busy", 3, 3, 0, 0, busy, w, 1, 2); }

    /* whole distance-2 ring busy -> next ring (distance 3: 5, 7) */
    busy[6] = 15.0;
    { int w[] = {5, 7}; expect("3x3 node 0, ring 2 busy -> ring 3", 3, 3, 0, 0, busy, w, 2, 3); }

    /* expired knowledge (busy_until <= now) counts as available again */
    for (int i = 0; i < 64; i++) busy[i] = 5.0;
    busy[4] = busy[1] = busy[3] = busy[5] = busy[7] = 15.0;
    { int w[] = {0, 2, 6, 8}; expect("3x3 centre node 4, expired entries", 3, 3, 1, 1, busy, w, 4, 2); }

    /* everything busy -> no candidate */
    for (int i = 0; i < 64; i++) busy[i] = 15.0;
    expect("3x3 all busy", 3, 3, 1, 1, busy, NULL, 0, 0);

    /* line topology 1x6, node 2 with neighbours 1 and 3 busy -> {0, 4} */
    for (int i = 0; i < 64; i++) busy[i] = -1e30;
    busy[1] = busy[2] = busy[3] = 15.0;
    { int w[] = {0, 4}; expect("1x6 line node 2", 1, 6, 0, 2, busy, w, 2, 2); }

    /* rectangular 4x5, centre (1,2) = 7 -> ring 2 in row-major order */
    for (int i = 0; i < 64; i++) busy[i] = -1e30;
    busy[7] = busy[2] = busy[12] = busy[6] = busy[8] = 15.0;
    { int w[] = {1, 3, 5, 9, 11, 13, 17}; expect("4x5 node 7 ring 2", 4, 5, 1, 2, busy, w, 7, 2); }

    printf("%s: %d failure(s)\n", failures ? "FAILED" : "ALL PASSED", failures);
    return failures != 0;
}
