#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Generate the experiment case files.

usage: python3 -I scripts/make_cases.py [--best "-N 2 -B 3"]

Each line: <experiment> <label> <rep> <m> <n> <program options...>
Repetitions are interleaved (rep 1 of every case, then rep 2, ...) so slow drift on
a shared machine does not bias one configuration.

Two targets:
  jobs/cases/local_*.txt  the 30+ grid sizes, 30+ line sizes and parameter sweeps,
                          run on the local machine (scripts/run_cases.sh, LAUNCHER=mpirun)
  jobs/cases/caas_*.txt   a compact set (5 grid sizes, 5 line sizes, version check)
                          that fits in ONE 20-minute SLURM job on 5 x 16 cores
"""
import argparse, os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "jobs", "cases")
REPS = 3
RUN_S = 10            # wall time per run (-t)

# 200 ms cycle (time-compressed: one cycle stands for a 10 s update period),
# moderately congested load so every message type occurs in every run.
COMMON = f"-c 200 -t {RUN_S} -a 0.6 -d 0.25 -l 0 -q"

# 33 grid shapes, 4 .. 64 charging nodes, square and elongated
GRID = [(2, 2), (2, 3), (3, 3), (2, 5), (3, 4), (2, 7), (4, 4), (3, 6), (4, 5), (3, 7), (2, 11),
        (4, 6), (5, 5), (3, 9), (4, 7), (5, 6), (2, 15), (4, 8), (3, 11), (5, 7), (6, 6), (4, 9),
        (5, 8), (6, 7), (3, 14), (4, 11), (6, 8), (7, 7), (5, 10), (6, 9), (7, 8), (4, 14), (8, 8)]
# 35 line sizes: every N up to 31, then the node counts of the CAAS square grids
LINE = list(range(2, 32)) + [36, 42, 49, 56, 64]
CAAS_GRID = [(3, 3), (4, 4), (5, 5), (6, 6), (8, 8)]


def local_experiments(best):
    e = {}
    e["E1_grid"] = [(f"{m}x{n}", m, n, f"{COMMON} {best}") for (m, n) in GRID]
    e["E2_line"] = [(f"1x{n}", 1, n, f"{COMMON} {best}") for n in LINE]
    e["E3_versions"] = [(f"{m}x{n}_N{nv}B{bv}", m, n, f"{COMMON} -N {nv} -B {bv}")
                        for (m, n) in [(3, 3), (4, 5), (5, 6), (1, 16)]
                        for nv in (1, 2, 3) for bv in (1, 2, 3)]
    e["E4_window"] = [(f"5x6_W{w}", 5, 6, f"{COMMON} {best} -W {w}") for w in (0.2, 0.5, 1, 2, 5)]
    e["E5_poll"] = [(f"5x6_p{p}", 5, 6, f"{COMMON} {best} -p {p}") for p in (20, 100, 200, 1000, 5000)]
    e["E6_load"] = [(f"5x6_a{a}", 5, 6, f"-c 200 -t {RUN_S} -a {a} -d 0.25 -l 0 -q {best}")
                    for a in (0.2, 0.4, 0.6, 0.8, 0.95)]
    e["E7_recheck"] = [(f"5x6_B{bv}_R{r}", 5, 6, f"{COMMON} -N 2 -B {bv} -R {r}")
                       for bv in (1, 2, 3) for r in (0, 1)]
    return e


def caas_experiments(best):
    e = {}
    e["C1_grid"] = [(f"{m}x{n}", m, n, f"{COMMON} {best}") for (m, n) in CAAS_GRID]
    e["C2_line"] = [(f"1x{m * n}", 1, m * n, f"{COMMON} {best}") for (m, n) in CAAS_GRID]
    e["C3_versions"] = [(f"6x6_N{nv}B{bv}", 6, 6, f"{COMMON} -N {nv} -B {bv}")
                        for nv in (1, 2, 3) for bv in (1, 2, 3)]
    e["C4_recheck"] = [(f"6x6_B{bv}_R0", 6, 6, f"{COMMON} -N 2 -B {bv} -R 0") for bv in (1, 2, 3)]
    return e


def write(name, exps, order, reps):
    runs = [f"{x} {label} {rep} {m} {n} {args}"
            for rep in range(1, reps + 1) for x in order for (label, m, n, args) in exps[x]]
    path = os.path.join(OUT, name + ".txt")
    with open(path, "w") as f:
        f.write("# experiment label rep m n options\n" + "\n".join(runs) + "\n")
    print(f"{os.path.relpath(path)}: {len(runs)} runs, ~{len(runs) * (RUN_S + 3) / 60:.0f} min")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--best", default="-N 2 -B 3", help="versions used by the scaling experiments")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    for f in os.listdir(OUT):
        if f.endswith(".txt"):
            os.remove(os.path.join(OUT, f))
    loc, caas = local_experiments(a.best), caas_experiments(a.best)
    write("local_1_versions", loc, ["E3_versions", "E7_recheck"], REPS)
    write("local_2_scaling", loc, ["E1_grid", "E2_line"], REPS)
    write("local_3_params", loc, ["E4_window", "E5_poll", "E6_load"], REPS)
    # CAAS: scaling with 3 repetitions, version checks with 2 -> one 20-minute job
    write("caas_1_scaling", caas, ["C1_grid", "C2_line"], REPS)
    write("caas_2_versions", caas, ["C3_versions", "C4_recheck"], 2)
    for name, ex in (("local", loc), ("caas", caas)):
        print(name + ": " + ", ".join(f"{k} {len(v)} configs" for k, v in ex.items()))


if __name__ == "__main__":
    main()
