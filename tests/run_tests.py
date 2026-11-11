#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Correctness tests for the WSN simulator.

Each case launches the real MPI program and checks: exit code, timeout (deadlock),
graceful shutdown of every node, no unexpected messages or shutdown warnings,
and scenario-specific expectations read from the base station summary.

usage: python3 -I tests/run_tests.py [--launcher "mpirun --oversubscribe"] [--quick]
"""
import argparse, csv, os, shlex, shutil, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(ROOT, "bin", "wsn")
OUT = os.path.join(ROOT, "results", "tests")


def run(launcher, m, n, args, name, timeout=90, nprocs=None):
    out = os.path.join(OUT, name)
    shutil.rmtree(out, ignore_errors=True)
    np_ = nprocs if nprocs is not None else m * n + 1
    cmd = shlex.split(launcher) + ["-n", str(np_), BIN, str(m), str(n), "-o", out, "-q"] + args
    t = time.time()
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           universal_newlines=True, timeout=timeout)   # Python 3.6 compatible
        rc, err = p.returncode, p.stderr
    except subprocess.TimeoutExpired:
        return {"name": name, "rc": None, "err": "TIMEOUT (deadlock?)", "out": out, "secs": timeout}
    return {"name": name, "rc": rc, "err": err, "out": out, "secs": time.time() - t}


def summary(res):
    path = os.path.join(res["out"], "summary.csv")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return next(csv.DictReader(f))


def common_checks(res, nodes):
    problems = []
    if res["rc"] is None:
        return [res["err"]]
    if res["rc"] != 0:
        problems.append(f"exit code {res['rc']}: {res['err'].strip()[-300:]}")
    s = summary(res)
    if s is None:
        return problems + ["summary.csv missing"]
    log = open(os.path.join(res["out"], "base_log.txt")).read()
    if f"all {nodes} nodes confirmed shutdown" not in log:
        problems.append("not every node confirmed shutdown")
    if int(s["unexpected"]) != 0:
        problems.append(f"{s['unexpected']} unexpected messages")
    for fn in os.listdir(res["out"]):
        if fn.startswith("node_") and "warning" in open(os.path.join(res["out"], fn)).read():
            problems.append(f"warning in {fn}")
    if "rror" in res["err"]:
        problems.append("stderr: " + res["err"].strip()[-300:])
    return problems


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--launcher", default="mpirun --oversubscribe")
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    fast = ["-c", "100", "-t", "3"]
    cases = []

    # 1. every node version x base version, congested grid so all paths run
    for nv in (1, 2, 3):
        for bv in (1, 2, 3):
            cases.append((f"v_N{nv}_B{bv}", 3, 3, fast + ["-N", str(nv), "-B", str(bv), "-a", "0.8", "-d", "0.2"],
                          lambda s: int(s["alerts"]) > 0 and int(s["local_redirects"]) >= 0, "alerts > 0"))

    # 2. topologies: single node, lines, rectangles
    shapes = [(1, 1), (1, 2), (1, 8), (2, 2), (2, 5), (4, 5)] if not a.quick else [(1, 1), (1, 4), (2, 3)]
    for (m, n) in shapes:
        for nv in (1, 2, 3):
            cases.append((f"shape_{m}x{n}_N{nv}", m, n, fast + ["-N", str(nv), "-a", "0.7"],
                          lambda s: True, "runs and shuts down"))

    # 3. low availability, neighbours free: only the centre is hot -> local redirects, no base alert
    cases.append(("scenario_neighbour_free", 3, 3,
                  ["-c", "100", "-t", "3", "-H", "4", "-A", "1.0", "-D", "0.0", "-a", "0.0", "-d", "1.0", "-C", "0"],
                  lambda s: int(s["local_redirects"]) > 0 and int(s["alerts"]) == 0,
                  "local redirects > 0 and no alerts"))

    # 4. hot corner block, rest free -> base must suggest candidates (reference example: node 0, 1, 3)
    cases.append(("scenario_base_candidates", 3, 3,
                  ["-c", "100", "-t", "3", "-H", "0,1,3", "-A", "1.0", "-D", "0.0", "-a", "0.0", "-d", "1.0", "-W", "0.15"],
                  lambda s: int(s["cand"]) > 0, "base replies with candidates"))

    # 5. everything saturated, long window -> neighbours keep reporting -> 'none nearby'
    cases.append(("scenario_none_nearby", 3, 3,
                  ["-c", "100", "-t", "3", "-H", "0,1,2,3,4,5,6,7,8", "-A", "1.0", "-D", "0.0", "-W", "10", "-C", "1"],
                  lambda s: int(s["none_nearby"]) > 0, "base replies 'none nearby'"))

    results, failed = [], 0
    for (name, m, n, args, check, desc) in cases:
        res = run(a.launcher, m, n, args, name)
        probs = common_checks(res, m * n)
        s = summary(res)
        if s is not None and not probs and not check(s):
            probs.append(f"expectation failed: {desc}")
        ok = not probs
        failed += not ok
        extra = f"alerts={s['alerts']} cand={s['cand']} none={s['none_nearby']} local={s['local_redirects']}" if s else ""
        print(f"{'PASS' if ok else 'FAIL'}  {name:28s} {res['secs']:5.1f}s  {extra}  {'; '.join(probs)}", flush=True)

    # 6. wrong process count must fail cleanly (own message, exit != 0, no output)
    res = run(a.launcher, 3, 3, [], "neg_wrong_np", nprocs=5, timeout=30)
    ok = res["rc"] not in (None, 0) and "needs 10 processes" in res["err"] and summary(res) is None
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {'neg_wrong_np':28s} rc={res['rc']}")

    res = run(a.launcher, 2, 2, ["-k", "0"], "neg_bad_k", timeout=30)
    ok = res["rc"] not in (None, 0) and "k must be" in res["err"] and summary(res) is None
    failed += not ok
    print(f"{'PASS' if ok else 'FAIL'}  {'neg_bad_k':28s} rc={res['rc']}")

    print(f"\n{'ALL PASSED' if failed == 0 else f'{failed} FAILED'}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
