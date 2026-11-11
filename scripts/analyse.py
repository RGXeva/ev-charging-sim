#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Turn raw run directories into the tables and figures under results/.

usage: python3 -I scripts/analyse.py [raw_dir] [out_dir]
       (defaults: results/raw -> results/analysis)

Every number comes from the raw files written by the base station:
  summary.csv  one row per run (counts, per-run medians, shutdown time, base CPU)
  samples.csv  every timing sample reported by the nodes at shutdown
Delays are pooled over the repetitions of a configuration (median + IQR).
"""
import csv, glob, math, os, sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results", "raw")
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "results", "analysis")
FIG, TAB = os.path.join(OUT, "figures"), os.path.join(OUT, "tables")

# reference categorical palette (fixed order) and text/grid tokens, light mode
S1, S2, S3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 10,
    "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlecolor": INK,
    "legend.frameon": False, "lines.linewidth": 2, "lines.markersize": 6,
})

NBR, ALERT, SYNC = 1, 2, 3


# ------------------------------------------------------------------ loading
def load_runs():
    runs = []
    for sp in sorted(glob.glob(os.path.join(RAW, "*", "*", "summary.csv"))):
        d = os.path.dirname(sp)
        exp, run = os.path.basename(os.path.dirname(d)), os.path.basename(d)
        label, rep = run.rsplit("_r", 1)
        with open(sp) as f:
            s = next(csv.DictReader(f))
        smp = defaultdict(list)
        with open(os.path.join(d, "samples.csv")) as f:
            for row in csv.DictReader(f):
                t = int(row["type"])
                smp[t].append(float(row["v1"]) * 1e3)
                if t == ALERT:
                    smp["net"].append(float(row["v2"]) * 1e3)
        runs.append({"exp": exp, "label": label, "rep": int(rep), "s": s, "smp": smp})
    return runs


def pct(v, p):
    if not v:
        return float("nan")
    v = sorted(v)
    return v[max(0, min(len(v) - 1, math.ceil(p * len(v)) - 1))]


def group(runs, exp):
    g = defaultdict(list)
    for r in runs:
        if r["exp"] == exp:
            g[r["label"]].append(r)
    return g


def stats(rs):
    """pooled delay statistics + mean counts for one configuration"""
    s0 = rs[0]["s"]
    pool = defaultdict(list)
    for r in rs:
        for k, v in r["smp"].items():
            pool[k].extend(v)
    num = lambda key: [float(r["s"][key]) for r in rs if r["s"][key] not in ("nan", "")]
    mean = lambda xs: sum(xs) / len(xs) if xs else float("nan")
    nodes = int(s0["nodes"])
    cycles = float(s0["sim_s"]) * 1000 / int(s0["cycle_ms"])
    out = {"m": int(s0["m"]), "n": int(s0["n"]), "nodes": nodes, "reps": len(rs),
           "node_ver": int(s0["node_ver"]), "base_ver": int(s0["base_ver"]),
           "window_s": float(s0["window_s"]), "poll_us": int(s0["poll_us"]), "p_arrive": float(s0["p_arrive"])}
    for key, name in ((NBR, "nbr"), (ALERT, "alert"), ("net", "net"), (SYNC, "sync")):
        v = pool.get(key, [])
        out[f"{name}_n"] = len(v)
        out[f"{name}_p25"], out[f"{name}_med"] = pct(v, 0.25), pct(v, 0.5)
        out[f"{name}_p75"], out[f"{name}_p95"] = pct(v, 0.75), pct(v, 0.95)
    sd = num("shutdown_ms")
    out["shutdown_med"], out["shutdown_min"], out["shutdown_max"] = pct(sd, 0.5), min(sd), max(sd)
    for key in ("alerts", "cand", "none_nearby", "none_grid", "local_redirects", "exchanges",
                "total_msgs", "total_bytes", "svc_mean_us", "base_cpu_s"):
        out[key] = mean(num(key))
    out["msgs_per_node_cycle"] = out["total_msgs"] / (nodes * cycles)
    return out


def write_table(rows, name, cols):
    os.makedirs(TAB, exist_ok=True)
    with open(os.path.join(TAB, name + ".csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([f"{r[c]:.4g}" if isinstance(r[c], float) else r[c] for c in cols])


def fit(xs, ys):
    """least-squares line y = a + b x and R^2 (ignores NaN)"""
    pts = [(x, y) for x, y in zip(xs, ys) if not math.isnan(y)]
    if len(pts) < 3:
        return float("nan"), float("nan"), float("nan")
    n = len(pts)
    mx, my = sum(p[0] for p in pts) / n, sum(p[1] for p in pts) / n
    sxx = sum((p[0] - mx) ** 2 for p in pts)
    b = sum((p[0] - mx) * (p[1] - my) for p in pts) / sxx if sxx else 0.0
    a = my - b * mx
    ss_tot = sum((p[1] - my) ** 2 for p in pts)
    ss_res = sum((p[1] - a - b * p[0]) ** 2 for p in pts)
    return a, b, (1 - ss_res / ss_tot) if ss_tot else float("nan")


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=180)
    plt.close(fig)


def err_series(ax, rows, x, key, color, label=None, host_line=None):
    xs = [r[x] for r in rows]
    med = [r[f"{key}_med"] for r in rows]
    lo = [max(0, r[f"{key}_med"] - r[f"{key}_p25"]) for r in rows]
    hi = [max(0, r[f"{key}_p75"] - r[f"{key}_med"]) for r in rows]
    ax.errorbar(xs, med, yerr=[lo, hi], fmt="o", color=color, ecolor=color, elinewidth=1,
                capsize=2, markersize=5, markeredgecolor=SURF, markeredgewidth=1, label=label)
    a, b, r2 = fit(xs, med)
    if not math.isnan(b):
        xx = [min(xs), max(xs)]
        ax.plot(xx, [a + b * v for v in xx], color=color, linewidth=1.2, alpha=0.6, linestyle="--")
    for h in host_line or []:
        ax.axvline(h, color=INK2, linewidth=1, linestyle=":")
    return b, r2


# ------------------------------------------------------------- experiments
DELAY_COLS = ["label", "m", "n", "nodes", "reps", "nbr_n", "nbr_med", "nbr_p25", "nbr_p75", "nbr_p95",
              "alert_n", "alert_med", "alert_p95", "net_med", "net_p25", "net_p75", "net_p95",
              "shutdown_med", "shutdown_min", "shutdown_max", "alerts", "local_redirects",
              "total_msgs", "msgs_per_node_cycle", "svc_mean_us"]


def host_lines(rows, hosts):
    """x positions where a run starts to need one more 16-core host (block placement)"""
    if not hosts:
        return []
    top = max(r["nodes"] for r in rows) + 1
    return [16 * k - 0.5 for k in range(1, top // 16 + 1) if 16 * k - 0.5 < top]


def scaling(runs, exp, title, fname, hosts=False):
    g = group(runs, exp)
    if not g:
        return None
    rows = []
    for label, rs in g.items():
        st = stats(rs)
        st["label"] = label
        st["shutdown_p25"], st["shutdown_p75"] = st["shutdown_min"], st["shutdown_max"]
        rows.append(st)
    rows.sort(key=lambda r: (r["nodes"], r["m"]))
    write_table(rows, fname, DELAY_COLS)
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
    panels = [("nbr", "Neighbour request/response\n(round trip, all neighbours)"),
              ("net", "Alert -> redirect\n(round trip minus base processing)"),
              ("shutdown", "Termination\n(TERMINATE sent -> last FINAL received)")]
    fits = {}
    for ax, (key, sub) in zip(axs, panels):
        b, r2 = err_series(ax, rows, "nodes", key, S1, host_line=host_lines(rows, hosts))
        fits[key] = (b, r2)
        ax.set_title(sub, loc="left")
        ax.set_xlabel("charging nodes (MPI processes = nodes + 1)")
        ax.set_ylabel("delay (ms)")
        ax.set_ylim(0, ax.get_ylim()[1] * 1.2)
        if not math.isnan(b):
            ax.text(0.02, 0.96, f"trend {b * 1e3:+.1f} us/node, R²={r2:.2f}", transform=ax.transAxes,
                    va="top", color=INK2, fontsize=9)
    if hosts:
        axs[0].text(0.98, 0.04, "dotted: run needs one more\n16-core CAAS host", transform=axs[0].transAxes,
                    ha="right", color=INK2, fontsize=8)
    fig.suptitle(title, x=0.01, ha="left", fontweight="bold", color=INK)
    save(fig, fname)
    return rows, fits


def compare_grid_line(grid, line, fname, hosts=False):
    if not grid or not line:
        return
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    for ax, key, sub in ((axs[0], "nbr", "Neighbour request/response"), (axs[1], "net", "Alert -> redirect (network part)")):
        for rows, color, name in ((grid, S1, "2-D grid"), (line, S2, "line (1 x N)")):
            xs = [r["nodes"] for r in rows]
            ys = [r[f"{key}_med"] for r in rows]
            ax.plot(xs, ys, "o", color=color, label=name, markeredgecolor=SURF)
            a, b, _ = fit(xs, ys)
            if not math.isnan(b):
                ax.plot([min(xs), max(xs)], [a + b * min(xs), a + b * max(xs)], color=color, linewidth=1.2, alpha=0.6)
                ax.annotate(name, (max(xs), a + b * max(xs)), xytext=(4, 0), textcoords="offset points",
                            color=INK2, fontsize=9, va="center")
        for h in host_lines(grid + line, hosts):
            ax.axvline(h, color=INK2, linewidth=1, linestyle=":")
        ax.set_title(sub + " - median delay", loc="left")
        ax.set_xlabel("charging nodes")
        ax.set_ylabel("delay (ms)")
        ax.set_ylim(0, ax.get_ylim()[1] * 1.15)
        ax.legend(loc="lower right")
    save(fig, fname)


def versions(runs, exp, fname, table):
    g = group(runs, exp)
    if not g:
        return
    rows = []
    for label, rs in g.items():
        st = stats(rs)
        st["label"], st["size"] = label, label.split("_")[0]
        rows.append(st)
    rows.sort(key=lambda r: (r["nodes"], r["node_ver"], r["base_ver"]))
    write_table(rows, table, ["label", "size", "node_ver", "base_ver", "reps", "nbr_med", "nbr_p95",
                                      "alert_med", "alert_p95", "net_med", "sync_med", "shutdown_med",
                                      "alerts", "total_msgs", "msgs_per_node_cycle", "svc_mean_us", "base_cpu_s"])
    sizes = sorted({r["size"] for r in rows}, key=lambda s: next(r["nodes"] for r in rows if r["size"] == s))

    def agg(metric, by, ver):
        out = []
        for s in sizes:
            v = [r[metric] for r in rows if r["size"] == s and r[by] == ver and not math.isnan(r[metric])]
            out.append(sum(v) / len(v) if v else float("nan"))
        return out

    fig, axs = plt.subplots(2, 2, figsize=(12, 7))
    specs = [(axs[0, 0], "nbr_med", "node_ver", "N", "Neighbour exchange median (ms)"),
             (axs[0, 1], "msgs_per_node_cycle", "node_ver", "N", "Messages per node per cycle"),
             (axs[1, 0], "net_med", "base_ver", "B", "Alert -> redirect median, network part (ms)"),
             (axs[1, 1], "base_cpu_s", "base_ver", "B", "Base station CPU time per run (s)")]
    names = {"N": {1: "N1 blocking p2p", 2: "N2 non-blocking on demand", 3: "N3 neighbourhood collective"},
             "B": {1: "B1 probe + blocking", 2: "B2 Irecv/Test", 3: "B3 persistent batch"}}
    w = 0.26
    for ax, metric, by, pre, title in specs:
        for i, (ver, color) in enumerate(((1, S1), (2, S2), (3, S3))):
            ys = agg(metric, by, ver)
            xs = [j + (i - 1) * (w + 0.02) for j in range(len(sizes))]
            ax.bar(xs, ys, width=w, color=color, label=names[pre][ver], edgecolor=SURF, linewidth=1)
        ax.set_xticks(range(len(sizes)))
        ax.set_xticklabels(sizes)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.45)     # headroom for the legend
        ax.set_title(title, loc="left")
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(axis="x", visible=False)
    fig.suptitle("Version comparison (top: node versions averaged over B1-B3; bottom: base versions averaged over N1-N3)",
                 x=0.01, ha="left", fontweight="bold", color=INK, fontsize=11)
    save(fig, fname)


def param_sweep(runs, exp, key, xlabel, fname, logx=False):
    g = group(runs, exp)
    if not g:
        return
    rows = sorted((stats(rs) for rs in g.values()), key=lambda r: r[key])
    for r, label in zip(rows, sorted(g, key=lambda l: stats(g[l])[key])):
        r["label"] = label
    write_table(rows, fname, ["label", key, "reps", "alerts", "cand", "none_nearby", "none_grid", "local_redirects",
                              "nbr_med", "nbr_p95", "alert_med", "net_med", "net_p95", "svc_mean_us", "base_cpu_s",
                              "total_msgs"])
    xs = [r[key] for r in rows]
    if exp == "E4_window":
        fig, ax = plt.subplots(figsize=(6.5, 3.8))
        cats = [("cand", "candidates suggested", S1), ("none_nearby", "'no node nearby'", S2), ("none_grid", "'none in grid'", S3)]
        bottom = [0.0] * len(rows)
        pos = range(len(rows))
        for k, name, color in cats:
            tot = [max(r["alerts"], 1e-9) for r in rows]
            share = [100 * r[k] / t for r, t in zip(rows, tot)]
            ax.bar(pos, share, bottom=bottom, color=color, label=name, edgecolor=SURF, linewidth=2, width=0.6)
            bottom = [b + s for b, s in zip(bottom, share)]
        ax.set_xticks(list(pos))
        ax.set_xticklabels([f"{x:g}" for x in xs])
        ax.set_xlabel(xlabel)
        ax.set_ylabel("share of base replies (%)")
        ax.set_title("Base station decision vs window W (5 x 6 grid)", loc="left")
        ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.0, 1.0))
        ax.grid(axis="x", visible=False)
        save(fig, fname)
        return
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.6))
    for ax, (k, t) in zip(axs, (("nbr_med", "Neighbour exchange median (ms)"), ("net_med", "Alert -> redirect median, network part (ms)"),
                                ("base_cpu_s" if exp == "E5_poll" else "alerts",
                                 "Base station CPU time (s)" if exp == "E5_poll" else "Alerts per run"))):
        ax.plot(xs, [r[k] for r in rows], "-o", color=S1, markeredgecolor=SURF)
        if logx:
            ax.set_xscale("log")
        ax.set_xlabel(xlabel)
        ax.set_title(t, loc="left")
        ax.set_ylim(bottom=0)
    titles = {"E5_poll": "Polling interval trade-off: detection latency vs CPU (5 x 6 grid, N2 / B2)",
              "E6_load": "EV load: alert rate and alert latency (5 x 6 grid, N2 / B2)"}
    fig.suptitle(titles.get(exp, exp), x=0.01, ha="left", fontweight="bold", color=INK)
    save(fig, fname)


def recheck(runs, exp, fname, table, size_label, r1_from=None):
    """R0 vs R1 per base version. r1_from: experiment whose N2 rows supply R1 (CAAS set)"""
    g = group(runs, exp)
    if not g:
        return
    rows = []
    for label, rs in g.items():
        st = stats(rs)
        st["label"], st["recheck"] = label, int(label.rsplit("_R", 1)[1])
        rows.append(st)
    if r1_from:
        for label, rs in group(runs, r1_from).items():
            st = stats(rs)
            if st["node_ver"] == 2:
                st["label"], st["recheck"] = label, 1
                rows.append(st)
    rows.sort(key=lambda r: (r["base_ver"], r["recheck"]))
    write_table(rows, table, ["label", "base_ver", "recheck", "reps", "nbr_med", "nbr_p95", "alert_med",
                                     "alert_p95", "net_med", "svc_mean_us", "base_cpu_s"])
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
    w = 0.38
    for ax, (k, t) in zip(axs, (("alert_med", "Alert -> redirect round trip, median (ms)"),
                                ("nbr_med", "Neighbour exchange, median (ms)"),
                                ("base_cpu_s", "Base station CPU time per run (s)"))):
        for i, (rv, color, name) in enumerate(((0, S1, "sleep after 1 empty poll (R0)"),
                                               (1, S2, "re-poll once, then sleep (R1)"))):
            ys = [next((r[k] for r in rows if r["base_ver"] == bv and r["recheck"] == rv), float("nan")) for bv in (1, 2, 3)]
            ax.bar([j + (i - 0.5) * (w + 0.02) for j in range(3)], ys, width=w, color=color, label=name,
                   edgecolor=SURF, linewidth=1)
        ax.set_xticks(range(3))
        ax.set_xticklabels(["B1\nIprobe", "B2\nIrecv+Test", "B3\npersistent Testsome"])
        ax.set_title(t, loc="left")
        ax.set_ylim(0, ax.get_ylim()[1] * 1.35)
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(axis="x", visible=False)
    fig.suptitle(f"Progress-engine optimisation: one immediate re-poll before sleeping ({size_label} grid, N2)",
                 x=0.01, ha="left", fontweight="bold", color=INK)
    save(fig, fname)


def best_label(runs, exp):
    r = next((r for r in runs if r["exp"] == exp), None)
    return f"N{r['s']['node_ver']} / B{r['s']['base_ver']}" if r else ""


def main():
    runs = load_runs()
    if not runs:
        sys.exit(f"no runs found under {RAW}")
    print(f"{len(runs)} runs loaded from {RAW}")
    fits = {}
    # local set (E*): 30+ grid and line sizes on one machine
    v = best_label(runs, "E1_grid")
    grid = scaling(runs, "E1_grid", f"Communication delay vs 2-D grid size, local machine ({v}, median and IQR)", "L1_grid_scaling")
    line = scaling(runs, "E2_line", f"Communication delay vs line topology 1 x N, local machine ({v})", "L2_line_scaling")
    if grid and line:
        compare_grid_line(grid[0], line[0], "L3_grid_vs_line")
    versions(runs, "E3_versions", "L4_versions", "L4_versions")
    param_sweep(runs, "E4_window", "window_s", "window W (s)", "L5_window")
    param_sweep(runs, "E5_poll", "poll_us", "polling interval (us, log scale)", "L6_poll", logx=True)
    param_sweep(runs, "E6_load", "p_arrive", "EV arrival probability per port per cycle", "L7_load")
    recheck(runs, "E7_recheck", "L8_recheck", "L8_recheck", "5 x 6")
    fits["local grid"], fits["local line"] = grid, line
    # CAAS set (C*): 5 grid sizes, same node counts as lines, across 1-5 hosts
    v = best_label(runs, "C1_grid")
    cg = scaling(runs, "C1_grid", f"CAAS: communication delay vs grid size ({v}, 16 cores per host)", "C1_grid_scaling", hosts=True)
    cl = scaling(runs, "C2_line", f"CAAS: communication delay vs line size 1 x N ({v})", "C2_line_scaling", hosts=True)
    if cg and cl:
        compare_grid_line(cg[0], cl[0], "C3_grid_vs_line", hosts=True)
    versions(runs, "C3_versions", "C4_versions", "C4_versions")
    recheck(runs, "C4_recheck", "C5_recheck", "C5_recheck", "6 x 6", r1_from="C3_versions")
    fits["caas grid"], fits["caas line"] = cg, cl
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "scaling_fits.txt"), "w") as f:
        for name, res in fits.items():
            if res:
                for key, (b, r2) in res[1].items():
                    f.write(f"{name:10s} {key:9s} slope {b * 1e3:+.2f} us/node  R2 {r2:.3f}\n")
    print(f"figures -> {FIG}\ntables  -> {TAB}")


if __name__ == "__main__":
    main()
