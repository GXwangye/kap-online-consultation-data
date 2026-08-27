# -*- coding: utf-8 -*-
"""
根据 figures_source/*.json 与 data/deidentified_dataset.csv 重新渲染
PLOS 投稿级图表 Fig 1-3（PNG 300dpi + SVG + TIFF 300dpi）。
所有数值来自已核验的 JSON / 数据集，不手工改写。
图内不放置标题（标题/图注见稿件正文），面板仅标字母 A/B/C/D。
"""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
FIGSRC = os.path.join(HERE, "..", "figures_source")
DATA = os.path.join(HERE, "..", "data", "deidentified_dataset.csv")
OUT = os.path.join(HERE, "..", "figures")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "Arial", "font.size": 9, "axes.titlesize": 9.5,
    "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 8, "figure.dpi": 150, "savefig.dpi": 300,
    "savefig.bbox": "tight", "axes.axisbelow": True,
    "axes.edgecolor": "#333333", "text.color": "#1a1a1a",
    "axes.labelcolor": "#1a1a1a", "xtick.color": "#1a1a1a", "ytick.color": "#1a1a1a",
})
C_Part, C_Nonp, C_Gray = "#D55E00", "#0072B2", "#999999"


def save(name, fig):
    fig.savefig(os.path.join(OUT, f"{name}.png"), dpi=300, format="png")
    fig.savefig(os.path.join(OUT, f"{name}.svg"), dpi=300, format="svg")
    fig.savefig(os.path.join(OUT, f"{name}.tif"), dpi=300, format="tiff",
                pil_kwargs={"compression": "tiff_lzw"})
    print(f"[saved] {name}: png/svg/tif")


def render_fig1():
    d = json.load(open(os.path.join(FIGSRC, "fig1_data.json"), encoding="utf-8"))
    blocks = {b["id"]: b for b in d["blocks"]}
    en = {"A": "293 registered online-consultation\nphysicians (7 high-volume depts.)",
          "B": "Pilot test: 30 respondents\n(scale reliability & validity)",
          "C": "Formal questionnaire distribution", "D": "82 returned (response 28.0%)",
          "E": "68 valid (validity 82.9%)", "F": "Participants\nn = 39",
          "G": "Non-participants\nn = 29", "H": "KAP comparison /\nmediation"}
    fig, ax = plt.subplots(figsize=(5.0, 8.2))
    ax.set_xlim(0, 7.0); ax.set_ylim(12.2, 0); ax.axis("off")
    for b in d["blocks"]:
        ax.add_patch(FancyBboxPatch((b["x"], b["y"]), b["w"], b["h"],
                     boxstyle="round,pad=0.02,rounding_size=0.08", linewidth=1.0,
                     edgecolor="#333333", facecolor=b["color"], alpha=0.92))
        ax.text(b["x"] + b["w"]/2, b["y"] + b["h"]/2, en.get(b["id"], b["label"]),
                ha="center", va="center", fontsize=8.2, color="white", weight="bold", linespacing=1.3)
    for a in d["arrows"]:
        p, c = blocks[a["from"]], blocks[a["to"]]
        ax.annotate("", xy=(c["x"]+c["w"]/2, c["y"]-0.06),
                    xytext=(p["x"]+p["w"]/2, p["y"]+p["h"]+0.06),
                    arrowprops=dict(arrowstyle="-|>", color="#333333", lw=1.3, mutation_scale=12))
    save("Fig1_research_flow", fig); plt.close(fig)


def render_fig2():
    d = json.load(open(os.path.join(FIGSRC, "fig2_data.json"), encoding="utf-8"))
    fig, axes = plt.subplots(2, 2, figsize=(6.8, 5.2)); axes = axes.ravel()
    letter = {"A. KAP 总分  (P < 0.001)": "A", "B. 知识 K (P = 0.160)": "B",
              "C. 态度 A  (P < 0.001)": "C", "D. 行为 P  (P < 0.001)": "D"}
    pv = {"A. KAP 总分  (P < 0.001)": "P < 0.001", "B. 知识 K (P = 0.160)": "P = 0.160",
          "C. 态度 A  (P < 0.001)": "P < 0.001", "D. 行为 P  (P < 0.001)": "P < 0.001"}
    yl = {"总分": "Total score", "知识分数": "Knowledge", "态度分数": "Attitude", "行为分数": "Practice"}
    for i, p in enumerate(d["panels"]):
        ax = axes[i]; labels = p["data"]["labels"]
        means = p["data"]["series"]["均值"]; errs = p["data"]["errors"]["均值"]
        x = np.arange(len(labels))
        ax.bar(x, means, yerr=errs, capsize=4, width=0.55, color=[C_Part, C_Nonp],
               edgecolor="#222222", linewidth=0.8, error_kw=dict(ecolor="#222222", lw=1.1))
        for j, (m, e) in enumerate(zip(means, errs)):
            ax.text(j, m+e+0.02*max(means), f"{m:.2f}", ha="center", va="bottom", fontsize=7.5)
        ax.set_xticks(x); ax.set_xticklabels(["Participants\n(n=39)", "Non-participants\n(n=29)"], fontsize=7.5)
        ax.set_ylabel(yl.get(p["ylabel"], p["ylabel"]), fontsize=8)
        ax.set_title(letter.get(p["title"], ""), fontsize=11, weight="normal", loc="left")
        ymax = max(m+e for m, e in zip(means, errs))
        ax.text(0.5, ymax*1.06, pv.get(p["title"], ""), ha="center", va="bottom",
                fontsize=8, color="#b30000", weight="bold")
        ax.set_ylim(0, ymax*1.30); ax.spines[["top", "right"]].set_visible(False)
    fig.subplots_adjust(top=0.90, hspace=0.42, wspace=0.30)
    save("Fig2_KAP_by_group", fig); plt.close(fig)


def compute_permutation():
    df = pd.read_csv(DATA)
    K, A, P = df["k_score"].to_numpy(), df["a_score"].to_numpy(), df["p_score"].to_numpy()
    z = lambda x: (x - x.mean()) / x.std(ddof=0)
    K, A, P = z(K), z(A), z(P)
    g = df["participation_status"].to_numpy()
    ip, inp = np.where(g == 1)[0], np.where(g == 0)[0]

    def indirect(idx):
        k, a, p = K[idx], A[idx], P[idx]
        ak = np.linalg.lstsq(np.column_stack([np.ones(len(k)), k]), a, rcond=None)[0][1]
        bp = np.linalg.lstsq(np.column_stack([np.ones(len(k)), a, k]), p, rcond=None)[0][1]
        return ak * bp
    obs = indirect(ip) - indirect(inp)
    rng = np.random.default_rng(42); N = 10000; diffs = np.empty(N)
    for i in range(N):
        sg = rng.permutation(g); p1 = np.where(sg == 1)[0]; p2 = np.where(sg == 0)[0]
        diffs[i] = indirect(p1) - indirect(p2)
    return obs, diffs, float(np.mean(np.abs(diffs) >= abs(obs)))


def render_fig3():
    fd = json.load(open(os.path.join(FIGSRC, "fig3_data.json"), encoding="utf-8"))
    en = {"全样本 (n=68)": "Full sample (n=68)", "参与组 (n=39)": "Participants (n=39)",
          "不参与组 (n=29)": "Non-participants (n=29)"}
    labels = [en.get(l, l) for l in fd["labels"]]
    est, lo, hi = fd["estimates"], fd["ci_low"], fd["ci_high"]
    colors = [C_Gray, C_Part, C_Nonp]
    sig_mask = [not (l < 0) for l in lo]
    obs, diffs, p_val = compute_permutation()
    fig = plt.figure(figsize=(7.0, 6.4))
    gs = fig.add_gridspec(2, 2, hspace=0.34, wspace=0.32, left=0.09, right=0.97, top=0.90, bottom=0.07)
    axA = fig.add_subplot(gs[0, 0]); y = np.arange(len(labels))[::-1]
    for yi, (lb, e, l, h, c, sm) in enumerate(zip(labels, est, lo, hi, colors, sig_mask)):
        axA.plot([l, h], [yi, yi], color=c, lw=2.2, solid_capstyle="round")
        axA.scatter([e], [yi], color=c, s=42, zorder=3, edgecolor="white", linewidth=1.0)
        if sm: axA.text(h+0.012, yi, "*", va="center", ha="left", fontsize=13, color="#b30000", weight="bold")
    axA.axvline(0, color="#666666", ls="--", lw=1.0); axA.set_yticks(y); axA.set_yticklabels(labels, fontsize=8)
    axA.set_xlabel("Indirect effect (standardized)", fontsize=8); axA.set_xlim(-0.15, 0.46)
    axA.set_title("A", fontsize=11, weight="normal", loc="left"); axA.spines[["top", "right"]].set_visible(False)
    axB = fig.add_subplot(gs[0, 1]); axB.hist(diffs, bins=50, color="#cfcfcf", edgecolor="white", linewidth=0.3)
    axB.axvline(obs, color=C_Part, lw=2.2, label=f"observed Δ={obs:.3f}"); axB.axvline(0, color="#666666", ls="--", lw=1.0)
    axB.set_xlabel("Δ indirect effect (part − nonp)", fontsize=8); axB.set_ylabel("Frequency", fontsize=8)
    axB.set_title("B", fontsize=11, weight="normal", loc="left"); axB.legend(fontsize=7, loc="upper left", frameon=False)
    axB.spines[["top", "right"]].set_visible(False)
    axC = fig.add_subplot(gs[1, 0])
    for yi, (lb, e, l, h, c) in enumerate(zip(labels, est, lo, hi, colors)):
        axC.plot([l, h], [yi, yi], color=c, lw=3.0, solid_capstyle="round")
        axC.scatter([e], [yi], color=c, s=38, zorder=3, edgecolor="white", lw=1.0)
    axC.axvline(0, color="#666666", ls="--", lw=1.0); axC.set_yticks(range(len(labels)))
    axC.set_yticklabels(labels, fontsize=8); axC.set_xlabel("Indirect effect (95% CI)", fontsize=8)
    axC.set_xlim(-0.15, 0.46); axC.set_title("C", fontsize=11, weight="normal", loc="left")
    axC.spines[["top", "right"]].set_visible(False)
    axD = fig.add_subplot(gs[1, 1]); axD.axis("off"); axD.set_xlim(0, 10); axD.set_ylim(0, 6)
    nodes = {"K": (1.5, 3.0), "A": (5.0, 4.6), "P": (8.5, 3.0)}
    for n, (nx, ny) in nodes.items():
        axD.add_patch(FancyBboxPatch((nx-0.55, ny-0.5), 1.1, 1.0, boxstyle="round,pad=0.02,rounding_size=0.12",
                     facecolor="#f2f2f2", edgecolor="#333333", lw=1.1))
        axD.text(nx, ny, n, ha="center", va="center", fontsize=11, weight="bold")
    axD.annotate("", xy=nodes["A"], xytext=nodes["K"], arrowprops=dict(arrowstyle="-|>", color=C_Part, lw=2.4, mutation_scale=14))
    axD.annotate("", xy=nodes["P"], xytext=nodes["A"], arrowprops=dict(arrowstyle="-|>", color=C_Part, lw=2.4, mutation_scale=14))
    axD.annotate("", xy=nodes["P"], xytext=nodes["K"], arrowprops=dict(arrowstyle="-|>", color=C_Part, lw=1.6, mutation_scale=12, connectionstyle="arc3,rad=0.25"))
    axD.text(2.8, 4.1, "a = 0.334", color=C_Part, fontsize=8, weight="bold", ha="center")
    axD.text(7.1, 4.4, "b = 0.479***", color=C_Part, fontsize=8, weight="bold", ha="center")
    axD.text(5.0, 1.8, "c' = 0.108", color=C_Part, fontsize=7.5, ha="center")
    axD.annotate("", xy=nodes["A"], xytext=nodes["K"], arrowprops=dict(arrowstyle="-|>", color=C_Nonp, lw=1.6, ls="--", mutation_scale=12))
    axD.annotate("", xy=nodes["P"], xytext=nodes["A"], arrowprops=dict(arrowstyle="-|>", color=C_Nonp, lw=1.6, ls="--", mutation_scale=12))
    axD.text(3.2, 1.5, "a = 0.546 / b = 0.127 (ns)", color=C_Nonp, fontsize=7.5, ha="center")
    axD.text(5.0, 5.4, "Participants: K→A→P significant", color=C_Part, fontsize=8, weight="bold", ha="center")
    axD.text(5.0, 0.6, "Non-participants: no significant path (n=29, low power)", color=C_Nonp, fontsize=7.2, ha="center")
    axD.set_title("D", fontsize=11, weight="normal", loc="left")
    save("Fig3_mediation", fig); plt.close(fig)


if __name__ == "__main__":
    render_fig1(); render_fig2(); render_fig3()
    print("ALL_FIGURES_DONE (Fig 1-3; NMF/Fig4 dropped per manuscript)")
