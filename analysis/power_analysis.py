# -*- coding: utf-8 -*-
"""
Post-hoc 统计功效分析（Monte Carlo 模拟）
口径与论文定稿脚本 14_结果统一重算.py 完全一致：
  - 路径系数 (a, b, c') 由数据用【带截距】OLS 取得（即中介分析的观测系数）；
  - 用各组观测系数生成数据，对每组重复 MC 次模拟，每次 Bootstrap 中介检验，
    统计显著率 = power；
  - gen_data 使用【模块级单一 rng】（仅 seed 一次 20250825），随机流跨场景连续，
    与 14 脚本逐字节一致，可复现稿件正文：参与组 48.5% / 未参与组 12.0% / 全样本 52.5%。
  - 同时扫描最小可检测效应量（80% power 对应的标准化间接效应）。
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(HERE, "..", "data", "deidentified_dataset.csv"))
df[["K_z", "A_z", "P_z"]] = StandardScaler().fit_transform(df[["k_score", "a_score", "p_score"]])


def coef_se(sub, x, y):
    X = np.column_stack([np.ones(len(sub))] + [sub[c].to_numpy() for c in x])
    yv = sub[y].to_numpy()
    n, k = X.shape
    beta, *_ = np.linalg.lstsq(X, yv, rcond=None)
    resid = yv - X @ beta
    sigma2 = (resid ** 2).sum() / (n - k)
    cov = sigma2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    return {c: beta[i + 1] for i, c in enumerate(x)}, {c: se[i + 1] for i, c in enumerate(x)}


def obs_paths(sub):
    K, A = sub["K_z"].to_numpy(), sub["A_z"].to_numpy()
    a = stats.linregress(K, A).slope
    b = coef_se(sub, ["A_z", "K_z"], "P_z")[0]["A_z"]      # A -> P | K
    cp = coef_se(sub, ["K_z", "A_z"], "P_z")[0]["K_z"]     # K -> P | A
    return a, b, cp


# 模块级单一 rng（仅 seed 一次），随机流跨所有 MC 与场景连续 —— 与 14 脚本一致
rng = np.random.default_rng(20250825)


def gen_data(n, a, b, cp):
    K = rng.normal(0, 1, n)
    A = a * K + rng.normal(0, np.sqrt(max(1 - a ** 2, 1e-9)), n)
    P = cp * K + b * A + rng.normal(0, np.sqrt(max(1 - cp ** 2 - b ** 2 - 2 * cp * b * a, 1e-9)), n)
    return K, A, P


def boot_indirect(K, A, P, n_boot=1000, seed=None):
    rr = np.random.default_rng(seed)
    n = len(K)
    inds = np.empty(n_boot)
    for i in range(n_boot):
        idx = rr.integers(0, n, n)
        Kb, Ab, Pb = K[idx], A[idx], P[idx]
        ab = np.cov(Kb, Ab, ddof=1)[0, 1] / np.var(Kb, ddof=1)
        # 间接效应 = a × b，其中 b 为 A→P|K 的偏回归系数（[Ab, Kb, ones] 的第 1 列）
        Xb = np.column_stack([Ab, Kb, np.ones(n)])
        bb = np.linalg.lstsq(Xb, Pb, rcond=None)[0][0]
        inds[i] = ab * bb
    return np.percentile(inds, [2.5, 97.5])


def power_at(n, a, b, cp, n_mc=200, n_boot=1000):
    sig = 0
    for m in range(n_mc):
        K, A, P = gen_data(n, a, b, cp)
        ci = boot_indirect(K, A, P, n_boot=n_boot, seed=m)
        if not (ci[0] <= 0 <= ci[1]):
            sig += 1
    return sig / n_mc


scenes = {
    "Full sample": df,
    "Participants": df[df["participation_status"] == 1],
    "Non-participants": df[df["participation_status"] == 0],
}

print("=" * 78)
print("Post-hoc power (Monte Carlo, 200 sims x Bootstrap 1000; rng seed 20250825)")
print("=" * 78)
for name, sub in scenes.items():
    a, b, cp = obs_paths(sub)
    pw = power_at(len(sub), a, b, cp)
    print(f"\n[{name}] n={len(sub)}  a={a:.3f} b={b:.3f} c'={cp:.3f}")
    print(f"  observed indirect a*b = {a*b:.4f}")
    print(f"  post-hoc power        = {pw*100:.1f}%")

print("\n" + "=" * 78)
print("Minimum detectable effect (80% power, alpha=0.05, two-sided)")
print("=" * 78)
for label, sub in [("Full sample n=68", df), ("Participants n=39", df[df["participation_status"] == 1])]:
    a0, _, _ = obs_paths(sub)
    best = None
    for b in np.arange(0.10, 1.01, 0.10):
        pw = power_at(len(sub), a0, b, 0.0, n_mc=100, n_boot=1000)
        if pw >= 0.80 and best is None:
            best = (b, a0 * b, pw)
    if best:
        print(f"  {label}: min detectable indirect effect ≈ {best[1]:.3f} (b={best[0]:.2f}, power={best[2]*100:.0f}%)")
    else:
        print(f"  {label}: 80% power not reached within scanned range")
print()
