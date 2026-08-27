# -*- coding: utf-8 -*-
"""
KAP 中介效应分析（Hayes PROCESS Model 4, Bootstrap 5000 次）
数据：../data/deidentified_dataset.csv （n=68；参与组 39，未参与组 29）
方法（与论文定稿口径 14_结果统一重算.py 一致）：
      - 标准化 K/A/P 在【全样本】上完成，再按组切片；
      - 路径 a=K→A（linregress，含截距）；b=A→P|K、c'=K→P|A 用【带截距】多元 OLS；
      - 间接效应 = a×b，Bootstrap 5000 次 (seed=42)，百分位 95% CI。
输出应与稿件正文一致：
      参与组 间接效应 0.160 (95% CI 0.048–0.294)
      全样本 0.179 (95% CI −0.029–0.383)
      未参与组 0.069 (95% CI −0.094–0.262)
"""
import os
import pandas as pd
import numpy as np
from scipy import stats
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(HERE, "..", "data", "deidentified_dataset.csv"))

# 全样本标准化（与定稿脚本一致）
df[["K_z", "A_z", "P_z"]] = StandardScaler().fit_transform(df[["k_score", "a_score", "p_score"]])


def coef_se(sub, x, y):
    """带截距多元 OLS 的回归系数与标准误。"""
    X = np.column_stack([np.ones(len(sub))] + [sub[c].to_numpy() for c in x])
    yv = sub[y].to_numpy()
    n, k = X.shape
    beta, *_ = np.linalg.lstsq(X, yv, rcond=None)
    resid = yv - X @ beta
    sigma2 = (resid ** 2).sum() / (n - k)
    cov = sigma2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    return {c: beta[i + 1] for i, c in enumerate(x)}, {c: se[i + 1] for i, c in enumerate(x)}


def mediation(sub, n_boot=5000, seed=42):
    rng = np.random.default_rng(seed)
    n = len(sub)
    K, A, P = sub["K_z"].to_numpy(), sub["A_z"].to_numpy(), sub["P_z"].to_numpy()
    a = stats.linregress(K, A).slope
    b = coef_se(sub, ["A_z", "K_z"], "P_z")[0]["A_z"]      # A -> P | K
    cp = coef_se(sub, ["K_z", "A_z"], "P_z")[0]["K_z"]     # K -> P | A
    indirect = a * b
    inds = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        Kb, Ab, Pb = K[idx], A[idx], P[idx]
        ab = stats.linregress(Kb, Ab).slope
        Xb = np.column_stack([np.ones(n), Ab, Kb])
        bb = np.linalg.lstsq(Xb, Pb, rcond=None)[0][1]     # A_z 系数（含截距）
        inds[i] = ab * bb
    ci = np.percentile(inds, [2.5, 97.5])
    return dict(a=a, b=b, cp=cp, indirect=indirect, ci=ci, sig=not (ci[0] <= 0 <= ci[1]))


participants = df[df["participation_status"] == 1]
non_participants = df[df["participation_status"] == 0]
groups = {
    "Full sample (n=68)": df,
    "Participants (n=39)": participants,
    "Non-participants (n=29)": non_participants,
}

print("=" * 78)
print("K -> A -> P mediation (standardized, intercept-included OLS; Bootstrap 5000)")
print("=" * 78)
for name, sub in groups.items():
    r = mediation(sub)
    print(f"\n[{name}]")
    print(f"  a (K->A)        = {r['a']:.4f}")
    print(f"  b (A->P | K)    = {r['b']:.4f}")
    print(f"  c' (K->P | A)   = {r['cp']:.4f}")
    print(f"  indirect a*b    = {r['indirect']:.4f}")
    print(f"  95% CI          = [{r['ci'][0]:.4f}, {r['ci'][1]:.4f}]")
    print(f"  significant     = {r['sig']}")
print()
