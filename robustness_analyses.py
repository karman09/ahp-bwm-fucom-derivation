# -*- coding: utf-8 -*-
"""
Robustness analyses for the derivation of BWM and FUCOM weights from AHP matrices.

Companion to the main script (MAIN_SCRIPT below), from which the MCDM_Extractor class,
the agreement measures and the seven published cases are loaded.

Analyses
  A. Input-based consistency ratio (Liang, Brunelli & Rezaei, 2020) of the extracted
     Best-to-Others / Others-to-Worst vectors, with the corresponding thresholds.
  B. Sensitivity to the choice of the best and worst criteria: each case is re-solved
     with the second-best criterion as B or the second-worst criterion as W.
  C. Monte Carlo simulation: random AHP matrices of known weights are perturbed,
     rounded to the Saaty scale, and the derived BWM (linear) and FUCOM (closed form,
     Eq. 8) weights are compared with the AHP weights and with the generating weights.

Outputs
  Console report, figure_simulation.svg/.pdf/.png and Robustness_Report.xlsx
  (sheets: InputCR, BW_Sensitivity, Simulation_Summary, Simulation_Raw, Simulation_Ties).
"""
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import rankdata

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------- #
# 0) Load the extractor, the agreement measures and the cases from the main script.
#    Only the part before the console report is executed, so that the main script's
#    figures and Excel report are not regenerated.
# --------------------------------------------------------------------------- #
MAIN_SCRIPT = "bwm_fucom_derivation.py"   # file name of the main script (if not run beforehand)
MARKER = "# 2) CONSOLE REPORT"

if all(k in globals() for k in ("MCDM_Extractor", "all_metrics", "EXAMPLES")):
    # Colab / Jupyter: the main script has already been run in an earlier cell.
    pass
else:
    import os
    if not os.path.exists(MAIN_SCRIPT):
        raise FileNotFoundError(
            f"'{MAIN_SCRIPT}' not found. Either run the main script in an earlier cell, "
            "or upload it to the working directory and set MAIN_SCRIPT to its file name.")
    with open(MAIN_SCRIPT, encoding="utf-8") as fh:
        _src = fh.read()
    _ns = {"__name__": "main_part"}
    exec(compile(_src.split(MARKER)[0], MAIN_SCRIPT, "exec"), _ns)
    MCDM_Extractor = _ns["MCDM_Extractor"]
    all_metrics = _ns["all_metrics"]
    EXAMPLES = _ns["EXAMPLES"]

SEED = 2026
N_REP = 500                                   # matrices per (n, sigma) cell
N_VALUES = [4, 5, 6, 7, 8, 9]
SIGMAS = [0.05, 0.15, 0.25, 0.35, 0.50, 0.70]  # s.d. of the log-normal perturbation
CR_BINS = [(0.00, 0.05), (0.05, 0.10), (0.10, 0.20), (0.20, np.inf)]
CR_LABELS = ["< 0.05", "0.05–0.10", "0.10–0.20", "≥ 0.20"]


# ============================================================================ #
# A) Input-based consistency ratio (Liang et al., 2020)
# ============================================================================ #
# Thresholds (Liang et al., 2020, Table 3): rows = scale a_BW (3..9), columns = n (3..9).
CRI_THRESHOLDS = {
    3: [0.1667, 0.1667, 0.1667, 0.1667, 0.1667, 0.1667, 0.1667],
    4: [0.1121, 0.1529, 0.1898, 0.2206, 0.2527, 0.2577, 0.2683],
    5: [0.1354, 0.1994, 0.2306, 0.2546, 0.2716, 0.2844, 0.2960],
    6: [0.1330, 0.1990, 0.2643, 0.3044, 0.3144, 0.3221, 0.3262],
    7: [0.1294, 0.2457, 0.2819, 0.3029, 0.3144, 0.3251, 0.3403],
    8: [0.1309, 0.2521, 0.2958, 0.3154, 0.3408, 0.3620, 0.3657],
    9: [0.1359, 0.2681, 0.3062, 0.3337, 0.3517, 0.3620, 0.3662],
}


def input_based_cr(bo, ow, b_idx, w_idx):
    """CR^I = max_j |a_Bj a_jW - a_BW| / (a_BW (a_BW - 1)); CR^I = 0 if a_BW = 1."""
    a_bw = bo[w_idx]
    if abs(a_bw - 1.0) < 1e-9:
        return 0.0
    return float(np.max(np.abs(bo * ow - a_bw)) / (a_bw * (a_bw - 1.0)))


def cri_threshold(a_bw, n):
    """Threshold for the scale ceil(a_BW) (non-integer a_BW arise from aggregated matrices)."""
    scale = int(np.ceil(a_bw - 1e-9))
    if scale < 3 or not 3 <= n <= 9:
        return np.nan
    return CRI_THRESHOLDS[min(scale, 9)][n - 3]


rows_cri = []
for k, ex in enumerate(EXAMPLES, 1):
    b = ex.names.index(ex.best_criteria[0])
    w = ex.names.index(ex.worst_criteria[0])
    cri = input_based_cr(ex.bo_vector, ex.ow_vector, b, w)
    a_bw = float(ex.bo_vector[w])
    thr = cri_threshold(a_bw, ex.n)
    rows_cri.append({"Case": k, "n": ex.n, "B": ex.names[b], "W": ex.names[w], "a_BW": a_bw,
                     "CR_I": cri, "Threshold": thr,
                     "Acceptable": (cri <= thr) if not np.isnan(thr) else None,
                     "xi_L": ex.bwm_xi_linear})
df_cri = pd.DataFrame(rows_cri)


# ============================================================================ #
# B) Sensitivity to the choice of the best and worst criteria
# ============================================================================ #
def solve_with(ex, b_idx, w_idx):
    """Re-solve the extracted BWM and FUCOM models with a prescribed (B, W)."""
    e = MCDM_Extractor(ex.matrix, ex.names, ex.reported_weights, bwm_model=ex.bwm_model)
    e.calculate_ahp()
    e.best_criteria, e.worst_criteria = [ex.names[b_idx]], [ex.names[w_idx]]
    e.bo_vector, e.ow_vector = e.matrix[b_idx, :], e.matrix[:, w_idx]
    e.assumption_ok = bool(np.all(e.bo_vector >= 1 - 1e-9) and np.all(e.ow_vector >= 1 - 1e-9))
    e.solve_bwm()
    e.extract_fucom()
    # FUCOM weights from the closed form (Eq. 8, Proposition 1), as in the simulation, so that
    # the result does not depend on the convergence of the numerical solver when a_Bj < 1.
    e.fucom_weights = (1.0 / e.bo_vector) / np.sum(1.0 / e.bo_vector)
    return e


rows_sens = []
for k, ex in enumerate(EXAMPLES, 1):
    w_ahp = ex.ahp_weights
    order = np.argsort(-w_ahp, kind="stable")          # descending, ties in index order
    # Baseline B and W as in the main analysis (Section 3.2: lowest index among tied criteria)
    b0 = ex.names.index(ex.best_criteria[0])
    w0 = ex.names.index(ex.worst_criteria[0])
    # Second-best / second-worst: the criterion with the next largest / smallest AHP weight
    b1 = next(i for i in order if i != b0)
    w1 = next(i for i in order[::-1] if i != w0)
    variants = [("Baseline", b0, w0, np.nan),
                ("Second-worst as W", b0, w1, abs(w_ahp[w1] - w_ahp[w0])),
                ("Second-best as B", b1, w0, abs(w_ahp[b0] - w_ahp[b1]))]
    for tag, bi, wi, gap in variants:
        e = solve_with(ex, bi, wi)
        rec = {"Case": k, "Selection": tag, "B": ex.names[bi], "W": ex.names[wi], "Gap": gap,
               "Assumption_ok": e.assumption_ok,
               "xi_L": e.bwm_xi_linear,
               "r_BWM_AHP": float(np.corrcoef(e.bwm_weights, w_ahp)[0, 1]),
               "r_FUCOM_AHP": float(np.corrcoef(e.fucom_weights, w_ahp)[0, 1])}
        for m, wv in (("BWM", e.bwm_weights), ("FUCOM", e.fucom_weights)):
            if m in ex.reported_weights:
                mt = all_metrics(wv, ex.reported_weights[m])
                rec[f"r_{m}_rep"], rec[f"WS_{m}_rep"] = mt["r"], mt["ws"]
            else:
                rec[f"r_{m}_rep"] = rec[f"WS_{m}_rep"] = np.nan
        rows_sens.append(rec)
df_sens = pd.DataFrame(rows_sens)


# ============================================================================ #
# C) Monte Carlo simulation
# ============================================================================ #
SAATY = np.array([1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3, 1/2, 1, 2, 3, 4, 5, 6, 7, 8, 9])
LOG_SAATY = np.log(SAATY)


def random_matrix(rng, n, sigma):
    """Generating weights ~ Dirichlet(1); a_ij = (w_i/w_j) exp(e_ij), e_ij ~ N(0, sigma^2),
    rounded to the nearest Saaty value on the log scale; the lower triangle is reciprocal."""
    w_true = rng.dirichlet(np.ones(n))
    A = np.ones((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            v = np.log(w_true[i] / w_true[j]) + rng.normal(0.0, sigma)
            a = SAATY[np.argmin(np.abs(LOG_SAATY - v))]
            A[i, j], A[j, i] = a, 1.0 / a
    return A, w_true


def safe_r(a, b):
    """Pearson r; equal to 1 when both vectors are constant (all weights equal)."""
    sa, sb = np.std(a), np.std(b)
    if sa < 1e-12 and sb < 1e-12:
        return 1.0
    if sa < 1e-12 or sb < 1e-12:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def derive(A, b_idx, w_idx):
    """Linear BWM and closed-form FUCOM (Eq. 8) weights for a prescribed (B, W)."""
    e = MCDM_Extractor(A)
    e.calculate_ahp()
    e.extract_bwm()
    e.best_criteria, e.worst_criteria = [e.names[b_idx]], [e.names[w_idx]]
    e.bo_vector, e.ow_vector = e.matrix[b_idx, :], e.matrix[:, w_idx]
    w_bwm, _ = e.solve_bwm_linear()
    w_fuc = (1.0 / e.bo_vector) / np.sum(1.0 / e.bo_vector)
    return w_bwm, w_fuc


def ws_coef(calc, ref, decimals=4):
    """WS coefficient (Eq. 11) with `ref` as the reference ranking."""
    Rd = rankdata(-np.round(calc, decimals), method="average")
    Rr = rankdata(-np.round(ref, decimals), method="average")
    n = len(Rr)
    return 1.0 - sum(2.0 ** (-Rr[j]) * abs(Rr[j] - Rd[j]) / max(abs(1 - Rr[j]), abs(n - Rr[j]))
                     for j in range(n))


rng = np.random.default_rng(SEED)
sim, tie_rows = [], []
for n in N_VALUES:
    for sigma in SIGMAS:
        for _ in range(N_REP):
            A, w_true = random_matrix(rng, n, sigma)
            e = MCDM_Extractor(A)
            e.calculate_ahp()
            e.extract_bwm()
            w_bwm, xi = e.solve_bwm_linear()
            w_fuc = (1.0 / e.bo_vector) / np.sum(1.0 / e.bo_vector)   # Eq. (8), Proposition 1
            w_ahp = e.ahp_weights
            b_idx = e.names.index(e.best_criteria[0])
            w_idx = e.names.index(e.worst_criteria[0])
            a_bw = float(e.bo_vector[w_idx])
            cri = input_based_cr(e.bo_vector, e.ow_vector, b_idx, w_idx)
            thr = cri_threshold(a_bw, n) if a_bw >= 1 - 1e-9 else np.nan
            sim.append({
                "n": n, "sigma": sigma, "CR": e.ahp_cr, "xi_L": xi, "CR_I": cri, "CR_I_thr": thr,
                "CR_I_ok": (cri <= thr) if not np.isnan(thr) else np.nan,
                "tie_B": len(e.best_criteria) > 1, "tie_W": len(e.worst_criteria) > 1,
                "assumption_ok": e.assumption_ok,
                "r_BWM": safe_r(w_bwm, w_ahp), "r_FUCOM": safe_r(w_fuc, w_ahp),
                "WS_BWM": ws_coef(w_bwm, w_ahp), "WS_FUCOM": ws_coef(w_fuc, w_ahp),
                "MAE_BWM": np.mean(np.abs(w_bwm - w_ahp)), "MAE_FUCOM": np.mean(np.abs(w_fuc - w_ahp)),
                "MAEtrue_AHP": np.mean(np.abs(w_ahp - w_true)),
                "MAEtrue_BWM": np.mean(np.abs(w_bwm - w_true)),
                "MAEtrue_FUCOM": np.mean(np.abs(w_fuc - w_true)),
            })
            # Ties: re-derive the weights with every admissible (B, W) pair among the tied criteria
            if len(e.best_criteria) > 1 or len(e.worst_criteria) > 1:
                bs = [e.names.index(c) for c in e.best_criteria]
                ws_ = [e.names.index(c) for c in e.worst_criteria]
                sols = [derive(A, bi, wi) for bi in bs for wi in ws_ if bi != wi]
                if len(sols) > 1:
                    dev_bwm = max(np.max(np.abs(s1[0] - s2[0])) for s1 in sols for s2 in sols)
                    dev_fuc = max(np.max(np.abs(s1[1] - s2[1])) for s1 in sols for s2 in sols)
                    tie_rows.append({"n": n, "CR": e.ahp_cr, "n_B": len(bs), "n_W": len(ws_),
                                     "maxdev_BWM": dev_bwm, "maxdev_FUCOM": dev_fuc})
df_sim = pd.DataFrame(sim)
df_sim["CR_bin"] = pd.cut(df_sim["CR"], [b[0] for b in CR_BINS] + [np.inf],
                          right=False, labels=CR_LABELS)


def p10(x):
    return np.nanpercentile(x, 10)


def summarise(df, by):
    g = df.groupby(by, observed=True)
    out = g.agg(N=("r_BWM", "size"),
                r_BWM_med=("r_BWM", "median"), r_BWM_p10=("r_BWM", p10),
                r_FUCOM_med=("r_FUCOM", "median"), r_FUCOM_p10=("r_FUCOM", p10),
                WS_BWM_med=("WS_BWM", "median"), WS_FUCOM_med=("WS_FUCOM", "median"),
                MAE_BWM_med=("MAE_BWM", "median"), MAE_FUCOM_med=("MAE_FUCOM", "median"),
                MAEtrue_AHP_med=("MAEtrue_AHP", "median"),
                MAEtrue_BWM_med=("MAEtrue_BWM", "median"),
                MAEtrue_FUCOM_med=("MAEtrue_FUCOM", "median"),
                assumption_viol=("assumption_ok", lambda s: 1 - s.mean()),
                CRI_exceed=("CR_I_ok", lambda s: 1 - np.nanmean(s.astype(float))))
    return out.reset_index()


sum_bin = summarise(df_sim, ["CR_bin"])
sum_n_bin = summarise(df_sim, ["n", "CR_bin"])
sum_acc = summarise(df_sim[df_sim["CR"] < 0.10], ["n"])
df_tie = pd.DataFrame(tie_rows)

# --------------------------------------------------------------------------- #
# Console report
# --------------------------------------------------------------------------- #
pd.set_option("display.width", 160)
pd.set_option("display.float_format", lambda v: f"{v:0.4f}")
print("=" * 78, "\n A) Input-based consistency ratio of the extracted vectors\n", "=" * 78)
print(df_cri.to_string(index=False))
print("\n", "=" * 78, "\n B) Sensitivity to the choice of best and worst criteria\n", "=" * 78)
print(df_sens.to_string(index=False))
print("\n", "=" * 78, f"\n C) Simulation: {len(df_sim)} matrices "
      f"(n = {N_VALUES[0]}..{N_VALUES[-1]}, {N_REP} per n x sigma cell, seed {SEED})\n", "=" * 78)
print(" Exact ties in B: %d, in W: %d, in B or W: %d of %d matrices"
      % (df_sim.tie_B.sum(), df_sim.tie_W.sum(), (df_sim.tie_B | df_sim.tie_W).sum(), len(df_sim)))
if len(df_tie):
    print(" Tied matrices: max |dw| across admissible (B, W) choices  "
          "BWM median = %.4f, max = %.4f;  FUCOM median = %.4f, max = %.4f;  share with BWM |dw| < 1e-6 = %.3f"
          % (df_tie.maxdev_BWM.median(), df_tie.maxdev_BWM.max(), df_tie.maxdev_FUCOM.median(),
             df_tie.maxdev_FUCOM.max(), (df_tie.maxdev_BWM < 1e-6).mean()))
    acc = df_tie[df_tie.CR < 0.10]
    print(" Tied matrices with CR < 0.10 (N = %d): BWM median = %.4f, p90 = %.4f; FUCOM median = %.4f, p90 = %.4f"
          % (len(acc), acc.maxdev_BWM.median(), acc.maxdev_BWM.quantile(0.9),
             acc.maxdev_FUCOM.median(), acc.maxdev_FUCOM.quantile(0.9)))
print("\n By CR class (all n):\n" + sum_bin.to_string(index=False))
print("\n By n and CR class:\n" + sum_n_bin.to_string(index=False))
print("\n Acceptable matrices (CR < 0.10), by n:\n" + sum_acc.to_string(index=False))

# ============================================================================ #
# Figure (drawn at its printed size, 15 cm wide, so that fonts appear at their
# nominal size in the manuscript; saved as vector SVG/PDF and as a 600-dpi PNG)
# ============================================================================ #
COLOR_BWM, COLOR_FUCOM = "#d9531e", "#128a5e"
COLOR_GRID, COLOR_TEXT, COLOR_AXIS = "#dcdcd8", "#222222", "#7a7a76"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.edgecolor": COLOR_AXIS, "axes.labelcolor": COLOR_TEXT, "text.color": COLOR_TEXT,
    "xtick.color": COLOR_TEXT, "ytick.color": COLOR_TEXT,
    "axes.labelsize": 8.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "svg.fonttype": "path", "pdf.fonttype": 42,
})

fig, axes = plt.subplots(1, 2, figsize=(5.9, 2.75), sharey=True)
panels = [(axes[0], sum_bin, CR_LABELS, "AHP consistency ratio (CR)", "(a) All matrices, by CR class"),
          (axes[1], sum_acc, N_VALUES, "Number of criteria (n)", "(b) Matrices with CR < 0.10, by n")]
for ax, d, ticks, xlabel, title in panels:
    x = np.arange(len(d))
    for m, col, mk, off in (("BWM", COLOR_BWM, "o", -0.11), ("FUCOM", COLOR_FUCOM, "s", 0.11)):
        med, lo = d[f"r_{m}_med"].values, d[f"r_{m}_p10"].values
        ax.vlines(x + off, lo, med, color=col, linewidth=1.4, alpha=0.5, zorder=2)
        ax.plot(x + off, med, color=col, linewidth=1.4, marker=mk, markersize=4.5,
                markeredgecolor="white", markeredgewidth=0.8, label=f"{m}: median", zorder=3)
        ax.scatter(x + off, lo, marker="_", s=45, linewidths=1.4, color=col, zorder=3,
                   label=f"{m}: 10th percentile")
    ax.set_xticks(x)
    ax.set_xticklabels([str(t) for t in ticks])
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=8.5, fontweight="bold", loc="left", pad=6)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.yaxis.grid(True, color=COLOR_GRID, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=0, pad=3)
    ax.set_xlim(-0.5, len(d) - 0.5)
axes[0].set_ylabel("Pearson r with AHP weights")
axes[0].set_ylim(0.4, 1.0)
axes[0].set_yticks(np.arange(0.4, 1.01, 0.1))
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, fontsize=7.5,
           handlelength=1.6, columnspacing=1.2, bbox_to_anchor=(0.5, 0.0))
fig.tight_layout(rect=(0, 0.08, 1, 1), w_pad=1.5)
for ext, kw in (("svg", {}), ("pdf", {}), ("png", {"dpi": 600})):
    fig.savefig(f"figure_simulation.{ext}", facecolor="white", **kw)
plt.close(fig)
print("\nFigure saved: figure_simulation.svg / .pdf / .png (600 dpi)")

# ============================================================================ #
# Excel report
# ============================================================================ #
with pd.ExcelWriter("Robustness_Report.xlsx", engine="openpyxl") as xw:
    df_cri.to_excel(xw, sheet_name="InputCR", index=False)
    df_sens.to_excel(xw, sheet_name="BW_Sensitivity", index=False)
    sum_bin.to_excel(xw, sheet_name="Simulation_Summary", index=False, startrow=1)
    r2 = len(sum_bin) + 4
    sum_acc.to_excel(xw, sheet_name="Simulation_Summary", index=False, startrow=r2)
    r3 = r2 + len(sum_acc) + 3
    sum_n_bin.to_excel(xw, sheet_name="Simulation_Summary", index=False, startrow=r3)
    ws = xw.sheets["Simulation_Summary"]
    ws.cell(row=1, column=1, value="By CR class (all n)")
    ws.cell(row=r2, column=1, value="Acceptable matrices (CR < 0.10), by n")
    ws.cell(row=r3, column=1, value="By n and CR class")
    df_sim.to_excel(xw, sheet_name="Simulation_Raw", index=False)
    if len(df_tie):
        df_tie.to_excel(xw, sheet_name="Simulation_Ties", index=False)
print("Excel report saved: Robustness_Report.xlsx")