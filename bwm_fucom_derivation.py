# -*- coding: utf-8 -*-
"""
Deriving BWM and FUCOM weights from an AHP pairwise comparison matrix and
comparing them with the weights reported in published studies.

Procedure
  1. Input: the pairwise comparison matrices are used as published; deviations
     from exact reciprocity (a_ji * a_ij != 1, due to rounding) are reported but
     not modified. The only correction is the transposed entry pair in Case 1.
  2. AHP: principal eigenvector and consistency ratio (Eq. 2).
  3. Best (B) and worst (W) criteria: largest and smallest AHP weights (Eq. 3).
  4. BWM: row B and column W (Eqs. 4-5); linear (Rezaei, 2016) and non-linear
     (Rezaei, 2015) models, the variant of the source study used for validation.
  5. FUCOM: ranking and comparative priorities taken from row B
     (Pamucar, Stevic & Sremac, 2018; Eqs. 6-8).
  6. Agreement: Pearson r, Spearman rho, WS coefficient (Salabun & Urbaniak,
     2020), MAE and maximum absolute deviation (Eqs. 9-12). Rank-based measures
     are computed from weights rounded to four decimals, with tied weights
     assigned average ranks.

Cases
  Case 1 - Resource allocation in cloud computing (Pamucar, Stevic & Sremac, 2018)
  Case 2 - Car purchase (Pamucar, Stevic & Sremac, 2018)
  Case 3 - Smart-city components (Ekin & Sarul, 2022)
  Case 4 - Beekeeping land suitability, Mardin (Mercan & Acibuca, 2025)
  Case 5 - Ranking of Libyan airlines (Badi & Abdulshahed, 2019)
  Case 6 - Wind-farm site selection, Libya (Badi et al., 2023)
  Case 7 - Flash-flood susceptibility, Ibar basin (Vujovic et al., 2025)

Outputs
  Console report, figure_case1-7.png, BWM_FUCOM_Derivation_Report.xlsx
  (sheets: Summary, Manuscript_Tables and one sheet per case) and the manuscript
  figures Figure_1-3 (.svg, .pdf, .png, .tiff).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.optimize import minimize, linprog
from scipy.stats import pearsonr, rankdata

from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from openpyxl.chart import BarChart, Reference, Series
from openpyxl.formatting.rule import CellIsRule, ColorScaleRule
from openpyxl.worksheet.page import PageMargins


# ============================================================================ #
# 0) MCDM_Extractor: derivation of BWM and FUCOM weights from an AHP matrix
# ============================================================================ #
class MCDM_Extractor:
    RI_dict = {1: 0.00, 2: 0.00, 3: 0.58, 4: 0.90, 5: 1.12,
               6: 1.24, 7: 1.32, 8: 1.41, 9: 1.45, 10: 1.49}

    def __init__(self, matrix, criteria_names=None, reported_weights=None, label="", reported_bw=None,
                 bwm_model="linear"):
        self.matrix = np.array(matrix, dtype=float)
        self.n = self.matrix.shape[0]
        self.names = criteria_names if criteria_names else [f"C{i+1}" for i in range(self.n)]
        self.reported_weights = reported_weights or {}
        self.label = label
        # reported_bw: best/worst selections reported in the source study,
        # given as [{"source": ..., "best": "C4", "worst": "C2"}, ...].
        self.reported_bw = reported_bw or []
        # bwm_model: BWM variant used for the comparison with the source study --
        # "linear" (Rezaei, 2016) or "nonlinear" (Rezaei, 2015). It should match the
        # variant used in the source; the other variant is reported as a sensitivity check.
        if bwm_model not in ("linear", "nonlinear"):
            raise ValueError("bwm_model must be 'linear' or 'nonlinear'")
        self.bwm_model = bwm_model
        self.bw_matches = []
        self.nonreciprocal_cells = []   # (j, i, a_ji, 1/a_ij) where a_ji != 1/a_ij

        # Input check: positive entries, unit diagonal; the matrix is used as published.
        if np.any(self.matrix <= 0):
            raise ValueError("All pairwise comparison entries must be positive.")
        if np.any(np.abs(np.diag(self.matrix) - 1.0) > 1e-9):
            raise ValueError("Diagonal entries must equal 1.")
        for i in range(self.n):
            for j in range(i + 1, self.n):
                expected_val = 1.0 / self.matrix[i, j]
                if abs(self.matrix[j, i] - expected_val) > 1e-5:
                    self.nonreciprocal_cells.append((j, i, float(self.matrix[j, i]), float(expected_val)))
        self.max_reciprocity_dev = float(np.max(np.abs(self.matrix * self.matrix.T - 1.0)))

        self.ahp_weights = None
        self.ahp_cr = None
        self.best_criteria, self.worst_criteria = [], []
        self.bo_vector = self.ow_vector = None
        self.fucom_seq = self.fucom_phi = None
        self.bwm_weights = self.fucom_weights = None

    def calculate_ahp(self):
        eigvals, eigvecs = np.linalg.eig(self.matrix)
        idx = np.argmax(np.real(eigvals))
        lambda_max = np.real(eigvals[idx])
        vec = np.abs(np.real(eigvecs[:, idx]))
        self.ahp_weights = vec / np.sum(vec)
        if self.n > 2:
            ci = (lambda_max - self.n) / (self.n - 1)
            self.ahp_cr = max(0.0, ci / self.RI_dict.get(self.n, 1.49))
        else:
            self.ahp_cr = 0.0
        return self.ahp_weights

    def extract_bwm(self):
        if self.ahp_weights is None:
            self.calculate_ahp()
        max_w, min_w = np.max(self.ahp_weights), np.min(self.ahp_weights)
        best_idx = np.where(np.abs(self.ahp_weights - max_w) < 1e-5)[0]
        worst_idx = np.where(np.abs(self.ahp_weights - min_w) < 1e-5)[0]
        self.best_criteria = [self.names[i] for i in best_idx]
        self.worst_criteria = [self.names[i] for i in worst_idx]
        b_idx, w_idx = best_idx[0], worst_idx[0]
        self.bo_vector = self.matrix[b_idx, :]
        self.ow_vector = self.matrix[:, w_idx]
        # BWM assumption: a_Bj >= 1 and a_jW >= 1
        self.assumption_ok = bool(np.all(self.bo_vector >= 1 - 1e-9) and np.all(self.ow_vector >= 1 - 1e-9))
        return self.best_criteria, self.worst_criteria, self.bo_vector, self.ow_vector

    def extract_fucom(self):
        if self.bo_vector is None:
            self.extract_bwm()
        ranked = np.argsort(self.bo_vector, kind="stable")  # ties broken by index order
        self.fucom_seq = [self.names[i] for i in ranked]
        phi, s = [], self.fucom_seq[0]
        for i in range(self.n - 1):
            ck, ck1 = ranked[i], ranked[i + 1]
            p = self.bo_vector[ck1] / self.bo_vector[ck]
            phi.append(p)
            s += (" = " if abs(p - 1.0) < 1e-5 else " > ") + self.fucom_seq[i + 1]
        self.fucom_phi = np.array(phi)
        self.fucom_ranking_str = s
        return self.fucom_seq, self.fucom_phi

    def _bwm_vectors(self):
        if self.bo_vector is None:
            self.extract_bwm()
        b_idx = self.names.index(self.best_criteria[0])
        w_idx = self.names.index(self.worst_criteria[0])
        return b_idx, w_idx, self.bo_vector, self.ow_vector

    def solve_bwm_linear(self):
        """Linear BWM (Rezaei, 2016): min xi s.t. |w_B - a_Bj w_j| <= xi, |w_j - a_jW w_W| <= xi."""
        b_idx, w_idx, bo, ow = self._bwm_vectors()
        nv = self.n + 1
        c = np.zeros(nv); c[-1] = 1.0
        A_ub, b_ub = [], []
        for j in range(self.n):
            for rs in (1.0, -1.0):
                row = np.zeros(nv)
                row[b_idx] += rs; row[j] -= rs * bo[j]; row[-1] = -1.0
                A_ub.append(row); b_ub.append(0.0)
            for rs in (1.0, -1.0):
                row = np.zeros(nv)
                row[j] += rs; row[w_idx] -= rs * ow[j]; row[-1] = -1.0
                A_ub.append(row); b_ub.append(0.0)
        A_eq = np.zeros((1, nv)); A_eq[0, :-1] = 1.0
        bounds = [(0.0, 1.0)] * self.n + [(0.0, None)]
        res = linprog(c, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                      A_eq=A_eq, b_eq=[1.0], bounds=bounds, method='highs')
        if not res.success:
            raise RuntimeError(f"{self.label}: linear BWM failed ({res.message})")
        self.bwm_weights_linear = res.x[:-1]
        self.bwm_xi_linear = float(res.x[-1])
        return self.bwm_weights_linear, self.bwm_xi_linear

    def _nl_constraints(self, xi):
        """For a fixed xi, the constraints of the non-linear BWM are linear in w:
        (a_Bj - xi) w_j <= w_B <= (a_Bj + xi) w_j  and  (a_jW - xi) w_W <= w_j <= (a_jW + xi) w_W."""
        b_idx, w_idx, bo, ow = self._bwm_vectors()
        A = []
        for j in range(self.n):
            r = np.zeros(self.n); r[b_idx] += 1; r[j] -= (bo[j] + xi); A.append(r)
            r = np.zeros(self.n); r[j] += (bo[j] - xi); r[b_idx] -= 1; A.append(r)
            r = np.zeros(self.n); r[j] += 1; r[w_idx] -= (ow[j] + xi); A.append(r)
            r = np.zeros(self.n); r[w_idx] += (ow[j] - xi); r[j] -= 1; A.append(r)
        return np.array(A), np.zeros(len(A))

    def _nl_feasible(self, xi):
        A, b = self._nl_constraints(xi)
        res = linprog(np.zeros(self.n), A_ub=A, b_ub=b, A_eq=np.ones((1, self.n)), b_eq=[1.0],
                      bounds=[(1e-9, 1.0)] * self.n, method='highs')
        return res.status == 0

    def solve_bwm_nonlinear(self, tol=1e-10):
        """Non-linear BWM (Rezaei, 2015): min xi s.t. |w_B/w_j - a_Bj| <= xi, |w_j/w_W - a_jW| <= xi.
        Solved exactly by bisection on xi (feasibility LPs). Because the optimum may not be unique,
        the interval [min, max] of each weight at xi* is computed; the point solution is the
        normalised interval centre (Rezaei, 2015)."""
        _, _, bo, ow = self._bwm_vectors()
        lo, hi = 0.0, float(max(np.max(np.abs(bo - 1)), np.max(np.abs(ow - 1)))) + 1.0
        if self._nl_feasible(0.0):
            hi = 0.0
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            if self._nl_feasible(mid):
                hi = mid
            else:
                lo = mid
        xi_star = hi
        A, b = self._nl_constraints(xi_star + 1e-7)  # small slack for LP tolerance
        lower, upper = np.empty(self.n), np.empty(self.n)
        for j in range(self.n):
            for sgn, store in ((1.0, lower), (-1.0, upper)):
                c = np.zeros(self.n); c[j] = sgn
                res = linprog(c, A_ub=A, b_ub=b, A_eq=np.ones((1, self.n)), b_eq=[1.0],
                              bounds=[(1e-9, 1.0)] * self.n, method='highs')
                if res.status != 0:
                    # retry with a slightly larger slack (LP tolerance at the optimum)
                    A2_, b2_ = self._nl_constraints(xi_star * (1 + 1e-6) + 1e-6)
                    res = linprog(c, A_ub=A2_, b_ub=b2_, A_eq=np.ones((1, self.n)), b_eq=[1.0],
                                  bounds=[(1e-9, 1.0)] * self.n, method='highs')
                if res.status != 0:
                    raise RuntimeError(f"{self.label}: non-linear BWM interval LP failed ({res.message})")
                store[j] = res.x[j]
        centre = 0.5 * (lower + upper)
        self.bwm_weights_nonlinear = centre / centre.sum()
        self.bwm_xi_nonlinear = float(xi_star)
        self.bwm_nl_interval = np.vstack([lower, upper]).T
        self.bwm_nl_max_width = float(np.max(upper - lower))
        return self.bwm_weights_nonlinear, self.bwm_xi_nonlinear

    def solve_bwm(self):
        """Solve both BWM variants; the one matching the source (self.bwm_model) is used for validation."""
        self.solve_bwm_linear()
        self.solve_bwm_nonlinear()
        if self.bwm_model == "linear":
            self.bwm_weights, self.bwm_xi = self.bwm_weights_linear, self.bwm_xi_linear
        else:
            self.bwm_weights, self.bwm_xi = self.bwm_weights_nonlinear, self.bwm_xi_nonlinear
        return self.bwm_weights, self.bwm_xi

    def solve_fucom(self):
        if self.fucom_phi is None:
            self.extract_fucom()
        seq_idx = [self.names.index(nm) for nm in self.fucom_seq]

        def obj(x):
            return x[-1]

        cons = [{'type': 'eq', 'fun': lambda x: np.sum(x[:-1]) - 1}]
        for i in range(self.n - 1):
            ck, ck1, p = seq_idx[i], seq_idx[i + 1], self.fucom_phi[i]
            cons.append({'type': 'ineq', 'fun': lambda x, ck=ck, ck1=ck1, p=p: x[-1] - (x[ck] / x[ck1] - p)})
            cons.append({'type': 'ineq', 'fun': lambda x, ck=ck, ck1=ck1, p=p: x[-1] - (p - x[ck] / x[ck1])})
        for i in range(self.n - 2):
            ck, ck2 = seq_idx[i], seq_idx[i + 2]
            pm = self.fucom_phi[i] * self.fucom_phi[i + 1]
            cons.append({'type': 'ineq', 'fun': lambda x, ck=ck, ck2=ck2, p=pm: x[-1] - (x[ck] / x[ck2] - p)})
            cons.append({'type': 'ineq', 'fun': lambda x, ck=ck, ck2=ck2, p=pm: x[-1] - (p - x[ck] / x[ck2])})
        bounds = [(1e-5, 1)] * self.n + [(0, None)]

        # Starting points: closed-form solution (Eq. 8) and the uniform vector.
        u = np.ones(self.n)
        for i in range(self.n - 1):
            u[i + 1] = u[i] / self.fucom_phi[i]
        w0 = u / np.sum(u)
        x0_seq = np.empty(self.n)
        for i, idx in enumerate(seq_idx):
            x0_seq[idx] = w0[i]
        x0 = np.append(x0_seq, 0.05)

        best = None
        runs = {}
        for tag, xs in (("closed_form", x0), ("uniform", np.append(np.ones(self.n) / self.n, 0.1))):
            r = minimize(obj, xs, method='SLSQP', bounds=bounds, constraints=cons,
                         options={'maxiter': 1000, 'ftol': 1e-14})
            runs[tag] = r
            if best is None or r.fun < best.fun:
                best = r
        self.fucom_weights = best.x[:-1]
        self.fucom_chi = float(max(best.x[-1], 0.0))
        # The uniform-start solution is an independent numerical check of Eq. (8).
        self.fucom_chi_uniform = float(max(runs["uniform"].x[-1], 0.0))
        self.fucom_closed_form_gap = float(np.max(np.abs(runs["uniform"].x[:-1] - w0_by_index(self, w0, seq_idx))))
        return self.fucom_weights, self.fucom_chi

    def check_bw_match(self):
        """Compare the derived best/worst criteria with those reported in the source."""
        auto_best, auto_worst = self.best_criteria[0], self.worst_criteria[0]
        self.bw_matches = []
        for entry in self.reported_bw:
            best_ok = entry["best"] == auto_best
            worst_ok = entry["worst"] == auto_worst
            self.bw_matches.append({
                "source": entry["source"],
                "reported_best": entry["best"], "reported_worst": entry["worst"],
                "auto_best": auto_best, "auto_worst": auto_worst,
                "best_ok": best_ok, "worst_ok": worst_ok, "full_match": best_ok and worst_ok,
            })
        return self.bw_matches

    def solve_all(self):
        self.calculate_ahp()
        self.extract_bwm()
        self.solve_bwm()
        self.extract_fucom()
        self.solve_fucom()
        self.check_bw_match()
        return self


def w0_by_index(ex, w0, seq_idx):
    """Map weights given in FUCOM-sequence order back to criterion-index order."""
    out = np.empty(ex.n)
    for i, idx in enumerate(seq_idx):
        out[idx] = w0[i]
    return out


def _ranks(weights, decimals=4):
    """Ranks of weights rounded to `decimals`; ties receive average ranks (1 = most important)."""
    return rankdata(-np.round(np.asarray(weights, dtype=float), decimals), method="average")


def rank_metrics(calc_weights, rep_weights, decimals=4):
    """Rank-based agreement measures (Section 3.6).
    rho: Pearson correlation between rank vectors (equal to Eq. 10 without ties).
    WS:  Eq. (11), with the reported ranking as the reference."""
    Rd, Rr = _ranks(calc_weights, decimals), _ranks(rep_weights, decimals)
    n = len(Rr)
    rho = 1.0 if np.allclose(Rd, Rr) else float(np.corrcoef(Rd, Rr)[0, 1])
    ws = 1.0 - sum(2.0 ** (-Rr[j]) * abs(Rr[j] - Rd[j]) / max(abs(1 - Rr[j]), abs(n - Rr[j]))
                   for j in range(n))
    return rho, float(ws)


def all_metrics(calc_weights, rep_weights):
    """Pearson r, Spearman rho, WS, MAE and maximum absolute deviation (Eqs. 9-12)."""
    c, r = np.asarray(calc_weights, float), np.asarray(rep_weights, float)
    rho, ws = rank_metrics(c, r)
    return {"r": float(pearsonr(c, r)[0]), "rho": rho, "ws": ws,
            "mae": float(np.mean(np.abs(c - r))), "maxdev": float(np.max(np.abs(c - r)))}


# ============================================================================ #
# 1) DATA: matrices and reported weights transcribed from the source tables
# ============================================================================ #
# --- Case 1: Pamucar et al. (2018), Example 3 -- Tables 4 and 6 ------------- #
# Exception to the use of published matrices: Table 4 of the source prints
# a23 = 4 and a32 = 1/4, which contradicts the source itself. With these entries
# the matrix gives CR = 0.105 and AHP weights that differ from the reported ones
# (r = 0.918), whereas the source reports CR = 0.000, identical weights for C3,
# C5 and C7, and the Best-to-Others vector (2, 4, 1, 4, 1, 4, 1, 4) read from the
# row of C3 (or the identical rows of C5 and C7). The two entries are therefore a
# transposition error and are entered as a23 = 1/4 and a32 = 4, which makes rows
# C3, C5 and C7 identical, as stated in the source.
A1 = [
    [1,   2,   1/2, 2,   1/2, 2,   1/2, 2],
    [1/2, 1,   1/4, 1,   1/4, 1,   1/4, 1],
    [2,   4,   1,   4,   1,   4,   1,   4],
    [1/2, 1,   1/4, 1,   1/4, 1,   1/4, 1],
    [2,   4,   1,   4,   1,   4,   1,   4],
    [1/2, 1,   1/4, 1,   1/4, 1,   1/4, 1],
    [2,   4,   1,   4,   1,   4,   1,   4],
    [1/2, 1,   1/4, 1,   1/4, 1,   1/4, 1]
]
ex1 = MCDM_Extractor(
    matrix=A1, criteria_names=["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"],
    reported_weights={
        "AHP":   [0.1111, 0.0556, 0.2222, 0.0556, 0.2222, 0.0556, 0.2222, 0.0556],
        "BWM":   [0.1111, 0.0556, 0.2222, 0.0556, 0.2222, 0.0556, 0.2222, 0.0556],
        "FUCOM": [0.1111, 0.0556, 0.2222, 0.0556, 0.2222, 0.0556, 0.2222, 0.0556],
    },
    label="Case 1 - Resource allocation in cloud computing (Pamucar et al. 2018, 8 criteria, Tables 4/6)"
).solve_all()

# --- Case 2: Pamucar et al. (2018), Example 4 -- Tables 7 and 9 ------------- #
A2 = [
    [1,    1/3,  3,    1,    5],
    [3,    1,    5,    3,    7],
    [1/3,  1/5,  1,    1/3,  3],
    [1,    1/3,  3,    1,    5],
    [1/5,  1/7,  1/3,  1/5,  1],
]
ex2 = MCDM_Extractor(
    matrix=A2, criteria_names=["C1", "C2", "C3", "C4", "C5"],
    reported_weights={
        "AHP":   [0.2017, 0.4641, 0.0888, 0.2017, 0.0436],
        "BWM":   [0.2105, 0.4211, 0.1053, 0.2105, 0.0526],
        "FUCOM": [0.2105, 0.4211, 0.1053, 0.2105, 0.0526],
    },
    # Reported BO = (2,1,4,2,8) and OW = (4,8,2,4,1), taken from Rezaei (2016)
    reported_bw=[
        {"source": "BWM vectors (Rezaei 2016, separate source)", "best": "C2", "worst": "C5"},
    ],
    label="Case 2 - Car purchase (Pamucar et al. 2018, 5 criteria, Tables 7/9)"
).solve_all()

# --- Case 3: Ekin & Sarul (2022) -- Smart-city components ------------------- #
# Table 8: AHP matrix (geometric mean of three experts), entered as published.
# Both triangles are rounded to two decimals in the source, so lower-triangle
# entries deviate slightly from the reciprocals of the upper-triangle entries.
# Criteria: C1 = smart economy, C2 = smart people, C3 = smart governance,
#           C4 = smart mobility, C5 = smart environment, C6 = smart living
A3 = [
    [1.00, 1.59, 0.44, 1.06, 0.58, 0.75],
    [0.63, 1.00, 0.44, 0.64, 0.38, 0.35],
    [2.29, 2.29, 1.00, 1.59, 1.31, 1.49],
    [0.94, 1.55, 0.63, 1.00, 0.71, 1.06],
    [1.71, 2.62, 0.76, 1.41, 1.00, 1.44],
    [1.34, 2.88, 0.67, 0.94, 0.69, 1.00],
]
ex3 = MCDM_Extractor(
    matrix=A3, criteria_names=["C1", "C2", "C3", "C4", "C5", "C6"],
    reported_weights={
        # Table 10 -- reported AHP priorities
        "AHP":   [0.13, 0.08, 0.25, 0.15, 0.22, 0.17],
        # Table 13 -- arithmetic mean of the three experts' BWM weights
        "BWM":   [0.1077, 0.0693, 0.2677, 0.1903, 0.2193, 0.1460],
        # Table 15 -- arithmetic mean of the three experts' FUCOM weights
        "FUCOM": [0.0983, 0.0727, 0.2710, 0.2010, 0.2230, 0.1340],
    },
    # Tables 13 and 15 -- best/worst criteria selected by each expert (same in both methods)
    reported_bw=[
        {"source": "Expert 1 (BWM/FUCOM)", "best": "C4", "worst": "C2"},
        {"source": "Expert 2 (BWM/FUCOM)", "best": "C5", "worst": "C4"},
        {"source": "Expert 3 (BWM/FUCOM)", "best": "C3", "worst": "C2"},
    ],
    label="Case 3 - Smart-city components (Ekin & Sarul 2022, 6 criteria, Tables 8/10/13/15)"
).solve_all()

# --- Case 4: Mercan & Acibuca (2025) -- Beekeeping land suitability -------- #
# PLOS ONE 20(11), e0335784. Table 4: AHP matrix (9 criteria, consensus panel).
# Criteria: C1 = land use/land cover, C2 = distance to water sources,
#           C3 = annual precipitation, C4 = slope, C5 = aspect, C6 = elevation,
#           C7 = temperature, C8 = distance to settlements, C9 = distance to roads
A4 = [
    [1,   4,   5,   6,   7,   7,   8,   9,   9],
    [1/4, 1,   2,   2,   3,   4,   6,   6,   8],
    [1/5, 1/2, 1,   2,   3,   4,   4,   5,   7],
    [1/6, 1/2, 1/2, 1,   2,   3,   4,   5,   6],
    [1/7, 1/3, 1/3, 1/2, 1,   2,   3,   4,   5],
    [1/7, 1/4, 1/4, 1/3, 1/2, 1,   2,   3,   5],
    [1/8, 1/6, 1/4, 1/4, 1/3, 1/2, 1,   2,   3],
    [1/9, 1/6, 1/5, 1/5, 1/4, 1/3, 1/2, 1,   2],
    [1/9, 1/8, 1/7, 1/6, 1/5, 1/5, 1/3, 1/2, 1],
]
ex4 = MCDM_Extractor(
    matrix=A4, criteria_names=["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"],
    reported_weights={
        # Table 5 -- reported AHP and FUCOM weights (BWM not applied in the source)
        "AHP":   [0.401, 0.171, 0.134, 0.101, 0.069, 0.050, 0.033, 0.024, 0.017],
        "FUCOM": [0.410, 0.137, 0.102, 0.082, 0.068, 0.059, 0.051, 0.046, 0.046],
    },
    # Table 3 -- FUCOM ranking of the consensus panel
    reported_bw=[
        {"source": "Consensus panel (FUCOM order)", "best": "C1", "worst": "C9"},
    ],
    label="Case 4 - Beekeeping land suitability, Mardin (Mercan & Acibuca 2025, 9 criteria, Tables 4/5)"
).solve_all()

# --- Case 5: Badi & Abdulshahed (2019) -- Ranking of Libyan airlines ------- #
# Operational Research in Engineering Sciences: Theory and Applications 2(1), 1-14.
# Table 1: AHP matrix (5 criteria, decision-maker group).
# Criteria: C1 = reliability, C2 = employees, C3 = management,
#           C4 = customer satisfaction, C5 = tangibles
A5 = [
    [1,   5,   4,   3,   7],
    [1/5, 1,   1/2, 1/3, 1],
    [1/4, 2,   1,   1/2, 2],
    [1/3, 3,   2,   1,   3],
    [1/7, 1,   1/2, 1/3, 1],
]
ex5 = MCDM_Extractor(
    matrix=A5, criteria_names=["C1", "C2", "C3", "C4", "C5"],
    reported_weights={
        # Table 3 (differs from the weights in Table 1 by rounding)
        "AHP":   [0.506, 0.076, 0.131, 0.216, 0.071],
        "FUCOM": [0.520, 0.094, 0.104, 0.192, 0.090],
    },
    # Table 2 -- FUCOM ranking (Step 1): C1 > C4 > C3 > C2 > C5
    reported_bw=[
        {"source": "FUCOM order (Step 1)", "best": "C1", "worst": "C5"},
    ],
    label="Case 5 - Ranking of Libyan airlines (Badi & Abdulshahed 2019, 5 criteria, Tables 1/2/3)"
).solve_all()

# --- Case 6: Badi, Pamucar, Stevic & Muhammad (2023) -- Wind-farm site ----- #
# selection, Libya. Scientific African 19, e01511.
# Table 4: AHP matrices of four experts (6 criteria). The analysis uses their
# element-wise geometric mean (aggregation of individual judgements).
# Criteria: C1 = safety/quality, C2 = economy/benefit, C3 = social impression,
#           C4 = environment/ecology, C5 = regulation, C6 = policy
# Entries printed as 0.333 in the source are entered as 1/3.
A6_e1 = [
    [1,     4,     2,     1,     3,     3],
    [0.25,  1,     1,     1,     2,     2],
    [0.5,   1,     1,     2,     3,     2],
    [1,     1,     0.5,   1,     3,     2],
    [1/3,   0.5,   1/3,   1/3,   1,     2],
    [1/3,   0.5,   0.5,   0.5,   0.5,   1],
]
A6_e2 = [
    [1,     3,     1,     1,     3,     3],
    [1/3,   1,     2,     1,     3,     2],
    [1,     0.5,   1,     2,     2,     3],
    [1,     1,     0.5,   1,     3,     2],
    [1/3,   1/3,   0.5,   1/3,   1,     2],
    [1/3,   0.5,   1/3,   0.5,   0.5,   1],
]
A6_e3 = [
    [1,     4,     2,     2,     3,     3],
    [0.25,  1,     1,     2,     3,     2],
    [0.5,   1,     1,     2,     2,     3],
    [0.5,   0.5,   0.5,   1,     1,     3],
    [1/3,   1/3,   0.5,   1,     1,     1],
    [1/3,   0.5,   1/3,   1/3,   1,     1],
]
A6_e4 = [
    [1,     2,     2,     3,     3,     3],
    [0.5,   1,     1,     1,     2,     2],
    [0.5,   1,     1,     2,     2,     3],
    [1/3,   1,     0.5,   1,     1,     3],
    [1/3,   0.5,   0.5,   1,     1,     1],
    [1/3,   0.5,   1/3,   1/3,   1,     1],
]
A6 = (np.array(A6_e1) * np.array(A6_e2) * np.array(A6_e3) * np.array(A6_e4)) ** 0.25
ex6 = MCDM_Extractor(
    matrix=A6, criteria_names=["C1", "C2", "C3", "C4", "C5", "C6"],
    reported_weights={
        # Table 6 -- reported aggregate AHP weights
        "AHP": [0.307, 0.172, 0.199, 0.157, 0.091, 0.074],
        # Table 3 -- BWM weights obtained from the averaged BO/OW vectors (Table 2)
        "BWM": [0.422, 0.177, 0.230, 0.081, 0.048, 0.043],
    },
    # Table 2 -- common selection of the four experts
    reported_bw=[
        {"source": "BWM survey (common choice of 4 experts, Table 2)", "best": "C1", "worst": "C6"},
    ],
    # The source solved the non-linear (ratio-based) BWM model (model printed after Table 4).
    bwm_model="nonlinear",
    label="Case 6 - Wind-farm site selection, Libya (Badi et al. 2023, 6 criteria, Tables 2/3/4/6)"
).solve_all()

# Individual expert matrices (Table 4): CR recomputed and compared with Table 5 of
# the source (CR = 0.084, 0.083, 0.046, 0.036). The published values could not be
# reproduced from the printed matrices (eigenvector, normalised-column and
# geometric-mean approximations give the same result); all recomputed values are
# below 0.10, so the source's conclusion that every matrix is acceptable holds.
ex6.expert_cr_reported = [0.084, 0.083, 0.046, 0.036]
_tmp = []
for _m in (A6_e1, A6_e2, A6_e3, A6_e4):
    _e = MCDM_Extractor(_m)
    _e.calculate_ahp()
    _tmp.append(float(_e.ahp_cr))
ex6.expert_cr = _tmp

# --- Case 7: Vujovic et al. (2025) -- Flash-flood susceptibility ----------- #
# Ibar River Basin, Montenegro. Water 17(6), 844.
# Table 6: AHP matrix (6 criteria).
# Criteria: C1 = slope, C2 = distance to drainage network, C3 = geology,
#           C4 = land cover, C5 = drainage density, C6 = bare soil index
A7 = [
    [1,    2,    3,    3,    3,    3],
    [0.5,  1,    2,    3,    3,    2],
    [1/3,  0.5,  1,    2,    3,    2],
    [1/3,  1/3,  0.5,  1,    3,    2],
    [1/3,  1/3,  1/3,  1/3,  1,    2],
    [1/3,  0.5,  0.5,  0.5,  0.5,  1],
]
ex7 = MCDM_Extractor(
    matrix=A7, criteria_names=["C1", "C2", "C3", "C4", "C5", "C6"],
    reported_weights={
        # Table 6 -- weight (%) / 100
        "AHP": [0.328, 0.226, 0.159, 0.125, 0.085, 0.077],
        # Table 7 -- weight (%) / 100 (B = C1, W = C6, input-based CR = 0.20)
        "BWM": [0.3794, 0.2213, 0.1475, 0.1106, 0.0885, 0.0527],
    },
    reported_bw=[
        {"source": "Joint AHP/BWM assessment (Tables 6/7)", "best": "C1", "worst": "C6"},
    ],
    label="Case 7 - Flash-flood susceptibility, Ibar basin, Montenegro (Vujovic et al. 2025, 6 criteria, Tables 6/7)"
).solve_all()

EXAMPLES = [ex1, ex2, ex3, ex4, ex5, ex6, ex7]


def available_methods(ex):
    """Return the methods (AHP/BWM/FUCOM) for which the source reports weights."""
    return [m for m in ["AHP", "BWM", "FUCOM"] if m in ex.reported_weights]


# ============================================================================ #
# 2) CONSOLE REPORT
# ============================================================================ #
pd.set_option("display.float_format", lambda v: f"{v:0.4f}")
pd.set_option("display.width", 120)


def build_comparison_df(ex):
    rows = []
    for method in available_methods(ex):
        calc = {"AHP": ex.ahp_weights, "BWM": ex.bwm_weights, "FUCOM": ex.fucom_weights}[method]
        rep = np.array(ex.reported_weights[method])
        for i, name in enumerate(ex.names):
            rows.append({
                "Method": method, "Criterion": name,
                "Derived": calc[i], "Reported": rep[i],
                "AbsDiff": abs(calc[i] - rep[i]),
            })
    return pd.DataFrame(rows)


def print_report(ex):
    print("\n" + "=" * 78)
    print(f" {ex.label}")
    print("=" * 78)
    print(f" AHP consistency ratio (CR): {ex.ahp_cr:0.4f}  "
          f"{'(acceptable, CR < 0.10)' if ex.ahp_cr < 0.10 else '(CR >= 0.10)'}")
    print(f" Best / worst criterion    : {ex.best_criteria[0]} / {ex.worst_criteria[0]}")
    print(f" FUCOM priority order      : {ex.fucom_ranking_str}")
    print(f" BWM model used for validation: {ex.bwm_model}   xi^L* (linear) = {ex.bwm_xi_linear:0.4f}   "
          f"xi* (non-linear) = {ex.bwm_xi_nonlinear:0.4f}   max NL interval width = {ex.bwm_nl_max_width:.1e}")
    print(f" xi (used) = {ex.bwm_xi:0.4f}   chi* (FUCOM) = {ex.fucom_chi:0.4f}   "
          f"chi* (uniform start) = {ex.fucom_chi_uniform:0.2e}")
    if not ex.assumption_ok:
        print(" WARNING: a_Bj < 1 or a_jW < 1 for at least one criterion; the BWM assumption is violated.")
    if ex.nonreciprocal_cells:
        print(f" Entries deviating from exact reciprocity (used as published): {len(ex.nonreciprocal_cells)}"
              f"  (max |a_ij * a_ji - 1| = {ex.max_reciprocity_dev:.4f})")

    methods_here = available_methods(ex)
    df = build_comparison_df(ex)
    pivot = df.pivot(index="Criterion", columns="Method", values=["Derived", "Reported", "AbsDiff"])
    pivot = pivot.reindex(ex.names)
    pivot = pivot.reorder_levels([1, 0], axis=1).sort_index(axis=1, level=0)
    pivot = pivot[methods_here]
    with pd.option_context("display.max_columns", None):
        print("\n" + pivot.to_string())

    if len(methods_here) < 3:
        print(f"\n (Only {' and '.join(methods_here)} reported in this source; "
              f"other methods are not compared.)")

    print("\n Agreement (derived vs. reported):")
    for method in methods_here:
        calc = {"AHP": ex.ahp_weights, "BWM": ex.bwm_weights, "FUCOM": ex.fucom_weights}[method]
        rep = np.array(ex.reported_weights[method])
        mt = all_metrics(calc, rep)
        r_pearson, rho_spearman, ws = mt["r"], mt["rho"], mt["ws"]
        agreement = ("EXCELLENT" if r_pearson >= 0.99 else
                 "VERY HIGH" if r_pearson >= 0.95 else
                 "HIGH" if r_pearson >= 0.80 else
                 "MODERATE" if r_pearson >= 0.50 else "WEAK")
        print(f"   {method:6s}: Pearson r = {r_pearson:0.4f}   Spearman rho = {rho_spearman:0.4f}   "
              f"WS = {ws:0.4f}   MAE = {mt['mae']:0.4f}   max|dw| = {mt['maxdev']:0.4f}   [{agreement}]")

    if ex.reported_bw:
        print("\n Best/worst match (derived "
              f"'{ex.best_criteria[0]}'/'{ex.worst_criteria[0]}' vs. selections reported in the source):")
        for m in ex.bw_matches:
            match_label = "FULL MATCH" if m["full_match"] else (
                "PARTIAL (worst)" if m["worst_ok"] else
                "PARTIAL (best)" if m["best_ok"] else "MISMATCH")
            print(f"   {m['source']:24s}: reported best={m['reported_best']}, worst={m['reported_worst']}"
                  f"   ->  [{match_label}]")
        n_full = sum(m["full_match"] for m in ex.bw_matches)
        print(f"   Full matches: {n_full}/{len(ex.bw_matches)}")
    else:
        print("\n Best/worst match: no best/worst selection reported in the source; check skipped.")

    if getattr(ex, "expert_cr", None):
        print("\n Individual expert matrices (CR recomputed vs. reported in the source):")
        for k, (cr, cr_rep) in enumerate(zip(ex.expert_cr, ex.expert_cr_reported), start=1):
            print(f"   Expert {k}: CR = {cr:0.3f}   reported = {cr_rep:0.3f}   difference = {cr - cr_rep:+0.3f}")
        print("   The reported values are not reproduced by the printed matrices; all CR < 0.10.")


for _ex in EXAMPLES:
    print_report(_ex)


# ============================================================================ #
# 3) FIGURES
# ============================================================================ #
COLOR_AHP = "#2a78d6"
COLOR_BWM = "#eb6834"
COLOR_FUCOM = "#1baf7a"
COLOR_REPORTED = "#0b0b0b"
COLOR_GRID = "#d8d8d5"
COLOR_TEXT = "#3a3a38"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.edgecolor": "#9a9a96",
    "axes.labelcolor": COLOR_TEXT,
    "text.color": COLOR_TEXT,
    "xtick.color": COLOR_TEXT,
    "ytick.color": COLOR_TEXT,
})


METHOD_COLORS = {"AHP": COLOR_AHP, "BWM": COLOR_BWM, "FUCOM": COLOR_FUCOM}


def plot_professional(ex, save_path):
    methods = available_methods(ex)
    m = len(methods)
    calc_map = {"AHP": ex.ahp_weights, "BWM": ex.bwm_weights, "FUCOM": ex.fucom_weights}
    x = np.arange(ex.n)
    width = 0.8 / m
    offsets = [-0.4 + width * (i + 0.5) for i in range(m)]
    fig, ax = plt.subplots(figsize=(11, 6), dpi=160)

    for method, off in zip(methods, offsets):
        ax.bar(x + off, calc_map[method], width, label=method, color=METHOD_COLORS[method],
               zorder=3, linewidth=0)

    # Reported weights are shown as horizontal markers.
    for i, (method, off) in enumerate(zip(methods, offsets)):
        rep = np.array(ex.reported_weights[method])
        ax.scatter(x + off, rep, marker="_", s=420, linewidths=2.4,
                   color=COLOR_REPORTED, zorder=4,
                   label="Reported (source study)" if i == 0 else None)

    ax.set_ylabel("Criterion weight", fontsize=11, fontweight="bold")
    ax.set_title(ex.label, fontsize=13, fontweight="bold", pad=14)
    ax.set_xticks(x)
    ax.set_xticklabels(ex.names, fontsize=10.5, fontweight="bold")
    ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.yaxis.grid(True, color=COLOR_GRID, linewidth=0.9, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", length=0)

    ax.legend(frameon=False, fontsize=10, loc="upper right", ncols=m + 1, bbox_to_anchor=(1.0, 1.12))

    footnote_parts = []
    for method in methods:
        r_val, _ = pearsonr(calc_map[method], np.array(ex.reported_weights[method]))
        ws_val = all_metrics(calc_map[method], ex.reported_weights[method])["ws"]
        footnote_parts.append(f"{method}: r={r_val:.3f} WS={ws_val:.3f}")
    footnote = "Pearson r / WS coefficient (vs. reported weights)  ·  " + "   ".join(footnote_parts)
    fig.text(0.5, -0.02, footnote, ha="center", fontsize=9, color="#6b6a66")

    fig.tight_layout()
    fig.savefig(save_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


plot_professional(ex1, "figure_case1.png")
plot_professional(ex2, "figure_case2.png")
plot_professional(ex3, "figure_case3.png")
plot_professional(ex4, "figure_case4.png")
plot_professional(ex5, "figure_case5.png")
plot_professional(ex6, "figure_case6.png")
plot_professional(ex7, "figure_case7.png")
print("\nFigures saved: figure_case1.png ... figure_case7.png")


# ============================================================================ #
# 4) EXCEL REPORT
# ============================================================================ #
FONT_NAME = "Arial"
NAVY = "1F3864"
LIGHT_BLUE = "D9E2F3"
HEADER_FILL = PatternFill("solid", fgColor=NAVY)
HEADER_FONT = Font(name=FONT_NAME, size=11, bold=True, color="FFFFFF")
TITLE_FONT = Font(name=FONT_NAME, size=14, bold=True, color=NAVY)
SUB_FONT = Font(name=FONT_NAME, size=10, italic=True, color="595959")
BODY_FONT = Font(name=FONT_NAME, size=10.5)
MAE_FONT = Font(name=FONT_NAME, size=10.5, bold=True)
THIN = Side(style="thin", color="BFBFBF")
BORDER = Border(left=THIN, right=THIN, top=THIN, bottom=THIN)
CENTER = Alignment(horizontal="center", vertical="center")


def set_landscape_fit(ws):
    ws.page_setup.orientation = "landscape"
    ws.page_setup.fitToWidth = 1
    ws.page_setup.fitToHeight = 0
    ws.sheet_properties.pageSetUpPr.fitToPage = True
    ws.page_margins = PageMargins(left=0.4, right=0.4, top=0.5, bottom=0.5)


wb = Workbook()

# ---- Sheet: Summary -------------------------------------------------------- #
ws0 = wb.active
ws0.title = "Summary"
ws0.sheet_view.showGridLines = False
set_landscape_fit(ws0)
ws0.column_dimensions["A"].width = 4
ws0.column_dimensions["B"].width = 46
ws0.column_dimensions["C"].width = 14
ws0.column_dimensions["D"].width = 20
ws0.column_dimensions["E"].width = 14
ws0.column_dimensions["F"].width = 62
ws0.column_dimensions["G"].width = 12

ws0["B2"] = "Derived vs. Reported Weights - Summary"
ws0["B2"].font = Font(name=FONT_NAME, size=16, bold=True, color=NAVY)
ws0["B3"] = ("Pamucar, Stevic & Sremac (2018)  +  Ekin & Sarul (2022)  +  Mercan & Acibuca (2025)  +  "
             "Badi & Abdulshahed (2019)  +  Badi et al. (2023)  +  Vujovic et al. (2025) -- "
             "derivation of BWM and FUCOM weights from AHP matrices")
ws0["B3"].font = SUB_FONT
ws0.merge_cells("B2:F2")
ws0.merge_cells("B3:F3")

headers = ["Case", "Method", "Pearson r (vs. reported)", "Agreement", "Source", "WS"]
r0 = 5
for j, h in enumerate(headers):
    col_idx = 2 + j if j < 5 else 7  # column F: source note, column G: WS
    c = ws0.cell(row=r0, column=col_idx, value=h)
    c.font = HEADER_FONT
    c.fill = HEADER_FILL
    c.alignment = CENTER
    c.border = BORDER

row = r0 + 1
summary_rows = []
SHEET_NAME_MAP = {id(ex1): "Case 1", id(ex2): "Case 2", id(ex3): "Case 3", id(ex4): "Case 4",
                   id(ex5): "Case 5", id(ex6): "Case 6", id(ex7): "Case 7"}

# Source of each case (column "Source")
CITATION_MAP = {
    id(ex1): "Pamucar, Stevic & Sremac (2018), Symmetry 10(9), 393 -- Example 3 (Tables 4/6)",
    id(ex2): "Pamucar, Stevic & Sremac (2018), Symmetry 10(9), 393 -- Example 4 (Tables 7/9)",
    id(ex3): "Ekin & Sarul (2022), Alphanumeric Journal 10(2) -- \"Investigation of Smart City "
             "Components by AHP-BWM-FUCOM and DEMATEL Methods\" (Tables 8/10/13/15)",
    id(ex4): "Mercan & Acibuca (2025), PLOS ONE 20(11) -- \"Utilizing FUCOM and AHP methods to "
             "identify the optimal beekeeping lands: A case study from Mardin, Turkiye\" (Tables 4/5)",
    id(ex5): "Badi & Abdulshahed (2019), Operational Research in Engineering Sciences: Theory and "
             "Applications 2(1), 1-14 -- \"Ranking the Libyan airlines by using FUCOM and AHP\" (Tables 1/2/3)",
    id(ex6): "Badi, Pamucar, Stevic & Muhammad (2023), Scientific African 19, e01511 -- "
             "\"Wind farm site selection using BWM-AHP-MARCOS method: A case study of Libya\" (Tables 2/3/4/6)",
    id(ex7): "Vujovic et al. (2025), Water 17(6), 844 -- \"A Comparison of the AHP and BWM Models "
             "for the Flash Flood Susceptibility Assessment: A Case Study of the Ibar River Basin "
             "in Montenegro\" (Tables 6/7)",
}

for ex in EXAMPLES:
    sheet_name = SHEET_NAME_MAP[id(ex)]
    for method in available_methods(ex):
        summary_rows.append((ex, sheet_name, method))

for ex, sheet_name, method in summary_rows:
    ws0.cell(row=row, column=2, value=ex.label).font = BODY_FONT
    ws0.cell(row=row, column=3, value=method).font = BODY_FONT
    ws0.cell(row=row, column=3).alignment = CENTER
    ws0.cell(row=row, column=4).number_format = "0.00000"
    ws0.cell(row=row, column=4).alignment = CENTER
    ws0.cell(row=row, column=5).alignment = CENTER
    ws0.cell(row=row, column=6, value=CITATION_MAP[id(ex)]).font = Font(name=FONT_NAME, size=9, color="595959")
    ws0.cell(row=row, column=6).alignment = Alignment(horizontal="left", vertical="center", wrap_text=True)
    ws0.row_dimensions[row].height = 30
    for col in range(2, 7):
        ws0.cell(row=row, column=col).border = BORDER
    row += 1

ws0.freeze_panes = "B6"
NOTE_ROW = row + 1
ws0.cell(row=NOTE_ROW, column=2,
         value=("Note: All matrices are used as published except Case 1, whose printed entries a23 = 4 and a32 = 1/4 are "
                "transposed (with them CR = 0.105 and the reported AHP weights are not reproduced); after correction the matrix "
                "is fully consistent (CR = 0.000), and the derived BWM/FUCOM weights coincide with the reported ones. In Case 2 (CR = 0.028), "
                "the reported BWM/FUCOM weights derive from the BO/OW vectors (2,1,4,2,8)/(4,8,2,4,1) of Rezaei (2016), not from "
                "the AHP matrix. In Case 3, the reported BWM and FUCOM weights are arithmetic means, computed in this study, of "
                "three experts' separately elicited weights; the source matrix is used as published (both triangles rounded "
                "to two decimals). In Cases 4 and 5, the source reports its FUCOM priorities directly "
                "((1,3,4,5,6,7,8,9,9) and (1; 2.7; 5; 5.5; 5.8)), whereas the AHP best rows are (1,4,5,6,7,7,8,9,9) and (1,5,4,3,7). "
                "The Case 5 matrix is reciprocal in the source (a25 = a52 = 1). In Case 6, the AHP matrix is the element-wise "
                "geometric mean of the four experts' matrices (Table 4), each transcribed cell by cell; the individual CRs recomputed from these matrices (0.056, 0.067, 0.050, 0.024) differ from those reported in Table 5 of the source (0.084, 0.083, 0.046, 0.036) but are all below 0.10. The source obtained its BWM weights from the averaged "
                "BO/OW vectors with the non-linear BWM. In Case 7, the source elicited its BWM inputs separately "
                "(BO = 1,2,3,4,5,6); the reported CR = 0.20 is the input-based BWM consistency ratio. Spearman rho and WS were "
                "computed from weights rounded to four decimals, with tied weights assigned average ranks."))
ws0.cell(row=NOTE_ROW, column=2).font = Font(name=FONT_NAME, size=9.5, italic=True, color="595959")
ws0.cell(row=NOTE_ROW, column=2).alignment = Alignment(wrap_text=True, vertical="top")
ws0.merge_cells(start_row=NOTE_ROW, start_column=2, end_row=NOTE_ROW + 4, end_column=6)


# ---- Sheets: one per case ------------------------------------------------- #
def build_example_sheet(wb, ex, sheet_name):
    ws = wb.create_sheet(sheet_name)
    ws.sheet_view.showGridLines = False
    set_landscape_fit(ws)
    ws.column_dimensions["A"].width = 10
    for col in "BCDEFGHIJK":
        ws.column_dimensions[col].width = 14

    ws["A1"] = ex.label
    ws["A1"].font = TITLE_FONT
    ws.merge_cells("A1:J1")

    ws["A2"] = (f"AHP consistency ratio (CR) = {ex.ahp_cr:.4f}    |    xi^L* = {ex.bwm_xi_linear:.4f}    |    xi* (NL) = {ex.bwm_xi_nonlinear:.4f}    |    BWM model: {ex.bwm_model}    |    "
                f"Best criterion: {ex.best_criteria[0]}    |    Worst criterion: {ex.worst_criteria[0]}")
    ws["A2"].font = SUB_FONT
    ws.merge_cells("A2:J2")

    ws["A3"] = f"FUCOM priority order: {ex.fucom_ranking_str}"
    ws["A3"].font = SUB_FONT
    ws.merge_cells("A3:J3")

    header_row = 5
    ws.cell(row=header_row, column=1, value="Criterion").font = HEADER_FONT
    ws.cell(row=header_row, column=1).fill = HEADER_FILL
    ws.cell(row=header_row, column=1).alignment = CENTER
    ws.cell(row=header_row, column=1).border = BORDER

    methods = available_methods(ex)
    calc_map = {"AHP": ex.ahp_weights, "BWM": ex.bwm_weights, "FUCOM": ex.fucom_weights}

    col = 2
    method_col_start = {}
    for method in methods:
        method_col_start[method] = col
        ws.merge_cells(start_row=header_row, start_column=col, end_row=header_row, end_column=col + 2)
        c = ws.cell(row=header_row, column=col, value=method)
        c.font = HEADER_FONT
        c.fill = HEADER_FILL
        c.alignment = CENTER
        for k in range(3):
            ws.cell(row=header_row, column=col + k).fill = HEADER_FILL
            ws.cell(row=header_row, column=col + k).border = BORDER
        col += 3

    sub_row = header_row + 1
    ws.cell(row=sub_row, column=1, value="").border = BORDER
    for method in methods:
        base = method_col_start[method]
        for k, label in enumerate(["Derived", "Reported", "|Diff|"]):
            c = ws.cell(row=sub_row, column=base + k, value=label)
            c.font = Font(name=FONT_NAME, size=9.5, bold=True, italic=True, color=NAVY)
            c.fill = PatternFill("solid", fgColor=LIGHT_BLUE)
            c.alignment = CENTER
            c.border = BORDER

    first_data_row = sub_row + 1
    for i, name in enumerate(ex.names):
        r = first_data_row + i
        cN = ws.cell(row=r, column=1, value=name)
        cN.font = Font(name=FONT_NAME, size=10.5, bold=True)
        cN.alignment = CENTER
        cN.border = BORDER
        for method in methods:
            base = method_col_start[method]
            calc_cell = ws.cell(row=r, column=base, value=float(calc_map[method][i]))
            calc_cell.number_format = "0.0000"
            calc_cell.font = BODY_FONT
            calc_cell.alignment = CENTER
            calc_cell.border = BORDER

            rep_cell = ws.cell(row=r, column=base + 1, value=float(ex.reported_weights[method][i]))
            rep_cell.number_format = "0.0000"
            rep_cell.font = BODY_FONT
            rep_cell.alignment = CENTER
            rep_cell.border = BORDER

            calc_letter = get_column_letter(base)
            rep_letter = get_column_letter(base + 1)
            diff_cell = ws.cell(row=r, column=base + 2, value=f"=ABS({calc_letter}{r}-{rep_letter}{r})")
            diff_cell.number_format = "0.0000"
            diff_cell.font = BODY_FONT
            diff_cell.alignment = CENTER
            diff_cell.border = BORDER

    last_data_row = first_data_row + ex.n - 1
    metric_row = last_data_row + 1
    ws.cell(row=metric_row, column=1, value="Pearson r").font = MAE_FONT
    ws.cell(row=metric_row, column=1).alignment = CENTER
    ws.cell(row=metric_row, column=1).border = BORDER
    ws.cell(row=metric_row, column=1).fill = PatternFill("solid", fgColor="FCE4D6")

    corr_cell_refs = {}
    for method in methods:
        base = method_col_start[method]
        calc_letter = get_column_letter(base)
        rep_letter = get_column_letter(base + 1)
        diff_letter = get_column_letter(base + 2)
        corr_cell = ws.cell(
            row=metric_row, column=base + 2,
            value=(f"=CORREL({calc_letter}{first_data_row}:{calc_letter}{last_data_row},"
                   f"{rep_letter}{first_data_row}:{rep_letter}{last_data_row})")
        )
        corr_cell.number_format = "0.0000"
        corr_cell.font = MAE_FONT
        corr_cell.alignment = CENTER
        corr_cell.border = BORDER
        corr_cell.fill = PatternFill("solid", fgColor="FCE4D6")
        corr_cell_refs[method] = f"'{sheet_name}'!{diff_letter}{metric_row}"

        # Fixed colour scale for absolute differences: 0 / 0.01 / 0.05
        rng = f"{diff_letter}{first_data_row}:{diff_letter}{last_data_row}"
        ws.conditional_formatting.add(
            rng,
            ColorScaleRule(start_type="num", start_value=0, start_color="C6E0B4",
                            mid_type="num", mid_value=0.01, mid_color="FFEB9C",
                            end_type="num", end_value=0.05, end_color="F8696B"),
        )

    # ---- WS, Spearman rho, MAE and max|dw| are computed in Python and written as values
    ws_row = metric_row + 1
    ws.cell(row=ws_row, column=1, value="WS").font = MAE_FONT
    ws.cell(row=ws_row, column=1).alignment = CENTER
    ws.cell(row=ws_row, column=1).border = BORDER
    ws.cell(row=ws_row, column=1).fill = PatternFill("solid", fgColor="E2EFDA")

    ws_cell_refs = {}
    for method in methods:
        base = method_col_start[method]
        diff_letter = get_column_letter(base + 2)
        mt = all_metrics(calc_map[method], ex.reported_weights[method])
        ws_val = mt["ws"]
        ws_cell = ws.cell(row=ws_row, column=base + 2, value=float(ws_val))
        for k, (lab, key) in enumerate([("Spearman rho", "rho"), ("MAE", "mae"), ("max|dw|", "maxdev")], start=1):
            ws.cell(row=ws_row + k, column=1, value=lab).font = MAE_FONT
            ws.cell(row=ws_row + k, column=1).alignment = CENTER
            ws.cell(row=ws_row + k, column=1).border = BORDER
            cc = ws.cell(row=ws_row + k, column=base + 2, value=float(mt[key]))
            cc.number_format = "0.0000"; cc.font = MAE_FONT; cc.alignment = CENTER; cc.border = BORDER
        ws_cell.number_format = "0.0000"
        ws_cell.font = MAE_FONT
        ws_cell.alignment = CENTER
        ws_cell.border = BORDER
        ws_cell.fill = PatternFill("solid", fgColor="E2EFDA")
        ws_cell_refs[method] = f"'{sheet_name}'!{diff_letter}{ws_row}"

    next_free_row = ws_row + 5

    # ---- Best/worst criterion match ------------------------------------------ #
    if ex.reported_bw:
        t = next_free_row
        ws.cell(row=t, column=1,
                value=(f"Best/worst match  (derived: best={ex.best_criteria[0]}, "
                       f"worst={ex.worst_criteria[0]})")).font = Font(name=FONT_NAME, size=10, bold=True, color=NAVY)
        ws.merge_cells(start_row=t, start_column=1, end_row=t, end_column=6)
        t += 1
        bw_headers = ["Source", "Reported best", "Reported worst", "Derived best", "Derived worst", "Result"]
        for j, h in enumerate(bw_headers):
            c = ws.cell(row=t, column=1 + j, value=h)
            c.font = Font(name=FONT_NAME, size=9.5, bold=True, italic=True, color=NAVY)
            c.fill = PatternFill("solid", fgColor=LIGHT_BLUE)
            c.alignment = CENTER
            c.border = BORDER
        t += 1
        for m in ex.bw_matches:
            vals = [m["source"], m["reported_best"], m["reported_worst"],
                    m["auto_best"], m["auto_worst"],
                    "Full match" if m["full_match"] else ("Partial" if (m["best_ok"] or m["worst_ok"]) else "Mismatch")]
            for j, v in enumerate(vals):
                c = ws.cell(row=t, column=1 + j, value=v)
                c.font = BODY_FONT
                c.alignment = CENTER
                c.border = BORDER
            fill_color = "C6E0B4" if m["full_match"] else ("FFEB9C" if (m["best_ok"] or m["worst_ok"]) else "F8696B")
            ws.cell(row=t, column=6).fill = PatternFill("solid", fgColor=fill_color)
            t += 1
        next_free_row = t + 1

    # ---- Individual expert CRs (Case 6) ----------------------------------- #
    if getattr(ex, "expert_cr", None):
        t = next_free_row
        ws.cell(row=t, column=1, value="Individual expert matrices (Table 4): CR recomputed vs. reported (Table 5)"
                ).font = Font(name=FONT_NAME, size=10, bold=True, color=NAVY)
        ws.merge_cells(start_row=t, start_column=1, end_row=t, end_column=6)
        t += 1
        for j, h in enumerate(["Expert", "CR (recomputed)", "CR (reported)", "Difference"]):
            c = ws.cell(row=t, column=1 + j, value=h)
            c.font = Font(name=FONT_NAME, size=9.5, bold=True, italic=True, color=NAVY)
            c.fill = PatternFill("solid", fgColor=LIGHT_BLUE)
            c.alignment = CENTER
            c.border = BORDER
        t += 1
        for k, (cr, cr_rep) in enumerate(zip(ex.expert_cr, ex.expert_cr_reported), start=1):
            for j, v in enumerate([f"E{k}", round(cr, 4), cr_rep, round(cr - cr_rep, 4)]):
                c = ws.cell(row=t, column=1 + j, value=v)
                c.font = BODY_FONT; c.alignment = CENTER; c.border = BORDER
                if j > 0:
                    c.number_format = "0.000"
            t += 1
        ws.cell(row=t, column=1, value=("The reported CRs are not reproduced by the printed matrices "
                                         "(eigenvector, normalised-column and geometric-mean approximations agree); "
                                         "all recomputed CRs are below 0.10.")).font = Font(name=FONT_NAME, size=9, italic=True, color="595959")
        ws.merge_cells(start_row=t, start_column=1, end_row=t, end_column=8)
        next_free_row = t + 2

    # Native Excel chart
    chart = BarChart()
    chart.type = "col"
    chart.grouping = "clustered"
    chart.title = f"{sheet_name}: Derived weights"
    chart.y_axis.title = "Weight"
    chart.x_axis.title = "Criterion"
    chart.height = 9
    chart.width = 24
    chart.style = 10

    cats = Reference(ws, min_col=1, min_row=first_data_row, max_row=last_data_row)
    for method in methods:
        base = method_col_start[method]
        data_calc = Reference(ws, min_col=base, min_row=first_data_row, max_row=last_data_row)
        chart.series.append(Series(data_calc, title=method))
    chart.set_categories(cats)
    chart.gapWidth = 60
    chart.overlap = -10
    ws.add_chart(chart, f"A{next_free_row + 1}")

    ws.freeze_panes = f"B{first_data_row}"
    return {"corr": corr_cell_refs, "ws": ws_cell_refs}, sheet_name


refs_by_sheet = {}
for ex in EXAMPLES:
    sheet_name = SHEET_NAME_MAP[id(ex)]
    refs, sname = build_example_sheet(wb, ex, sheet_name)
    refs_by_sheet[sname] = refs

row = r0 + 1
for ex, sheet_name, method in summary_rows:
    corr_ref = refs_by_sheet[sheet_name]["corr"][method]
    ws_ref = refs_by_sheet[sheet_name]["ws"][method]
    ws0.cell(row=row, column=4, value=f"={corr_ref}")
    ws0.cell(row=row, column=5,
             value=(f'=IF(D{row}>=0.99,"Excellent",IF(D{row}>=0.95,"Very high",'
                    f'IF(D{row}>=0.80,"High",IF(D{row}>=0.50,"Moderate","Weak"))))')).font = BODY_FONT
    ws0.cell(row=row, column=7, value=f"={ws_ref}").font = BODY_FONT
    ws0.cell(row=row, column=7).number_format = "0.0000"
    ws0.cell(row=row, column=7).alignment = CENTER
    ws0.cell(row=row, column=7).border = BORDER
    row += 1

first_row, last_row = r0 + 1, row - 1
ws0.conditional_formatting.add(
    f"D{first_row}:D{last_row}",
    ColorScaleRule(start_type="num", start_value=0, start_color="F8696B",
                    mid_type="num", mid_value=0.8, mid_color="FFEB9C",
                    end_type="num", end_value=1, end_color="C6E0B4"),
)
ws0.conditional_formatting.add(
    f"E{first_row}:E{last_row}",
    CellIsRule(operator="equal", formula=['"Excellent"'], fill=PatternFill("solid", fgColor="C6E0B4")),
)
ws0.conditional_formatting.add(
    f"E{first_row}:E{last_row}",
    CellIsRule(operator="equal", formula=['"Weak"'], fill=PatternFill("solid", fgColor="F8696B")),
)
ws0.conditional_formatting.add(
    f"G{first_row}:G{last_row}",
    ColorScaleRule(start_type="num", start_value=0, start_color="F8696B",
                    mid_type="num", mid_value=0.8, mid_color="FFEB9C",
                    end_type="num", end_value=1, end_color="C6E0B4"),
)


# ============================================================================ #
# 5) MANUSCRIPT TABLES (Tables 1-5)
# ============================================================================ #
def fmt_vec(v):
    return "(" + ", ".join(f"{x:g}" if float(x).is_integer() else f"{x:.3f}" for x in np.round(v, 3)) + ")"

T1, T2, T3, T4, T5, T6 = [], [], [], [], [], []
for k, ex in enumerate(EXAMPLES, start=1):
    T1.append([k, ex.n, ", ".join(available_methods(ex)), round(ex.ahp_cr, 4)])
    T2.append([k, ", ".join(ex.best_criteria), ", ".join(ex.worst_criteria), fmt_vec(ex.bo_vector),
               fmt_vec(ex.ow_vector), ex.fucom_ranking_str, round(ex.bwm_xi_linear, 4),
               round(ex.bwm_xi_nonlinear, 4), round(ex.fucom_chi, 4)])
    calc_map = {"AHP": ex.ahp_weights, "BWM": ex.bwm_weights, "FUCOM": ex.fucom_weights}
    for m in available_methods(ex):
        mt = all_metrics(calc_map[m], ex.reported_weights[m])
        T3.append([k, m, round(mt["r"], 4), round(mt["rho"], 4), round(mt["ws"], 4),
                   round(mt["mae"], 4), round(mt["maxdev"], 4)])
    T4.append([k, round(pearsonr(ex.ahp_weights, ex.bwm_weights)[0], 4),
               round(pearsonr(ex.ahp_weights, ex.fucom_weights)[0], 4),
               round(pearsonr(ex.bwm_weights, ex.fucom_weights)[0], 4)])
    if "BWM" in ex.reported_weights:
        for model, wv in (("linear", ex.bwm_weights_linear), ("non-linear", ex.bwm_weights_nonlinear)):
            mt = all_metrics(wv, ex.reported_weights["BWM"])
            used = "*" if model.replace("-", "") == ex.bwm_model else ""
            T6.append([k, model + used, round(mt["r"], 4), round(mt["rho"], 4), round(mt["ws"], 4),
                       round(mt["mae"], 4), round(mt["maxdev"], 4)])
    for mm in ex.bw_matches:
        res = "Full match" if mm["full_match"] else ("Partial (worst)" if mm["worst_ok"] else
              "Partial (best)" if mm["best_ok"] else "Mismatch")
        T5.append([k, mm["source"], mm["reported_best"], mm["reported_worst"], mm["auto_best"], mm["auto_worst"], res])

TABLES = {
    "Table 1 (CR)": (["Case", "n", "Methods", "CR"], T1),
    "Table 2": (["Case", "Best", "Worst", "a_B", "a_W", "FUCOM order", "xi^L*", "xi* (NL)", "chi*"], T2),
    "Table 3": (["Case", "Method", "r", "rho", "WS", "MAE", "max|dw|"], T3),
    "Table 4": (["Case", "AHP-BWM", "AHP-FUCOM", "BWM-FUCOM"], T4),
    "Table 5": (["Case", "Source", "Rep. best", "Rep. worst", "Der. best", "Der. worst", "Result"], T5),
    "Table 6 (BWM variant sensitivity; * = variant used in the source)":
        (["Case", "BWM model", "r", "rho", "WS", "MAE", "max|dw|"], T6),
}
for name, (hdr, rows) in TABLES.items():
    print("\n" + name)
    print(pd.DataFrame(rows, columns=hdr).to_string(index=False))

wst = wb.create_sheet("Manuscript_Tables", 1)
wst.sheet_view.showGridLines = False
rr = 1
for name, (hdr, rows) in TABLES.items():
    wst.cell(row=rr, column=1, value=name).font = TITLE_FONT
    rr += 1
    for j, h in enumerate(hdr, start=1):
        c = wst.cell(row=rr, column=j, value=h)
        c.font = HEADER_FONT; c.fill = HEADER_FILL; c.alignment = CENTER; c.border = BORDER
    rr += 1
    for row_vals in rows:
        for j, v in enumerate(row_vals, start=1):
            c = wst.cell(row=rr, column=j, value=v)
            c.font = BODY_FONT; c.border = BORDER
            if isinstance(v, float):
                c.number_format = "0.0000"
        rr += 1
    rr += 2
for col, wdt in zip("ABCDEFGHI", [8, 34, 14, 40, 40, 46, 12, 12, 12]):
    wst.column_dimensions[col].width = wdt

out_path = "BWM_FUCOM_Derivation_Report.xlsx"
wb.save(out_path)
print(f"Excel report saved: {out_path}")


# ============================================================================ #
# 6) MANUSCRIPT FIGURES (Fig. 2 and Fig. 3)
#    Drawn at their final printed size (text width 15 cm) with 7.5-9 pt text;
#    saved as vector (SVG, PDF) and 600 dpi raster (PNG, TIFF).
# ============================================================================ #
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

FIG_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.linewidth": 0.6,
    "axes.edgecolor": "#808080",
    "axes.labelcolor": "#1a1a1a",
    "text.color": "#1a1a1a",
    "xtick.color": "#1a1a1a",
    "ytick.color": "#1a1a1a",
    "svg.fonttype": "path",
    "pdf.fonttype": 42,
}
CASE_TITLES = {
    1: "Case 1 – Cloud resource allocation",
    2: "Case 2 – Car purchase",
    3: "Case 3 – Smart city components",
    4: "Case 4 – Beekeeping land suitability",
    5: "Case 5 – Libyan airlines",
    6: "Case 6 – Wind farm site selection",
    7: "Case 7 – Flash flood susceptibility",
}
CM = 1 / 2.54


def save_all_formats(fig, stem):
    fig.savefig(f"{stem}.svg", facecolor="white")
    fig.savefig(f"{stem}.pdf", facecolor="white")
    fig.savefig(f"{stem}.png", dpi=600, facecolor="white")
    fig.savefig(f"{stem}.tiff", dpi=600, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})


def manuscript_figure_2(stem="Figure_2"):
    with plt.rc_context(FIG_RC):
        fig, axes = plt.subplots(4, 2, figsize=(15.0 * CM, 19.0 * CM))
        for k, ex in enumerate(EXAMPLES, start=1):
            ax = axes.flat[k - 1]
            methods = available_methods(ex)
            calc_map = {"AHP": ex.ahp_weights, "BWM": ex.bwm_weights, "FUCOM": ex.fucom_weights}
            m = len(methods)
            x = np.arange(ex.n)
            width = 0.8 / m
            offsets = [-0.4 + width * (i + 0.5) for i in range(m)]
            for method, off in zip(methods, offsets):
                ax.bar(x + off, calc_map[method], width, color=METHOD_COLORS[method],
                       linewidth=0, zorder=3)
                rep = np.asarray(ex.reported_weights[method], dtype=float)
                ax.hlines(rep, x + off - width * 0.42, x + off + width * 0.42,
                          colors=COLOR_REPORTED, linewidth=1.4, zorder=4)
            ax.set_title(CASE_TITLES[k], fontsize=8.5, fontweight="bold", loc="left", pad=4)
            ax.set_xticks(x)
            ax.set_xticklabels(ex.names, fontsize=7.5)
            ax.tick_params(axis="y", labelsize=7.5, length=0, pad=2)
            ax.tick_params(axis="x", length=0, pad=2)
            ax.yaxis.set_major_formatter(mticker.PercentFormatter(xmax=1, decimals=0))
            ax.yaxis.set_major_locator(mticker.MaxNLocator(5))
            ax.set_ylim(0, max(np.max(np.r_[calc_map[mm], ex.reported_weights[mm]]) for mm in methods) * 1.12)
            ax.yaxis.grid(True, color="#d9d9d9", linewidth=0.5, zorder=0)
            ax.set_axisbelow(True)
            for side in ("top", "right", "left"):
                ax.spines[side].set_visible(False)
            if (k - 1) % 2 == 0:
                ax.set_ylabel("Weight", fontsize=8)
        leg_ax = axes.flat[7]
        leg_ax.axis("off")
        handles = [Line2D([0], [0], color=COLOR_REPORTED, linewidth=1.6, label="Reported")] + \
                  [Patch(facecolor=METHOD_COLORS[mm], label=(f"{mm} (recomputed)" if mm == "AHP" else f"{mm} (derived)")) for mm in ("AHP", "BWM", "FUCOM")]
        leg_ax.legend(handles=handles, loc="center", frameon=False, fontsize=9, handlelength=1.8,
                      handleheight=1.0, labelspacing=0.8)
        fig.tight_layout(h_pad=1.2, w_pad=1.0)
        save_all_formats(fig, stem)
        plt.close(fig)


def manuscript_figure_3(stem="Figure_3"):
    markers = {"AHP": ("o", 22), "BWM": ("s", 20), "FUCOM": ("^", 24)}
    with plt.rc_context(FIG_RC):
        fig, ax = plt.subplots(figsize=(11.5 * CM, 10.5 * CM))
        for method in ("AHP", "BWM", "FUCOM"):
            xs, ys = [], []
            for ex in EXAMPLES:
                if method in ex.reported_weights:
                    calc = {"AHP": ex.ahp_weights, "BWM": ex.bwm_weights, "FUCOM": ex.fucom_weights}[method]
                    xs.extend(ex.reported_weights[method])
                    ys.extend(calc)
            mk, sz = markers[method]
            ax.scatter(xs, ys, marker=mk, s=sz, color=METHOD_COLORS[method], edgecolors="white",
                       linewidths=0.4, label=method, zorder=3)
        ax.plot([0, 0.55], [0, 0.55], linestyle="--", color="#555555", linewidth=0.9, zorder=2)
        ax.set_xlim(0, 0.55)
        ax.set_ylim(0, 0.55)
        ax.set_aspect("equal")
        ax.set_xlabel("Reported weight (original study)", fontsize=9)
        ax.set_ylabel("Derived weight (this study)", fontsize=9)
        ax.tick_params(labelsize=8.5, length=2.5, width=0.6)
        ax.grid(True, color="#e3e3e3", linewidth=0.5, zorder=0)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.legend(loc="upper left", frameon=False, fontsize=9, handletextpad=0.3, borderaxespad=0.6)
        fig.tight_layout()
        save_all_formats(fig, stem)
        plt.close(fig)


manuscript_figure_2()
manuscript_figure_3()
print("Manuscript figures saved: Figure_2 and Figure_3 (.svg, .pdf, .png, .tiff)")


# ============================================================================ #
# 7) FIG. 1 -- framework diagram (final size 15 cm x 9.0 cm)
# ============================================================================ #
from matplotlib.patches import FancyBboxPatch


def manuscript_figure_1(stem="Figure_1"):
    rc = {"font.family": "serif",
          "font.serif": ["Times New Roman", "Liberation Serif", "DejaVu Serif"],
          "mathtext.fontset": "stix", "font.size": 9,
          "svg.fonttype": "path", "pdf.fonttype": 42}
    W, H = 15.0, 9.0                                   # cm
    with plt.rc_context(rc):
        fig = plt.figure(figsize=(W * CM, H * CM))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, W)
        ax.set_ylim(0, H)
        ax.axis("off")

        def box(x, y, w, h, title, lines, shaded=False, dashed=False):
            ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="square,pad=0",
                                        facecolor="#e6e6e6" if shaded else "white",
                                        edgecolor="black", linewidth=0.8,
                                        linestyle=(0, (3, 2)) if dashed else "solid"))
            ax.text(x + w / 2, y + h - 0.42, title, ha="center", va="center",
                    fontsize=9.5, fontweight="bold")
            for i, ln in enumerate(lines):
                ax.text(x + w / 2, y + h - 1.0 - i * 0.5, ln, ha="center", va="center", fontsize=9)

        def arrow(x0, y0, x1, y1, dashed=False):
            ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                        arrowprops=dict(arrowstyle="-|>", color="black", linewidth=0.8,
                                        mutation_scale=8, shrinkA=0, shrinkB=0,
                                        linestyle=(0, (3, 2)) if dashed else "solid"))

        bw, bh = 2.6, 2.3
        ymid = 3.35
        xs = [0.05, 2.95, 5.85]
        box(xs[0], ymid, bw, bh, "Step 1. Input",
            ["Published AHP", r"matrix $A=[a_{ij}]$", "check: Eq. (1)"], shaded=True)
        box(xs[1], ymid, bw, bh, "Step 2. AHP",
            [r"Eigenvector $\mathbf{w}^{AHP}$", r"$\lambda_{\max}$, CR", "Eq. (2)"])
        box(xs[2], ymid, bw, bh, "Step 3. B and W",
            [r"$B=\mathrm{arg\,max}\; w_j^{AHP}$", r"$W=\mathrm{arg\,min}\; w_j^{AHP}$", "Eq. (3)"])
        x4, w4, h4 = 9.0, 2.7, 2.8
        box(x4, 6.0, w4, h4, "Step 4a. BWM",
            [r"$\mathbf{a}_B$ = row $B$ of $A$", r"$\mathbf{a}_W$ = column $W$", "Eqs. (4)–(5)",
             r"$\rightarrow \mathbf{w}^{BWM},\ \xi^{L*}$"])
        box(x4, 0.8, w4, h4, "Step 4b. FUCOM",
            [r"order by $a_{Bj}$", r"$\varphi_{k/(k+1)}$: Eq. (6)", "Eqs. (7)–(8)",
             r"$\rightarrow \mathbf{w}^{FUCOM},\ \chi^{*}=0$"])
        x5 = 12.2
        box(x5, ymid, 2.75, bh, "Step 5. Validation",
            [r"$r$, $\rho$, WS, MAE", "best/worst", "match"], shaded=True)
        box(x5, 6.8, 2.75, 1.8, "Source study",
            [r"$\mathbf{w}^{rep}$; reported", r"$B$ and $W$"], dashed=True)

        yc = ymid + bh / 2
        arrow(xs[0] + bw, yc, xs[1], yc)
        arrow(xs[1] + bw, yc, xs[2], yc)
        arrow(xs[2] + bw, yc + 0.5, x4, 6.0 + h4 / 2)
        arrow(xs[2] + bw, yc - 0.5, x4, 0.8 + h4 / 2)
        arrow(x4 + w4, 6.0 + h4 / 2, x5, yc + 0.5)
        arrow(x4 + w4, 0.8 + h4 / 2, x5, yc - 0.5)
        arrow(x5 + 2.75 / 2, 6.8, x5 + 2.75 / 2, ymid + bh, dashed=True)
        # AHP weights passed directly to Step 5
        xa = xs[1] + bw / 2
        ax.plot([xa, xa, x5 + 2.75 / 2], [ymid, 0.3, 0.3], color="black", linewidth=0.8)
        arrow(x5 + 2.75 / 2, 0.3, x5 + 2.75 / 2, ymid)
        ax.text((xa + xs[2] + bw) / 2, 0.55, r"$\mathbf{w}^{AHP}$", ha="center", va="center", fontsize=9)

        save_all_formats(fig, stem)
        plt.close(fig)


manuscript_figure_1()
print("Manuscript figure saved: Figure_1 (.svg, .pdf, .png, .tiff)")