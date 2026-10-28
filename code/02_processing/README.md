# 02_processing — tables from the simulations

Reads the per-tumour files in `data/raw/<block>/` and writes tables to `data/processed/<block>/`.
The bash scripts in `code/01_datageneration/` call these scripts automatically after the
simulations; they can be re-run on their own at any time (minutes), e.g. with another `k_min`:

```bash
julia -p 8 --project=. code/02_processing/02_statistics.jl 01_dynamics --kmin-frac 0.05 --suffix _kmin0.05
```

Time is in cell cycles (1 cycle ≈ 4 days); `t_rel` is the time since detection (negative during
growth). Resistant = clones with α < 1. Functions live in `code/X_helpers/processing.jl`.

## Tables per block

### `outcomes.csv` — one row per tumour × arm (`01_collect.jl`)
`block, cohort, tumour, arm` · `outcome` (progressed / extinct / censored; the single round ends
`regrown` at N_det and is left out of survival summaries) · `ttp` (time to
progression; missing unless progressed) · `t_end` · `cum_dose` (∫ m dt) · `n_cycles` (completed
holiday cycles) · `n_blocks` · `nu, alpha, fK, cost` (tumour parameters) · `kind, level, tmax` (arm;
`tmax` = its follow-up) ·
`det_N, det_res_frac, det_n_res, det_largest, det_B, det_mean_burden` (at detection).

### `counts.csv.gz` — time series, one row per arm × time (every 1, scan arms every 2)
`block, cohort, tumour, arm, t_rel, N, m` · `nS, nR` (sensitive / resistant cells) · `n_res`
(resistant clones) · `res_pre, res_new` (resistant lineages — resistant clones whose parent is
sensitive — arisen before / after detection; second hits inside resistant lineages are clones but
not new lineages) · `largest`
(frequency of the largest resistant clone) · `mean_burden`.

### `snapshots.csv` — entire tree at every stored time point (`02_statistics.jl`)
`block, cohort, tumour, arm, tag` (growth / detection / blockend / grid / early / end) · `idx`
(block number or grid time) · `t_rel, N, m` · `B, B_pre, B_new` (total mutational burden; on
lineages born before / after detection) · `clonal` · resistant clones: `res_frac, res_D0, res_D1,
res_D2, res_Dinf` (Hill numbers of order 0, 1, 2, ∞), `res_pre, res_new, largest` · burden
`mean_burden, var_burden` · passengers from the SFS: `S` (number of mutations), `D1n, D2n, Dinfn,
H1n` (Hill numbers / Shannon entropy, normalised frequencies), `H1cf, D2cf` (cell-fraction
frequencies, v1 eq. 10), `gamma, gamma_se, n_tail` (tail exponent at `k_min = max(2, kmin_frac·N)`,
over k_min ≤ k < N: clonal mutations left out).

### `samples.csv` — subsamples at every stored time point
Same identifiers plus `n` (100, 1000, 10⁴; only if ≤ N) and `fmin` (0, 0.01, 0.05: mutations
in fewer than `fmin·n` cells removed; the tail fit starts at k_min ≥ fmin·n) · SFS statistics as
above · `mean_burden, var_burden` ·
`res_D0, res_D2, res_Dinf` of the subsample.

## Summaries (`03`–`08`)

| script | output | content |
|---|---|---|
| `03_sfs.jl` | `sfs_mean.csv` | `cohort, arm, tag, idx, k, S_mean, S_q25, S_q75, n_tumours`: log-binned SFS (mean per integer k in each bin, fixed bin edges) averaged over tumours; single round: only tumours that respond (reach 0.1 N_det) |
| `04_examples.jl` | `example_traj.csv`, `example_clones.csv` | example tumour of `full` (closest to the median resistant fraction among tumours 1–10): `N, nS, nR, m` over time for `mtd, hol, mod, met`; sizes of resistant clones over time |
| `05_tradeoff.jl` | `tradeoff.csv`, `km.csv` | per cohort × arm: `n, n_progressed, n_censored, ttp_km_median, ttp_rmst` (restricted mean up to the follow-up), resistant clones `resD0_*` and `resD2_med`, burden `B_*` at the end of the arm, `resD0_t50_med` at t = 50, paired gain in restricted TTP vs MTD `gain_med, gain_lo, gain_hi, gain_p` (Wilcoxon); Kaplan–Meier curves. Block 02 is checked against block 01 (identical tumours must give identical outcomes) |
| `06_trial.jl` | `trial.csv`, `trial_cor.csv`, `trial_policy.csv`, `trial_km.csv` | per patient: parameters, restricted TTP / event / dose per arm, `benefit_hol`, `benefit_mod`, diversity at detection (oracle = entire tree; sample = 1000 cells, f_min 1%). Correlations (Spearman, bootstrap CI): `setting` = within (block 01 `full`, one parameter set), across (trial), partial (trial, controlling for log ν, α, f_K); `time` = detection or early (end of the first full-dose phase); all sample sizes and thresholds. Policies: always MTD / holidays / modulation and "adaptive if statistic below/above threshold", fitted on one half, evaluated on the other |
| `08_powerlaw.jl` | `powerlaw.csv` | `cohort, tumour, arm, tag, idx, N, kmin, kmin_frac, gamma, se, n_tail, ks_kmin`: tail exponent of the entire-tree SFS over k_min = 10⁻⁴–10⁻¹ N (8 per decade) and the KS-selected k_min, at detection and after the last stored holiday cycle (cohorts `neutral`, `full`) |
| `07_validation.jl` | `validation_sfs.csv`, `validation_ld.csv` | mean SFS (k ≤ 200) vs $2\mu N\sum_l (d/b)^l/((k+l)(k+l+1))$; resistant clones and cells at detection |

Restricted time to progression: progression time, or the arm's follow-up `tmax` if the tumour did
not progress — censored, or extinct (cured: progression-free until `tmax`, not lost at the time of
extinction); `event` marks progression.
