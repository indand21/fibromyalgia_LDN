# Curated data for the LDN–fibromyalgia QSP model

Every CSV carries `source`, `location` and `verified` columns and is validated by
`tests/testthat/test-data.R` through `read_dataset()` in `R/data_io.R`.

`verified` values: `yes` = the number was read in the primary source; `secondary` = taken
from a review or digitized from a figure; `digitised` = read off a figure of the primary source
by a stated method (R4 fix round 1, Bested 2023 Figure 4; method in `location`); `no` = uncertain.
Missing values are `NA`.
No field contains the `#` character (it would truncate a row, because `read_dataset()`
calls `read.csv(comment.char = "#")`), and every row has a non-empty `location`.

Arm naming is `placebo` for placebo arms and `ldn_<dose>mg` for active arms. `schedule`
strings use `startday:dose;startday:dose` with strictly increasing start days beginning at 0.

---

## training/pk_summary.csv (22 rows; R4 and R4 fix round 1: 10 Bested 2023 low-dose rows)

Oral naltrexone (`NTX`) and 6β-naltrexol (`BN`) plasma PK summary statistics.
**Verebey 1976 (12 rows, 100 mg):** Cmax, Tmax, terminal half-life and 24 h concentrations after an
acute dose, and Cmax and half-life after chronic daily dosing, n = 4. All rows are `verified = yes`,
read in the paper's abstract, and every row carries a numeric `dose_mg`. `sd` is the between-subject SD
(mean ± SD of 4 subjects), never an SE.
**Bested 2023 (10 rows, 4.5 mg; R4, corrected in R4 fix round 1):** plasma NTX and 6β-naltrexol from a
randomised placebo-controlled crossover trial in fibromyalgia (PMID 38226027, Pain Rep 8(4):e1080;
PMC10789452). Verified against the full text and Supplemental Digital Content 4 (R4 fix round 1).
- Day-1 profile (8 rows, regimen `single`): stats `C0.25`, `C0.5`, `C0.75` and `C1` (concentration 15,
  30, 45 and 60 min after the first 4.5 mg dose) for each analyte. Figure 4 plots means with shaded SD
  bands and 95% CI error bars and prints N at each time (34, 37, 29 and 9-10). Means and SDs were
  digitised (`verified = digitised`): mean = midpoint of the plotted CI; SD = CI half-width × √N /
  t(0.975, N − 1), cross-checked against the SD band (agreement within about 10%). Where the text states
  a value it is used instead: NTX 0.66 and BN 4.61 µg/L at 30 min (Figure 4 caption). **Conflict:** the
  plotted BN mean at 30 min is 5.61 µg/L (the CI 4.09-7.14 is not centred on 4.61); the difference of
  exactly 1.00 suggests a typographical error in the caption, but the stated value is kept. The
  between-subject CVs are 0.6-1.5 (0.87 and 0.99 at 30 min). These are repeated measures on the same
  patients, so each row has `obj_weight = 0.25` and the profile counts as about one observation per
  analyte. n is the N printed at each time (37 at 30 min, not the 52 per-protocol patients).
- Day-21 (2 rows, regimen `daily`, stat `Css13h`): medians 0.33 (NTX) and 5.29 (BN) µg/L from samples
  taken 1-3 h after tablet intake (days 14 and 21; Supplemental Digital Content 4). No spread and no n
  are reported (`sd` and `n` NA). `cv_assumed = 0.8` is an assumption justified by the day-1 CVs, with a
  sensitivity at 0.5 and 1.2 in `scripts/01_fit_pk.R`. `Css13h` is the model mean over the 1-3 h window
  at steady state (the window is averaged because the time within it is not reported). A median is
  compared with a model mean: it is a typical value, and under a log-normal with CV 0.8 the mean is
  about 1.28 × the median, whereas the Verebey rows are arithmetic means.
- Assay: UPLC-MS/MS with LOQ 0.015 µg/L (signal to noise > 10) and intermediate precision 9.5% at
  0.35 µg/L (Supplemental Digital Content 4, http://links.lww.com/PR9/A195). Every NTX value is at
  least 20 times the LOQ, so near-LOQ bias is not a concern.
- Population: 46 of 52 patients were women (male/female 6/46), mean weight 86.5 kg, mean BMI 31.2;
  4.5 mg tablets compounded by a pharmacy (magistral production, Glostrup Apotek, Copenhagen).
  The trial is a crossover, so the LDN samples come from whichever period a patient received LDN.
  Verebey studied men at 100 mg: body size, sex and formulation compete with dose as explanations of
  any low-dose difference, and the model cannot separate them.
- Weighting (`R/pk.R`, `pk_obs_cv`, `pk_obs_weight`): CV = `sd / mean` where a dispersion is reported,
  else `cv_assumed`, else the CV floor (0.10; Tmax 0.25); each row contributes `obj_weight` × (log
  residual / CV)². The two optional columns `cv_assumed` and `obj_weight` (R4 fix round 1) are `NA` and
  `1` on every Verebey row, which are therefore weighted exactly as before.
`regimen` (added 2026-10-01, no value changed) is `single` for the acute rows and `daily` for the "Verebey 1976 chronic" rows: 100 mg orally once daily at steady state, per the abstract (PMID 954353, DOI 10.1002/cpt1976203315: levels 24 h after the daily doses showed steady state rapidly attained); the duration of chronic dosing is not stated in the abstract.
`C24` rows (added 2026-10-01): concentration 24 h after the acute 100 mg dose, naltrexone 2.1 ± 0.47 and 6β-naltrexol 17.6 ± 5.0 ng/mL (abstract sentence quoted in `location`; ± taken as SD, as for the Cmax rows). Not curated as numbers, and recorded here instead: (i) the chronic 24 h levels, which the abstract describes only qualitatively ("by 24 hr both compounds declined to levels of the same order as in the acute state at 24 hr"; "no accumulation of naltrexone and beta naltrexol in the plasma after chronic treatment"); (ii) "The plasma levels of naltrexone declined slowly from 24 through 72 hr from 2.4 to 1.7 ng/ml, with an apparent half-life of 96 hr", because the sentence does not say whether it follows the acute or the last chronic dose. Its two numbers are held separately in `training/pk_summary_sensitivity_terminal.csv` (NTX `C24` = 2.4 and `C72` = 1.7 ng/mL, regimen `daily` on the most natural reading, sd NA), used ONLY for the terminal-phase sensitivity fit in `scripts/01_fit_pk.R`, never for the primary PK fit.

**Dose-independent PK values held here as context, not as rows.** The statistics below were
verified in their primary sources but are either pooled across doses or stated without a
dose, so attaching them to a `dose_mg` would be an invention. They are recorded here so that
they remain traceable and usable as priors, and they are deliberately absent from the CSV
because the downstream PK code requires a dose on every row.

| Value | Source | Location |
|---|---|---|
| Naltrexone t½ ≈ 4 h; 6β-naltrexol t½ ≈ 12 h | Meyer 1984 (PMID 6469932) | Abstract; pooled over the 50/100/200 mg tablet arms and the syrup reference, n = 24 |
| Naltrexone terminal t½ 8.9 h after oral dosing (2.7 h after i.v.) | Wall 1981 (PMID 6114837) | Abstract; the abstract states no dose, and the full text was not obtainable, so no dose is asserted |
| Naltrexone t½ 4 h; 6β-naltrexol t½ 13 h | ReVia FDA label 2013 (NDA 018932) | CLINICAL PHARMACOLOGY, Pharmacokinetics; a dose-independent label statement |
| Peak plasma levels of both analytes within 1 h of dosing | ReVia FDA label 2013 (NDA 018932) | CLINICAL PHARMACOLOGY, Absorption; a dose-independent label statement |

## training/efficacy_training.csv (4 rows)

The two Younger crossover trials, which are the estimation (training) set. The `scale`
column is `0-100` on every row, because both trials use a 0–100 VAS while `sd_res_nrs` in
`priors/placebo_priors.csv` is on a 0–10 NRS; anything that combines them must rescale.

* **Younger 2013** (n = 28 analysed of 31 randomised) is a randomised, double-blind,
  counterbalanced crossover with **no washout**: a 2-week baseline with no capsules, then
  placebo (4 weeks) and low-dose naltrexone 4.5 mg (12 weeks) in randomised order, then a
  1-month follow-up. The endpoint of each condition is the mean of the **final 3 days** of
  that condition, expressed as percent reduction from the 14-day baseline mean. `day` is
  therefore the day within each condition: 28 for placebo and 84 for LDN. `baseline_nrs`
  is 50.8 on the trial's **0–100** pain VAS.
* **Younger 2009** (n = 10 completers) is a single-blind, fixed-order crossover: baseline
  2 weeks, placebo 2 weeks, drug 8 weeks, washout 2 weeks. `day` is 14 for placebo and 56
  for LDN. Its primary outcome is **overall fibromyalgia symptom severity** on a 0–100 VAS,
  not a pain-specific score; the 2.3 % and 32.5 % reductions are that outcome. The paper
  reports only F and p for daily pain, so `baseline_nrs`, `se` and `sd_pct` are `NA`.

**Dispersion for Younger 2013.** The paper states "18.0 ± 10.8 %" and "28.8 ± 9.3 %" as
**mean ± 95 % confidence interval**, not ± SD. The interval came from n = 28, so `se` is the
half-width divided by **t(27) = 2.052** (5.26 for placebo, 4.53 for LDN) and `sd_pct` is
se × √28 (27.9 and 24.0). *Caveat:* in a counterbalanced crossover the reported interval may
be a within-subject interval rather than a between-subject one. If it is, `sd_pct` overstates
the between-subject SD, and any model that treats it as between-subject variability will be
conservative.
*Conflict:* Bested et al. 2023 (Pain Reports crossover trial; sample-size section 2.12.2) quotes the
same two statistics, as FIQ pain-reduction figures, as "28.8 (12.5) %" and "18.0 (14.6) %" with the
bracketed values called SDs. Those SDs are not
reproducible from Younger's reported CI and n, so the primary CI was used.
*Sensitivity (R3a, review I10):* `scripts/05_calibrate_variability.R` repeats the variability
calibration with the quoted SDs (14.6 placebo, 12.5 LDN) in memory, without changing this file, and
writes `output/results/omegas_sd_sensitivity.csv`; `scripts/07_virtual_trials.R` carries those
omegas into `virtual_summary_sd_sensitivity.csv` and `power_curve_sd_sensitivity.csv`.

## params/binding_constants.csv (13 rows)

`analyte` ∈ {NTX, BN} identifies the molecule and `species` holds the organism the receptor
came from; the two together are what distinguishes, for example, the Wentland naltrexone
rows from the Pelotte 6β-naltrexol rows. `quantity` ∈ {Ki, koff, IC50, fu}.

(analyte, species, target, quantity) is **not** a unique key: three constants are reported by
more than one source on purpose — human MOR Ki for naltrexone appears three times (Wentland
0.11 nM, Cassel 1.1 nM and 0.86 nM) and the mouse TLR4 IC50 twice (Wang and Selfridge, the
same 105.5 µM measured in the same assay). Code that needs one value per constant must
select or pool by `source` rather than assume a single row.

* MOR/DOR/KOR Ki for naltrexone (0.11, 60, 0.19 nM) from Wentland 2009 Table 1 (compound 3,
  CHO membranes expressing the human receptors) and for 6β-naltrexol (2.12, 213, 7.42 nM)
  from Pelotte 2009 Table 1 (human receptors in CHO cells).
* Two further human MOR Ki values for naltrexone (1.1 and 0.86 nM) and the MOR dissociation
  rate from Cassel 2005 Tables 1 and 2. **koff is stored as 21.0 h⁻¹**, because the model
  integrates in hours and nothing in the pipeline converts units; Cassel reports it as
  0.35 min⁻¹ with a dissociation half-life of 2.0 min, and both the reported value and the
  conversion are stated in that row's `location`.
* TLR4 IC50 of (+)-naltrexone, 105.5 µM (95 % CI 95.4–115.6 µM) for inhibition of
  LPS-induced nitric oxide in BV-2 mouse microglia (Wang 2016), plus the same value and the
  (+)-β-naltrexol value (242.8 µM) from Selfridge 2015 Table 1. These are cell-based
  functional IC50 values, not TLR4/MD2 binding constants; no Kd of naltrexone at MD2 was found.
* `fu` = 0.79 for naltrexone in human plasma. The FDA label states 21 % plasma protein
  binding; the unbound fraction is one minus that, as recorded in `location`.

## params/dose_response.csv (2 rows)

The Bruun-Plesner 2020 up-and-down estimates ED50 3.88 mg and ED95 5.40 mg. The response
definition is a reduction in self-reported fibromyalgia symptoms judged over a fixed dosing
period in a single-blind up-and-down design with 0.75 mg steps over 0.75–6 mg, 25 evaluable
women; the paper gives no confidence intervals for either dose.

## test/pet_ed50.csv (1 row)

The Rabiner 2011 MOR-occupancy ED50 of 5.60 mg (95 % CI 3.65–7.54), same columns as
`params/dose_response.csv`. It summarises the held-out PET dataset, so it sits under
`test/` behind the stage gate rather than with the estimation parameters.

## params/baseline_bounds.csv (3 rows)

Added in R2 fix round 1 (Younger floor set in fix round 2) (an addition with provenance; no existing value changed).
The `final` and `innova` rows are eligibility criteria read in the sources (`verified = yes`). The
`younger2013` floor is **not** an eligibility criterion but a modelling choice made by the stability rule
below, so its row is `verified = secondary` and its `location` says so (R3a, R2 re-review; no bound value changed). The per-trial
bounds of the simulated **baseline** pain distribution in `simulate_individuals()`: the lower bound
is the trial's documented minimum-pain eligibility criterion, 0 where none applies; the upper bound
is always 10. The reported baseline mean and SD describe the ENROLLED population, which has already
passed the screen, so the code moment-matches a parent normal whose truncation to these bounds has
the reported moments (it does not truncate N(reported mean, reported SD), which would apply the
screen twice).

* **younger2013**: [0.5, 10]. The baseline is the 14-day mean of a 0-100 diary VAS rescaled by 10,
  which is not the instrument any NRS eligibility threshold was applied to, so no screen applies.
  The floor is **0.5**, chosen by a stability rule rather than asserted (R2 fix round 2): an
  enrolled fibromyalgia population has non-trivial pain, and percentage change is undefined at a
  baseline of zero, so with a floor of 0 a handful of near-zero simulated baselines dominate the SD
  of % reduction that `calibrate_omegas()` moment-matches. Rule: the smallest floor among
  {0.5, 1, 1.5, 2} at which, across 200 seeds, the 5-95% range of both (a) the SD of the
  measurement-noise term 100 e / pain0 of the simulated placebo cohort and (b) the calibrated
  `omega_p` lie within ±10% of their medians. Parent moment-matched on [floor, 10] to 5.08/1.50 in
  every case (H0 posterior, n = 4000 per calibration; `output/results/younger_floor_stability.csv`,
  written by the pipeline script `scripts/05b_younger_floor_study.R`, an optional step after script 05
  that uses `calibrate_omegas()` and reads the noise term from its own simulated cohort; the table below
  was first produced by an equivalent scratchpad script in R2 and is reproduced by 05b):

  | floor | parent (mu / sigma) | noise SD median [5%, 95%] | omega_p median [5%, 95%] | passes |
  |---|---|---|---|---|
  | 0 (reference) | 5.081 / 1.508 | 13.32 [-12.8%, +221%] | 1.922 [-96.1%, +4.1%] | no |
  | **0.5** | 5.077 / 1.515 | 11.77 [-4.1%, +4.3%] | 1.950 [-3.0%, +3.1%] | **yes (chosen)** |
  | 1 | 5.065 / 1.530 | 11.44 [-3.3%, +2.9%] | 1.965 [-2.9%, +2.8%] | yes |
  | 1.5 | 5.035 / 1.564 | 11.21 [-2.8%, +2.5%] | 1.975 [-3.0%, +2.7%] | yes |
  | 2 | 4.951 / 1.642 | 11.01 [-2.5%, +2.4%] | 1.983 [-2.9%, +2.7%] | yes |
  *Baseline SD is ASSUMED:* `efficacy_training.csv` has no `baseline_sd`, and none could be derived:
  the paper (Arthritis Rheum 2013;65:529-538, not in PMC) is not open access and its abstract reports
  no baseline dispersion. `calibrate_omegas()` therefore keeps 1.5 NRS points and labels it
  `assumed 1.5` in `omegas.csv` (`base_sd_source`). The same value is the baseline SD the FINAL
  protocol assumed for its own target population (Sample size {14}: "mean of 6.7 ... with a standard
  deviation of 1.5 NRS points"), which is an analogue from a different trial, not a Younger value.
* **final**: [4, 10]. ClinicalTrials.gov NCT04270877, Inclusion Criteria: "A minimum score of 4 in
  self-reported average pain during the last 7 days on a 0-10 numeric rating scale at baseline".
  The protocol paper (Trials 2021, PMID 34781989) lists the criteria under Eligibility criteria {10},
  but that list is absent from the PMC text extraction, so the registry wording is quoted.
* **innova**: [4, 10]. Protocol (Colomer-Carbonell 2022, BMJ Open 12:e055351, PMID 34992118),
  Abstract: "Eligibility criteria include being adult, having a diagnosis of FMS and experiencing
  pain of 4 or higher on a 10-point numerical rating scale".
  *Limitation (accepted by controller ruling, R2 fix round 2):* the moment-matched parent is
  N(9.95, 3.75), a near-flat distribution over [4, 10] centred near the ceiling, because INNOVA's
  reported baseline SD (1.605) is close to that of a uniform on [4, 10] (1.732). It reproduces the
  reported mean and SD exactly, which is all the data constrain; its shape is otherwise weakly
  determined. An enrolled cohort with mean 7.6 and SD 1.6 genuinely has members at the NRS
  ceiling, and in the virtual trials 0.06-1.9% of simulated subjects end at 10 on average (up to
  4.2% in a single posterior draw, INNOVA placebo); `virtual_pain_end_clip.csv` reports it.

## params/cytokine_params.csv (4 rows)

Human plasma half-lives: TNF-α ≈ 20 min (Selby 1987), IL-6 4.2 h (Weber 1993, terminal
half-life after subcutaneous dosing), and IL-10 as a range (`thalf_lower` 2.7 h,
`thalf_upper` 4.5 h; Huhn 1997 reports the mean terminal half-life as a range across dose
groups, not a single mean). See "not found" for IL-1β, IL-17A and TGF-β.

## priors/placebo_priors.csv (5 rows)

* `pbo_rel_mean` = 0.1114 on the relative-reduction scale. **Derivation:** Häuser 2011
  pooled weighted mean difference in fibromyalgia placebo arms is 7.69 points on a 0–100
  pain scale (72 trials, 9827 patients); Häuser does not report a pooled baseline pain
  value, so the pooled fibromyalgia baseline pain of 6.9 on 0–10 (= 69 on 0–100) from
  Byon 2010 Table II (4 phase 2/3 fibromyalgia trials, 2759 patients) was used as the
  denominator: 7.69 / 69 = 0.1114.
* `pbo_logit_sd` = 0.1206. **Derivation:** the same rescaling applied to Häuser's 95 % CI
  (6.10–9.29 → 0.0884–0.1346 relative reduction), transformed to the logit scale and divided
  by 2 × 1.96. *Caveat:* this is the uncertainty of the pooled mean, so it is a tight prior;
  it does not represent between-trial heterogeneity, which Häuser's abstract does not report.
* `kpbo_per_day` = 0.106 day⁻¹ (SE 0.00833), read directly from Byon 2010 Table III
  (placebo onset rate constant; the table also gives the corresponding half-time of 6.5 days).
  **Scale convention (R2, review C2).** The `scale` column says what `value` and `sd` are on, and
  `log_scale_prior()` in `R/hypotheses.R` reads it. `per_day` means the row stores the rate itself
  and its standard error per day; the sampler works on `log_kpbo`, so the code converts at use to
  mean = log(value) and SD = sd / value (delta method; 0.00833 / 0.106 = 0.0786). A row with
  `scale` = `log_per_day` (or `log`) must instead store **log(value)** in `value` and the SD of
  log(value) in `sd`; such a row is used as it stands. A missing or unknown `scale` stops the fit.
* `pbo_emax_logit` = 1.52 (SE 0.125), **reference only, used by no code** (controller ruling, R2;
  labelled in R3a): the placebo Emax on the logit pain scale, Byon 2010 Table III. It is curated for
  traceability; the placebo ceiling is parameterised by `mu_p` (from `pbo_rel_mean` and
  `pbo_logit_sd`), not by this Emax.
* `sd_res_nrs` = 0.48 NRS points, the within-subject measurement SD of the **0–10** NRS.
  Source: Alghadir 2018, a 24-hour test–retest study in 121 patients with knee osteoarthritis
  (ICC 0.95, SEM 0.48, MDC 1.33). No equivalent test–retest SEM measured in a fibromyalgia
  cohort was found, so this is an out-of-population estimate.

## test/pet_occupancy.csv (11 rows)

Held-out PET data. `time_h` is the time after the last dose, and two extra columns say how
far each row can be trusted:

* `time_basis` — `reported` when the paper states the sampling time, otherwise the assumption
  made and its justification.
* `n_doses` (added 2026-10-01, no value changed; **R4 correction**: Weerts row changed
  4 → 3 with reason recorded in `location`) — the number of doses given before the scan,
  derived from the already-curated `regimen` and `location` text: 1 for every single-dose row.
  **Weerts 2008:** `n_doses = 3` (R4 correction from previous value of 4). The internally
  consistent interpretation is that doses were given on days 15, 16 and 17, and the day-18
  scan was at the 24-hour trough after the day-17 dose (before the day-18 dose). The previous
  value of 4 assumed the day-18 dose also preceded the scan, which is inconsistent with the
  assumed 24-hour trough (`time_h = 24`). The 3- versus 4-dose predictions are computed by
  `scripts/03_test_pet.R` at the adopted fit and written to `pet_test_note.txt` (R4 fix round 1;
  the numbers quoted in R4 came from a scratch script on the pre-R4 fit and were removed). The
  4-dose variants (24 h, and a post-dose scan at 2–12 h on day 18) are explicitly labelled
  sensitivities in `pet_weerts_scan_time_sensitivity.csv`.
* `reference` — `absolute` when the occupancy is an absolute percentage, or
  `normalised to 1 h scan` for Lee 1988, whose percentage blockade is defined against that
  study's own 1-hour scan. Code comparing those rows to the model must divide by the model's
  1-hour occupancy rather than using the raw prediction.

Rows:

* **Lee 1988** (50 mg single oral dose, `verified = yes`, `time_basis = reported`):
  percentage blockade of [11C]carfentanil binding of 91, 80, 46 and 30 % at 48, 72, 120 and
  168 h. The 1-hour scan is 100 % by construction and is not a data row. Lee reports **SEM**,
  not SD (6.4, 6.3, 7.4 and 13.8), so `sd` is `NA` and the SEM is recorded in `location`.
  Blockade half-time 72–108 h.
* **Weerts 2008** (50 mg/day on days 15–19, `verified = yes`): mean MOR inhibition
  94.9 ± 4.9 % (SD), n = 21, scanned on day 18. The occupancy is a primary-source value; the
  interval between the last dose and the scan is not stated, so `time_h = 24` is an assumed
  once-daily trough and `time_basis` says so.
* **Rabiner 2011** (`verified = secondary`, `digitized = yes`): six single-dose points
  digitized from Figure 1d (2 mg → 27 %; 5 mg → 45 % and 58 %; 15 mg → 61 %; 50 mg → 91 %
  and 97 %). Doses were chosen adaptively over 2–50 mg, so the plotted x positions were read
  as the administered doses. `time_h = 4` is an assumption: the paper states only that the
  dose–occupancy fit used scans acquired under 8 h post dose. The panel is small
  (429 × 302 px in PMC), so both axes carry a few percent of digitization error.

## test/trials_test.csv (24 rows)

Held-out parallel-group trials. The `estimand` column says what each row's `mean` is, because
the two trials do not report the same quantity.

* **FINAL** (Due Bruun 2024, n = 99): baseline pain 6.3 (1.3) LDN and 6.2 (1.6) placebo on
  NRS 0–10 (Table 1); 12-week change −1.3 (95 % CI −1.7 to −0.8) versus −0.9 (−1.4 to −0.5);
  between-group −0.34 (−0.95 to 0.27), whose `estimand` is
  *difference in change from baseline*. Responder proportions at 12 weeks are from Table 2
  (intention-to-treat): ≥15 % 0.53 vs 0.42, ≥30 % 0.41 vs 0.26, ≥50 % 0.24 vs 0.14.
  *Conflict:* the Discussion quotes 45 % vs 28 % for the ≥30 % responder rate; that comes from
  a complete-case sensitivity analysis (20/44 and 13/46). The intention-to-treat Table 2
  values were used.
  **Titration:** the trial started at 1.5 mg once daily and added 1.5 mg every seventh day to
  6 mg at week 4, then held the highest tolerated dose for 8 weeks. The `schedule` string is
  `0:1.5;7:3;14:4.5;21:6`, confirmed in the paper, in the protocol (Trials 2021;22:693) and
  in the ClinicalTrials.gov record NCT04270877.
* **INNOVA** (Rodríguez-Freire 2026, 98 randomised, 96 in the intention-to-treat analysis,
  47 LDN and 49 placebo): fixed 4.5 mg/day with **no titration**, so `schedule` is `0:4.5`.
  The paper reports absolute NRS levels and an LMM coefficient B (LDN minus placebo **in NRS
  level**, so positive favours placebo), not change scores. Both forms are held here:
  * `diff_nrs` rows carry B unchanged, with `estimand` = *LMM-adjusted difference in level*:
    0.49 (−0.32 to 1.31) at 3 months, 0.13 (−0.72 to 0.97) at 6 months and 0.37
    (−0.53 to 1.26) at 12 months (days 90, 180, 365).
  * `change_nrs` rows per arm at the same days are **derived** as the Table 3 level minus the
    Table 3 baseline: LDN −0.33, −0.45, −0.28; placebo −0.64, −0.39, −0.37. `se` is `NA`,
    because no standard error of a change is derivable without the within-subject correlation,
    which the paper does not report. Each row's `location` names the two levels used.
  * The **unadjusted** change difference at 3 months is therefore +0.31 (−0.33 minus −0.64),
    against the LMM-adjusted +0.49. The gap reflects the baseline imbalance between the arms
    (7.66 LDN versus 7.47 placebo), which the mixed model adjusts for and the raw difference
    does not.

## test/innova_nrs_levels.csv (8 rows)

The INNOVA absolute NRS means and SDs per arm at baseline, 3, 6 and 12 months (Table 3), kept
separately because `trials_test.csv` has no endpoint for an absolute level. Baseline `n` is
the intention-to-treat 47 and 49 used in Table 4, not the 48 and 50 randomised. Per-visit n
is `NA`: the paper states only that numbers varied across visits because of dropout.

## test/cytokines_parkitny.csv (6 rows)

Parkitny & Younger 2017: a 10-week, single-blind pilot in 8 women (2-week baseline, then
8 weeks of 4.5 mg naltrexone; no placebo capsules were actually given). `day = 56` marks the
drug phase, the final two weeks of dosing. The paper reports baseline and drug-phase **means
(SE) in pg/mL**, not percent changes, so `pct_change` was computed from those two reported
means and each row's `location` states the pair used. Because the means are rounded to one
decimal, small percentages (TNF-α, −5.8 %) carry noticeable rounding error. `se` of the
percent change is not derivable and is `NA`. `p_value` is the linear-mixed-model p from
Table 2; the `p_lt` flag is `yes` only for IL-6, which the paper reports as p < 0.001, so
0.001 there is a bound rather than an estimate.

## test/qst_final.csv (6 rows)

The FINAL exploratory QST paper (Due Bruun 2025, CNS Drugs), Table 2 between-group
differences (LDN minus placebo) in the complete-case population (45 LDN vs 47 placebo):
pressure pain tolerance left leg −2.6 kPa (−5.8 to 0.5) and right leg −0.5 (−3.3 to 2.2),
temporal summation −0.1 VAS cm (−0.9 to 0.7), conditioned pain modulation +2.0 kPa
(0.4 to 3.7), plus the two muscular-endurance outcomes. The `favours` column gives the sign
convention per measure: a positive difference favours LDN for every measure except temporal
summation of pain, where less pain facilitation is better and a negative difference favours LDN.

## module_screen.csv (11 rows)

The spec §3.2 decision table with a real citation in every row. It is **not** a dataset of
measured values but a record of design decisions, so it has no `source`/`location`/`verified`
triple and is therefore excluded by name from the provenance test in
`tests/testthat/test-data.R`. Traceability is instead carried by the `citations` column,
which holds "Author Year (DOI ...)" for every row and was checked programmatically to contain
a DOI on all 11 rows.

---

## Not found

These items were searched for and are deliberately absent rather than estimated.

1. **Per-dose Cmax, Tmax and AUC from Meyer 1984** (50/100/200 mg). The paper is not in PMC
   or any open repository; only its abstract was obtainable, which gives half-lives,
   clearances and urinary recovery but no per-dose concentration statistics.
2. **AUC for naltrexone or 6β-naltrexol at any dose.** Rabiner 2011 Table 1 gives only
   min–max ranges of AUC(0,t) and Cmax across its adaptive dose range, not means per dose.
3. **The oral dose used in Wall 1981.** The abstract states none, and the full text was not
   obtainable, so the 8.9 h oral terminal half-life is kept as context above rather than
   being attached to an assumed 50 mg.
4. **A formal pharmacokinetic analysis of low-dose (≤ 6 mg) naltrexone.** Bested et al. 2023
   (Pain Reports) sampled plasma at 4.5 mg; the day-1 15-60 min profile (digitised from Figure 4) and
   the day-21 1–3 h medians are now curated as `C0.25`-`C1` and `Css13h` rows in
   `training/pk_summary.csv` (R4, R4 fix round 1). The paper states that
   "detailed pharmacokinetic analyses were not performed because of the short sampling period",
   so no NCA-derived AUC, Cmax (peak) or half-life values were available from Bested.
5. **Plasma unbound fraction of 6β-naltrexol.** The label and the PK papers report protein
   binding for naltrexone only.
6. **Binding affinity of naltrexone at OGFr.** OGFr was cloned and characterised with
   radiolabelled opioid growth factor ([Met5]-enkephalin); no Ki or Kd for naltrexone itself
   at OGFr was found. This makes hypothesis H3 non-identifiable from binding data.
7. **A TLR4/MD2 Kd or a binding IC50 for naltrexone.** Only cell-based functional IC50
   values were found (Wang 2016, Selfridge 2015). Hutchinson 2008 used a single 10 µM
   concentration in HEK-TLR4 reporter cells and reports no concentration–response constant.
8. **Human plasma half-life of IL-1β, IL-17A and TGF-β.** Phase I trials of recombinant
   human IL-1β report clinical and haematological effects but no plasma half-life in their
   abstracts, and no human pharmacokinetic study of IL-17A or of TGF-β was found. Those
   three cytokines therefore appear in `test/cytokines_parkitny.csv` but have no row in
   `params/cytokine_params.csv`.
9. **A pooled baseline pain value in Häuser 2011.** The published text states only that
   baseline pain did not differ between the two disease populations; hence the Byon 2010
   baseline was used for the placebo-prior derivation.
10. **Pain-specific percent change for Younger 2009.** Only F statistics and p values are
    reported for daily pain; the percentages available are for overall symptom severity.
11. **Time after the last dose for the Weerts 2008 scan**, and the exact post-dose times of
    the individual Rabiner 2011 scans. Both are recorded as assumptions in the `time_basis`
    column of `test/pet_occupancy.csv`.
12. **A within-subject correlation or SE of change for INNOVA**, without which no standard
    error can be attached to the derived `change_nrs` rows.
