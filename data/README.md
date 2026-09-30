# Curated data for the LDN–fibromyalgia QSP model

Every CSV carries `source`, `location` and `verified` columns and is validated by
`tests/testthat/test-data.R` through `read_dataset()` in `R/data_io.R`.

`verified` values: `yes` = the number was read in the primary source; `secondary` = taken
from a review or digitized from a figure; `no` = uncertain. Missing values are `NA`.
No field contains the `#` character (it would truncate a row, because `read_dataset()`
calls `read.csv(comment.char = "#")`), and every row has a non-empty `location`.

Arm naming is `placebo` for placebo arms and `ldn_<dose>mg` for active arms. `schedule`
strings use `startday:dose;startday:dose` with strictly increasing start days beginning at 0.

---

## training/pk_summary.csv (10 rows)

Oral naltrexone (`NTX`) and 6β-naltrexol (`BN`) plasma PK summary statistics from
Verebey 1976 at 100 mg: Cmax, Tmax and terminal half-life after an acute dose and after
chronic daily dosing, n = 4. All rows are `verified = yes`, read in the paper's abstract,
and every row carries a numeric `dose_mg`.

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
*Conflict:* Due Bruun 2023 (Pain Reports crossover trial) quotes the same two statistics as
"28.8 (12.5) %" and "18.0 (14.6) %" with the bracketed values called SDs. Those SDs are not
reproducible from Younger's reported CI and n, so the primary CI was used.

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
* `pbo_emax_logit` = 1.52 (SE 0.125), the placebo Emax on the logit pain scale, Byon 2010
  Table III.
* `sd_res_nrs` = 0.48 NRS points, the within-subject measurement SD of the **0–10** NRS.
  Source: Alghadir 2018, a 24-hour test–retest study in 121 patients with knee osteoarthritis
  (ICC 0.95, SEM 0.48, MDC 1.33). No equivalent test–retest SEM measured in a fibromyalgia
  cohort was found, so this is an out-of-population estimate.

## test/pet_occupancy.csv (11 rows)

Held-out PET data. `time_h` is the time after the last dose, and two extra columns say how
far each row can be trusted:

* `time_basis` — `reported` when the paper states the sampling time, otherwise the assumption
  made and its justification.
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
4. **Any published PK study of low-dose (≤ 6 mg) naltrexone.** The Danish crossover trial
   (Due Bruun 2023, Pain Reports) sampled plasma but states that detailed pharmacokinetic
   analyses were not performed because of the short sampling period.
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
