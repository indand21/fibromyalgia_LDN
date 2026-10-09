# Archived analyses cited in the manuscript

The primary analysis is in `output/` and is reproduced by `Rscript scripts/run_all.R`. This folder holds two secondary analyses that the manuscript quotes.

- `clean_run_2026-10-08_preR4_verified/`: the full pipeline with the PK fitted to the high-dose (Verebey 1976, 100 mg) data only, with linear first pass. This is the "high-dose-only PK" reference: occupancy at 4.5 mg is 0.762. It was verified independently.
- `linear_bested_pk_sensitivity_2026-10-08/`: the PK refitted with the Bested 2023 low-dose data but with the linear (dose-independent) structure forced, followed by scripts 02-06 and 08. This is the lower end of the occupancy bracket: 0.739 at 4.5 mg. The one-line code change is in `01_fit_pk_linear_forced.R`.

Each folder has an `ARCHIVE_NOTE.txt`. Files omitted to keep the repository small:
- the R binary caches (`*.rds`), which the code regenerates;
- the 635 MB per-subject `virtual_trials.csv`, which the manuscript does not cite.
