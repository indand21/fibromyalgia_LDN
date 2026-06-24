# Main Simulation Script: QSP Model of LDN in Fibromyalgia
# Requires: mrgsolve, dplyr, ggplot2

library(mrgsolve)
library(dplyr)
library(ggplot2)

# ============================================================
# 1. LOAD MODEL
# ============================================================

mod <- mread("qsp_model", file.path("..", "model"))

# ============================================================
# 2. DEFINE DOSING REGIMEN (LDN 4.5 mg QD x 28 days)
# ============================================================

dose_events <- ev(
  amt  = 4.5,       # mg
  cmt  = 1,         # A_gut compartment
  ii   = 24,        # dosing interval (hours)
  addl = 27,        # additional doses (28 total)
  rate = 0          # bolus (oral)
)

# Alternative dose regimens for dose-response
dose_levels <- c(1.5, 3.0, 4.5, 9.0, 25.0, 50.0)

# ============================================================
# 3. SIMULATION TIME
# ============================================================

sim_time <- seq(0, 28 * 24, by = 0.5)  # 28 days in hours

# ============================================================
# 4. RUN BASELINE FIBROMYALGIA (No Drug)
# ============================================================

cat("Running baseline fibromyalgia simulation...\n")
sim_baseline <- mod %>%
  ev(dose_events %>% mutate(amt = 0)) %>%  # zero dose
  mrgsim(end = max(sim_time), delta = 0.5)

baseline_df <- as.data.frame(sim_baseline) %>%
  mutate(time_days = time / 24, Scenario = "Baseline FMQ")

# ============================================================
# 5. RUN LDN 4.5 mg QD
# ============================================================

cat("Running LDN 4.5mg simulation...\n")
sim_ldn <- mod %>%
  ev(dose_events) %>%
  mrgsim(end = max(sim_time), delta = 0.5)

ldn_df <- as.data.frame(sim_ldn) %>%
  mutate(time_days = time / 24, Scenario = "LDN 4.5mg QD")

# ============================================================
# 6. DOSE-RESPONSE SIMULATION
# ============================================================

cat("Running dose-response simulations...\n")
dose_response_results <- list()

for (dose in dose_levels) {
  cat(sprintf("  Dose: %.1f mg\n", dose))

  dose_ev <- ev(
    amt  = dose,
    cmt  = 1,
    ii   = 24,
    addl = 27,
    rate = 0
  )

  sim <- mod %>%
    ev(dose_ev) %>%
    mrgsim(end = max(sim_time), delta = 0.5)

  df <- as.data.frame(sim) %>%
    mutate(time_days = time / 24, Dose = dose)

  dose_response_results[[as.character(dose)]] <- df
}

dose_response_df <- bind_rows(dose_response_results)

# ============================================================
# 7. EXTRACT STEADY-STATE RESULTS (Day 28)
# ============================================================

get_ss_results <- function(df, label) {
  df %>%
    filter(time_days >= 27.5) %>%
    summarise(
      Scenario     = label,
      Pain_VAS_ss  = mean(Pain_VAS),
      TNF_ss       = mean(TNF),
      IL1B_ss      = mean(IL1B),
      IL6_ss       = mean(IL6),
      IL10_ss      = mean(IL10),
      BDNF_ss      = mean(BDNF),
      NGF_ss       = mean(NGF),
      TAC1_ss      = mean(TAC1),
      C_CNS_ss     = mean(C_CNS),
      TLR4_occ_ss  = mean(TLR4_occupancy),
      Pain_red_ss  = mean(Pain_reduction),
      Inf_index_ss = mean(Inflammatory_index),
      Neuro_index_ss = mean(Neurotrophin_index)
    )
}

ss_baseline <- get_ss_results(baseline_df, "Baseline FMQ")
ss_ldn     <- get_ss_results(ldn_df, "LDN 4.5mg")

ss_summary <- bind_rows(ss_baseline, ss_ldn)
print(ss_summary)

# Dose-response at steady state
ss_dose_response <- dose_response_df %>%
  filter(time_days >= 27.5) %>%
  group_by(Dose) %>%
  summarise(
    Pain_VAS_ss  = mean(Pain_VAS),
    Pain_red_ss  = mean(Pain_reduction),
    TNF_ss       = mean(TNF),
    IL1B_ss      = mean(IL1B),
    IL6_ss       = mean(IL6),
    IL10_ss      = mean(IL10),
    TLR4_occ_ss  = mean(TLR4_occupancy),
    C_CNS_ss     = mean(C_CNS)
  )

print(ss_dose_response)

# ============================================================
# 8. VISUALIZATION
# ============================================================

# Create output directory
dir.create(file.path("..", "output", "figures"), showWarnings = FALSE, recursive = TRUE)

# --- Plot 1: Pain Score Over Time ---
p1 <- bind_rows(baseline_df, ldn_df) %>%
  ggplot(aes(x = time_days, y = Pain_VAS, color = Scenario)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Pain Score (VAS) Over Time",
    subtitle = "Fibromyalgia Baseline vs LDN 4.5mg QD",
    x = "Time (days)",
    y = "Pain VAS (0-10)",
    color = "Scenario"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Baseline FMQ" = "red", "LDN 4.5mg QD" = "blue")) +
  ylim(0, 10)

ggsave(file.path("..", "output", "figures", "pain_over_time.png"), p1, width = 10, height = 6)

# --- Plot 2: Cytokine Panel ---
cytokine_data <- bind_rows(baseline_df, ldn_df) %>%
  select(time_days, Scenario, TNF, IL1B, IL6, IL10, CXCL8) %>%
  tidyr::pivot_longer(cols = c(TNF, IL1B, IL6, IL10, CXCL8),
                      names_to = "Cytokine", values_to = "Concentration")

p2 <- cytokine_data %>%
  ggplot(aes(x = time_days, y = Concentration, color = Scenario)) +
  geom_line(linewidth = 1) +
  facet_wrap(~Cytokine, scales = "free_y", ncol = 2) +
  labs(
    title = "Cytokine Dynamics: Baseline vs LDN 4.5mg",
    x = "Time (days)",
    y = "Concentration (uM)"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Baseline FMQ" = "red", "LDN 4.5mg QD" = "blue"))

ggsave(file.path("..", "output", "figures", "cytokine_panel.png"), p2, width = 12, height = 10)

# --- Plot 3: Neurotrophin Panel ---
neuro_data <- bind_rows(baseline_df, ldn_df) %>%
  select(time_days, Scenario, BDNF, NGF, TAC1, TRPV1_act) %>%
  tidyr::pivot_longer(cols = c(BDNF, NGF, TAC1, TRPV1_act),
                      names_to = "Target", values_to = "Concentration")

p3 <- neuro_data %>%
  ggplot(aes(x = time_days, y = Concentration, color = Scenario)) +
  geom_line(linewidth = 1) +
  facet_wrap(~Target, scales = "free_y", ncol = 2) +
  labs(
    title = "Pain Mediator Dynamics: Baseline vs LDN 4.5mg",
    x = "Time (days)",
    y = "Concentration / Activity (uM)"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Baseline FMQ" = "red", "LDN 4.5mg QD" = "blue"))

ggsave(file.path("..", "output", "figures", "neurotrophin_panel.png"), p3, width = 12, height = 10)

# --- Plot 4: TLR4 Occupancy & CNS Concentration ---
p4 <- bind_rows(baseline_df, ldn_df) %>%
  ggplot(aes(x = time_days)) +
  geom_line(aes(y = TLR4_occupancy, color = "TLR4 Occupancy (%)"), linewidth = 1.2) +
  geom_line(aes(y = C_CNS * 1000, color = "CNS Conc (nM x1000)"), linewidth = 1.2) +
  facet_wrap(~Scenario) +
  labs(
    title = "Naltrexone PK/PD: TLR4 Occupancy",
    x = "Time (days)",
    y = "Value",
    color = "Metric"
  ) +
  theme_minimal()

ggsave(file.path("..", "output", "figures", "tlr4_occupancy.png"), p4, width = 12, height = 6)

# --- Plot 5: Dose-Response Curve ---
p5 <- ss_dose_response %>%
  ggplot(aes(x = Dose, y = Pain_red_ss)) +
  geom_line(linewidth = 1.2, color = "darkblue") +
  geom_point(size = 3, color = "darkblue") +
  scale_x_log10() +
  labs(
    title = "Dose-Response: Pain Reduction vs Naltrexone Dose",
    subtitle = "Steady-state (Day 28)",
    x = "Naltrexone Dose (mg, log scale)",
    y = "Pain Reduction (%)"
  ) +
  theme_minimal() +
  geom_hline(yintercept = 30, linetype = "dashed", color = "red", alpha = 0.5) +
  annotate("text", x = 2, y = 32, label = "MCID (30%)", color = "red")

ggsave(file.path("..", "output", "figures", "dose_response.png"), p5, width = 8, height = 6)

# --- Plot 6: Inflammatory vs Neurotrophin Index ---
p6 <- bind_rows(baseline_df, ldn_df) %>%
  ggplot(aes(x = time_days)) +
  geom_line(aes(y = Inflammatory_index, color = "Inflammatory Index"), linewidth = 1.2) +
  geom_line(aes(y = Neurotrophin_index, color = "Neurotrophin Index"), linewidth = 1.2) +
  facet_wrap(~Scenario) +
  labs(
    title = "Composite Indices Over Time",
    x = "Time (days)",
    y = "Index (normalized to baseline)",
    color = "Index"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("Inflammatory Index" = "firebrick", "Neurotrophin Index" = "steelblue"))

ggsave(file.path("..", "output", "figures", "composite_indices.png"), p6, width = 12, height = 6)

# ============================================================
# 9. SAVE RESULTS
# ============================================================

dir.create(file.path("..", "output", "results"), showWarnings = FALSE, recursive = TRUE)

write.csv(ss_summary, file.path("..", "output", "results", "steady_state_summary.csv"), row.names = FALSE)
write.csv(ss_dose_response, file.path("..", "output", "results", "dose_response_ss.csv"), row.names = FALSE)
write.csv(baseline_df, file.path("..", "output", "results", "baseline_simulation.csv"), row.names = FALSE)
write.csv(ldn_df, file.path("..", "output", "results", "ldn_simulation.csv"), row.names = FALSE)
write.csv(dose_response_df, file.path("..", "output", "results", "dose_response_full.csv"), row.names = FALSE)

cat("\n=== Simulation Complete ===\n")
cat("Results saved to: output/results/\n")
cat("Figures saved to: output/figures/\n")
cat("\n--- Steady-State Summary ---\n")
print(ss_summary)
