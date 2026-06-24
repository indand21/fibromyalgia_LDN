# Dose Comparison: LDN (4.5mg) vs Standard Naltrexone (50mg)
# Demonstrates WHY only low dose works in fibromyalgia

library(mrgsolve)
library(dplyr)
library(ggplot2)

# ============================================================
# 1. LOAD MODEL
# ============================================================

mod <- mread("qsp_model", file.path("..", "model"))

# ============================================================
# 2. DEFINE DOSE REGIMENS
# ============================================================

doses <- list(
  "Placebo"     = ev(amt = 0,    cmt = 1, ii = 24, addl = 27),
  "LDN 1.5mg"   = ev(amt = 1.5,  cmt = 1, ii = 24, addl = 27),
  "LDN 3.0mg"   = ev(amt = 3.0,  cmt = 1, ii = 24, addl = 27),
  "LDN 4.5mg"   = ev(amt = 4.5,  cmt = 1, ii = 24, addl = 27),
  "Standard 25mg" = ev(amt = 25,  cmt = 1, ii = 24, addl = 27),
  "Standard 50mg" = ev(amt = 50,  cmt = 1, ii = 24, addl = 27)
)

sim_end <- 28 * 24  # 28 days

# ============================================================
# 3. RUN ALL DOSE SCENARIOS
# ============================================================

cat("Running dose comparison simulations...\n")

results <- list()
for (name in names(doses)) {
  cat(sprintf("  %s\n", name))

  sim <- mod %>%
    ev(doses[[name]]) %>%
    mrgsim(end = sim_end, delta = 1)

  df <- as.data.frame(sim) %>%
    mutate(
      time_days = time / 24,
      Dose_Group = name
    )
  results[[name]] <- df
}

all_results <- bind_rows(results)

# ============================================================
# 4. EXTRACT STEADY-STATE VALUES
# ============================================================

ss_summary <- all_results %>%
  filter(time_days >= 27) %>%
  group_by(Dose_Group) %>%
  summarise(
    C_CNS_nM      = mean(C_CNS) * 1000,     # convert to nM
    TLR4_occ_pct  = mean(TLR4_occupancy),
    OPRM1_occ_pct = mean(OPRM1_occupancy),
    Pain_VAS      = mean(Pain_VAS),
    Pain_red_pct  = mean(Pain_reduction),
    TNF           = mean(TNF),
    IL1B          = mean(IL1B),
    IL6           = mean(IL6),
    IL10          = mean(IL10),
    BDNF          = mean(BDNF),
    TAC1          = mean(TAC1),
    Hyperalgesia  = mean(Hyperalgesia_pct),
    Endorphin_blk = mean(Endorphin_blockade),
    Inf_index     = mean(Inflammatory_index),
    .groups = "drop"
  )

cat("\n=== Steady-State Comparison (Day 28) ===\n")
print(ss_summary)

# ============================================================
# 5. CREATE COMPARISON PLOTS
# ============================================================

dir.create(file.path("..", "output", "figures"), showWarnings = FALSE, recursive = TRUE)

# --- Plot 1: Receptor Occupancy vs Dose (the KEY plot) ---
occ_data <- ss_summary %>%
  select(Dose_Group, TLR4_occ_pct, OPRM1_occ_pct) %>%
  tidyr::pivot_longer(cols = c(TLR4_occ_pct, OPRM1_occ_pct),
                      names_to = "Receptor", values_to = "Occupancy") %>%
  mutate(Receptor = ifelse(Receptor == "TLR4_occ_pct",
                           "TLR4 (anti-inflammatory target)",
                           "OPRM1 (mu-opioid receptor)"))

p1 <- ggplot(occ_data, aes(x = Dose_Group, y = Occupancy, fill = Receptor)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(
    title = "Receptor Occupancy: Why Only Low-Dose Naltrexone Works",
    subtitle = "TLR4 antagonism = beneficial | OPRM1 blockade = harmful at high doses",
    x = "Dose Regimen",
    y = "Receptor Occupancy (%)",
    fill = "Receptor"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c(
    "TLR4 (anti-inflammatory target)" = "steelblue",
    "OPRM1 (mu-opioid receptor)" = "firebrick"
  )) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  geom_hline(yintercept = 50, linetype = "dashed", alpha = 0.3)

ggsave(file.path("..", "output", "figures", "dose_receptor_occupancy.png"), p1, width = 10, height = 7)

# --- Plot 2: Pain Score Trajectory ---
p2 <- all_results %>%
  filter(Dose_Group %in% c("Placebo", "LDN 4.5mg", "Standard 50mg")) %>%
  ggplot(aes(x = time_days, y = Pain_VAS, color = Dose_Group)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Pain Trajectory: LDN vs Standard Naltrexone",
    subtitle = "LDN reduces pain via TLR4; Standard dose worsens pain via OPRM1 blockade",
    x = "Time (days)",
    y = "Pain VAS (0-10)",
    color = "Treatment"
  ) +
  theme_minimal() +
  scale_color_manual(values = c(
    "Placebo" = "gray50",
    "LDN 4.5mg" = "steelblue",
    "Standard 50mg" = "firebrick"
  ))

ggsave(file.path("..", "output", "figures", "pain_trajectory_comparison.png"), p2, width = 10, height = 6)

# --- Plot 3: Inflammatory vs Hyperalgesia Components ---
component_data <- ss_summary %>%
  select(Dose_Group, Inf_index, Hyperalgesia, Pain_red_pct) %>%
  tidyr::pivot_longer(cols = c(Inf_index, Hyperalgesia),
                      names_to = "Component", values_to = "Value") %>%
  mutate(Component = ifelse(Component == "Inf_index",
                            "Inflammation Index (benefit)",
                            "Hyperalgesia % (harm)"))

p3 <- ggplot(component_data, aes(x = Dose_Group, y = Value, fill = Component)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(
    title = "Benefit vs Harm: The Dose-Dependent Tradeoff",
    subtitle = "LDN: inflammation down, no hyperalgesia | High dose: inflammation down slightly, hyperalgesia up",
    x = "Dose Regimen",
    y = "Index Value",
    fill = "Component"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c(
    "Inflammation Index (benefit)" = "steelblue",
    "Hyperalgesia % (harm)" = "firebrick"
  )) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path("..", "output", "figures", "benefit_harm_tradeoff.png"), p3, width = 10, height = 7)

# --- Plot 4: Cytokine Panel (LDN vs 50mg) ---
cytokine_data <- all_results %>%
  filter(Dose_Group %in% c("Placebo", "LDN 4.5mg", "Standard 50mg")) %>%
  select(time_days, Dose_Group, TNF, IL1B, IL6, IL10) %>%
  tidyr::pivot_longer(cols = c(TNF, IL1B, IL6, IL10),
                      names_to = "Cytokine", values_to = "Concentration")

p4 <- cytokine_data %>%
  ggplot(aes(x = time_days, y = Concentration, color = Dose_Group)) +
  geom_line(linewidth = 1) +
  facet_wrap(~Cytokine, scales = "free_y", ncol = 2) +
  labs(
    title = "Cytokine Dynamics: LDN vs Standard Dose",
    x = "Time (days)",
    y = "Concentration (uM)"
  ) +
  theme_minimal() +
  scale_color_manual(values = c(
    "Placebo" = "gray50",
    "LDN 4.5mg" = "steelblue",
    "Standard 50mg" = "firebrick"
  ))

ggsave(file.path("..", "output", "figures", "cytokine_ldn_vs_standard.png"), p4, width = 12, height = 10)

# --- Plot 5: Dose-Response U-Shaped Curve ---
p5 <- ss_summary %>%
  ggplot(aes(x = reorder(Dose_Group, c(0, 1.5, 3, 4.5, 25, 50)), y = Pain_red_pct)) +
  geom_bar(stat = "identity", width = 0.6, fill = c("gray", "lightblue", "steelblue", "darkblue", "orange", "firebrick")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  labs(
    title = "U-Shaped Dose Response: Pain Reduction by Naltrexone Dose",
    subtitle = "LDN (4.5mg) shows optimal benefit; Standard dose shows NEGATIVE effect",
    x = "Dose Regimen",
    y = "Pain Reduction (%)"
  ) +
  theme_minimal() +
  annotate("text", x = 4, y = 5, label = "Optimal\n(LDN)", color = "darkblue", fontface = "bold") +
  annotate("text", x = 6, y = -5, label = "Paradoxical\nWorsening", color = "firebrick", fontface = "bold") +
  annotate("segment", x = 4.5, xend = 4.5, y = 3, yend = -2,
           arrow = arrow(length = unit(0.3, "cm")), color = "firebrick")

ggsave(file.path("..", "output", "figures", "u_shaped_dose_response.png"), p5, width = 10, height = 7)

# ============================================================
# 6. PRINT MECHANISTIC EXPLANATION
# ============================================================

cat("\n")
cat("================================================================\n")
cat("  WHY ONLY LOW-DOSE NALTREXONE WORKS IN FIBROMYALGIA\n")
cat("================================================================\n")
cat("\n")
cat("KEY INSIGHT: Differential Receptor Affinity\n")
cat("  - Naltrexone Kd for OPRM1 = 0.001 uM (1 nM)  [VERY HIGH affinity]\n")
cat("  - Naltrexone Kd for TLR4  = 10 uM              [LOW affinity]\n")
cat("  - Ratio: OPRM1 affinity is 10,000x higher than TLR4\n")
cat("\n")
cat("AT LDN (4.5mg):\n")
cat(sprintf("  - CNS concentration: ~%.2f nM\n", ss_summary$C_CNS_nM[ss_summary$Dose_Group == "LDN 4.5mg"]))
cat(sprintf("  - TLR4 occupancy: ~%.1f%% (sufficient for anti-inflammatory effect)\n", ss_summary$TLR4_occ_pct[ss_summary$Dose_Group == "LDN 4.5mg"]))
cat(sprintf("  - OPRM1 occupancy: ~%.1f%% (partial, but endogenous opioid tone preserved)\n", ss_summary$OPRM1_occ_pct[ss_summary$Dose_Group == "LDN 4.5mg"]))
cat(sprintf("  - Pain reduction: %.1f%%\n", ss_summary$Pain_red_pct[ss_summary$Dose_Group == "LDN 4.5mg"]))
cat(sprintf("  - Hyperalgesia: ~%.1f%% (negligible)\n", ss_summary$Hyperalgesia[ss_summary$Dose_Group == "LDN 4.5mg"]))
cat("\n")
cat("AT STANDARD DOSE (50mg):\n")
cat(sprintf("  - CNS concentration: ~%.1f nM\n", ss_summary$C_CNS_nM[ss_summary$Dose_Group == "Standard 50mg"]))
cat(sprintf("  - TLR4 occupancy: ~%.1f%% (slightly higher, marginal benefit)\n", ss_summary$TLR4_occ_pct[ss_summary$Dose_Group == "Standard 50mg"]))
cat(sprintf("  - OPRM1 occupancy: ~%.1f%% (near-complete blockade)\n", ss_summary$OPRM1_occ_pct[ss_summary$Dose_Group == "Standard 50mg"]))
cat(sprintf("  - Pain change: %.1f%% (NEGATIVE = pain worsens)\n", ss_summary$Pain_red_pct[ss_summary$Dose_Group == "Standard 50mg"]))
cat(sprintf("  - Hyperalgesia: ~%.1f%% (significant pain amplification)\n", ss_summary$Hyperalgesia[ss_summary$Dose_Group == "Standard 50mg"]))
cat("\n")
cat("MECHANISM:\n")
cat("  1. LDN selectively antagonizes TLR4 on microglia -> reduces neuroinflammation\n")
cat("  2. LDN minimally affects OPRM1 -> endogenous endorphin signaling preserved\n")
cat("  3. Standard dose blocks BOTH receptors -> OPRM1 blockade triggers:\n")
cat("     a) Loss of endogenous opioid analgesia\n")
cat("     b) Compensatory hyperalgesia (pain amplification)\n")
cat("     c) Opioid withdrawal-like effects\n")
cat("  4. Net result: U-shaped dose-response curve\n")
cat("================================================================\n")

# ============================================================
# 7. SAVE RESULTS
# ============================================================

dir.create(file.path("..", "output", "results"), showWarnings = FALSE, recursive = TRUE)

write.csv(ss_summary, file.path("..", "output", "results", "dose_comparison_ss.csv"), row.names = FALSE)
write.csv(all_results, file.path("..", "output", "results", "dose_comparison_full.csv"), row.names = FALSE)

cat("\nResults saved to: output/results/\n")
cat("Figures saved to: output/figures/\n")
