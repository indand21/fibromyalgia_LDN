# Monte Carlo Uncertainty + Virtual Patient Population
# QSP Model v2: Naltrexone in Fibromyalgia

source("qsp_model_functions.R")
library(ggplot2)

cat("=== Monte Carlo Uncertainty Analysis ===\n")

# ============================================================
# 1. PARAMETER UNCERTAINTY
# ============================================================

# Define uncertainty ranges (CV = coefficient of variation)
uncertain_params <- data.frame(
  param = c("kon","koff","kpTNF","kdTNF","kpIL6","kdIL6","kpIL17","kdIL17",
            "k_up","k_reb","Kd_OPRM1","k_hyper","kIL6_BDNF","kB_TAC1",
            "kpBDNF","kdBDNF","Em_BDNF","Em_TNF","Em_TAC1","kpCRH","kdCRH",
            "kpNPY","kdNPY","ksHT","kdHT"),
  nominal = c(0.1, 1.0, 0.5, 0.693, 0.4, 0.116, 0.08, 0.116,
              0.05, 0.02, 0.001, 0.5, 0.5, 1.0,
              0.1, 0.023, 0.5, 0.3, 0.4, 0.1, 0.693,
              0.05, 0.116, 0.1, 0.5),
  CV = c(0.5, 0.5, 0.3, 0.3, 0.3, 0.3, 0.5, 0.3,
         0.5, 0.5, 0.5, 0.5, 0.3, 0.3,
         0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
         0.3, 0.3, 0.3, 0.3)
)

n_mc <- 200
cat(sprintf("Running %d Monte Carlo samples...\n", n_mc))

set.seed(42)
mc_results <- list()

for (i in 1:n_mc) {
  if(i%%50==0) cat(sprintf("  %d/%d\n", i, n_mc))

  # Sample parameters from log-normal distribution
  p_mod <- get_params()
  for (j in 1:nrow(uncertain_params)) {
    nm <- uncertain_params$param[j]
    nom <- uncertain_params$nominal[j]
    cv <- uncertain_params$CV[j]
    sigma <- sqrt(log(1 + cv^2))
    mu <- log(nom) - sigma^2/2
    p_mod[[nm]] <- exp(rnorm(1, mu, sigma))
  }

  # Run for LDN 4.5mg (Ccns=0.01)
  df <- run_qsp(0.01, p=p_mod, t_end=28*24)
  ss <- get_ss(df)
  mc_results[[i]] <- data.frame(
    sample=i, Pain_VAS=ss$Pain, Pain_red=ss$Pain_red,
    TNF=ss$TNF, IL6=ss$IL6, IL10=ss$IL10, IL17=ss$IL17,
    BDNF=ss$BDNF, TAC1=ss$TAC1, CPT=ss$CPT,
    OPRM1_fold=ss$OPRM1_fold, Endo_fold=ss$Endo_fold
  )
}

mc_df <- bind_rows(mc_results)

# Statistics
cat("\n=== Monte Carlo Results (LDN 4.5mg, N=200) ===\n")
mc_stats <- mc_df %>%
  summarise(across(Pain_VAS:Endo_fold, list(
    mean=~mean(.), sd=~sd(.),
    q05=~quantile(., 0.05), q50=~quantile(., 0.5), q95=~quantile(., 0.95)
  ))) %>% tidyr::pivot_longer(everything(), names_to=c("Variable","Stat"), names_sep="_") %>%
  tidyr::pivot_wider(names_from=Stat, values_from=value)

print(as.data.frame(mc_stats))

# Plots
p1 <- ggplot(mc_df, aes(x=Pain_red)) +
  geom_histogram(bins=30, fill="steelblue", alpha=0.7) +
  geom_vline(xintercept=mean(mc_df$Pain_red), color="red", linewidth=1.2) +
  geom_vline(xintercept=quantile(mc_df$Pain_red, 0.05), color="red", linetype="dashed") +
  geom_vline(xintercept=quantile(mc_df$Pain_red, 0.95), color="red", linetype="dashed") +
  labs(title="Monte Carlo: Pain Reduction Distribution (LDN 4.5mg)",
       subtitle=sprintf("N=%d, Mean=%.1f%%, 90%% CI: [%.1f%%, %.1f%%]",
                        n_mc, mean(mc_df$Pain_red),
                        quantile(mc_df$Pain_red,0.05), quantile(mc_df$Pain_red,0.95)),
       x="Pain Reduction (%)", y="Count") +
  theme_minimal()
ggsave("../output/figures/mc_pain_distribution.png", p1, width=10, height=6)

p2 <- mc_df %>%
  select(TNF, IL6, IL10, BDNF) %>%
  tidyr::pivot_longer(everything(), names_to="Var", values_to="Val") %>%
  ggplot(aes(x=Val)) +
  geom_histogram(bins=25, fill="steelblue", alpha=0.7) +
  facet_wrap(~Var, scales="free", ncol=2) +
  labs(title="Monte Carlo: Biomarker Distributions", x="Value", y="Count") +
  theme_minimal()
ggsave("../output/figures/mc_biomarkers.png", p2, width=12, height=8)

write.csv(mc_df, "../output/results/monte_carlo_results.csv", row.names=FALSE)

# ============================================================
# 2. VIRTUAL PATIENT POPULATION
# ============================================================

cat("\n=== Virtual Patient Population ===\n")

# Create diverse virtual patients varying baseline disease state
n_patients <- 100
cat(sprintf("Generating %d virtual patients...\n", n_patients))

set.seed(123)
patient_results <- list()

for (i in 1:n_patients) {
  if(i%%25==0) cat(sprintf("  Patient %d/%d\n", i, n_patients))

  # Vary baseline disease severity
  disease_severity <- runif(1, 0.5, 1.5)  # multiplier

  y0_patient <- get_y0(
    TNF = 0.08 * disease_severity * runif(1, 0.7, 1.3),
    IL1B = 0.05 * disease_severity * runif(1, 0.7, 1.3),
    IL6 = 0.10 * disease_severity * runif(1, 0.7, 1.3),
    IL10 = 0.015 * runif(1, 0.5, 1.5),
    IL17 = 0.03 * disease_severity * runif(1, 0.7, 1.3),
    BDNF = 0.15 * disease_severity * runif(1, 0.7, 1.3),
    NGF = 0.10 * disease_severity * runif(1, 0.7, 1.3),
    TAC1 = 0.08 * disease_severity * runif(1, 0.7, 1.3),
    Pain = min(10, 7.0 * disease_severity * runif(1, 0.8, 1.2)),
    CPT = max(5, 25.0 / disease_severity * runif(1, 0.7, 1.3)),
    OPRM1d = runif(1, 0.7, 1.3),
    Endo = 0.002 * runif(1, 0.5, 1.5)
  )

  # Vary OPRM1 sensitivity
  p_patient <- get_params(
    Kd_OPRM1 = 0.001 * runif(1, 0.5, 2.0),
    k_up = 0.05 * runif(1, 0.5, 2.0),
    k_reb = 0.02 * runif(1, 0.5, 2.0)
  )

  # Run LDN 4.5mg
  df <- run_qsp(0.01, p=p_patient, y0=y0_patient, t_end=28*24)
  ss <- get_ss(df)

  # Also run placebo
  df_p <- run_qsp(0, p=p_patient, y0=y0_patient, t_end=28*24)
  ss_p <- get_ss(df_p)

  patient_results[[i]] <- data.frame(
    patient=i,
    severity=disease_severity,
    baseline_pain=y0_patient["Pain"],
    placebo_pain=ss_p$Pain,
    ldn_pain=ss$Pain,
    # Pain reduction from INITIAL baseline (not placebo)
    pain_red_from_baseline=(y0_patient["Pain"] - ss$Pain)/y0_patient["Pain"] * 100,
    # LDN benefit over placebo (negative = LDN worse than natural resolution)
    ldn_vs_placebo=(ss_p$Pain - ss$Pain)/ss_p$Pain * 100,
    baseline_cpt=y0_patient["CPT"],
    ldn_cpt=ss$CPT,
    cpt_improvement=(ss$CPT - y0_patient["CPT"])/y0_patient["CPT"] * 100,
    OPRM1_fold=ss$OPRM1_fold,
    Endo_fold=ss$Endo_fold,
    responder=ifelse((y0_patient["Pain"] - ss$Pain)/y0_patient["Pain"] > 0.15, "Responder", "Non-responder")
  )
}

patient_df <- bind_rows(patient_results)

cat("\n=== Virtual Patient Summary ===\n")
cat(sprintf("  Responders (>15%% reduction from baseline): %d/%d (%.0f%%)\n",
            sum(patient_df$responder=="Responder"), n_patients,
            mean(patient_df$responder=="Responder")*100))
cat(sprintf("  Mean pain reduction from baseline: %.1f%% (SD: %.1f%%)\n",
            mean(patient_df$pain_red_from_baseline), sd(patient_df$pain_red_from_baseline)))
cat(sprintf("  Mean CPT improvement: %.1f%%\n", mean(patient_df$cpt_improvement)))
cat(sprintf("  Mean OPRM1 fold change: %.2fx\n", mean(patient_df$OPRM1_fold)))
cat(sprintf("  Mean Endorphin fold change: %.2fx\n", mean(patient_df$Endo_fold)))

p3 <- ggplot(patient_df, aes(x=pain_red_from_baseline, fill=responder)) +
  geom_histogram(bins=25, alpha=0.7) +
  labs(title="Virtual Patient Population: Pain Response to LDN 4.5mg",
       subtitle=sprintf("N=%d patients, %d%% responders", n_patients,
                        round(mean(patient_df$responder=="Responder")*100)),
       x="Pain Reduction (%)", y="Count", fill="Response") +
  theme_minimal() +
  scale_fill_manual(values=c("Responder"="steelblue","Non-responder"="gray70"))
ggsave("../output/figures/virtual_patient_response.png", p3, width=10, height=6)

p4 <- ggplot(patient_df, aes(x=severity, y=pain_red_from_baseline, color=responder)) +
  geom_point(size=2, alpha=0.7) +
  geom_smooth(method="lm", se=TRUE, alpha=0.2) +
  labs(title="Disease Severity vs LDN Response",
       x="Disease Severity (multiplier)", y="Pain Reduction (%)") +
  theme_minimal() +
  scale_color_manual(values=c("Responder"="steelblue","Non-responder"="gray70"))
ggsave("../output/figures/severity_vs_response.png", p4, width=10, height=6)

write.csv(patient_df, "../output/results/virtual_patients.csv", row.names=FALSE)

cat("\nAll analyses complete.\n")
