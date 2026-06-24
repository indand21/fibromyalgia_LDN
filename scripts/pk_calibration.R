# Naltrexone PK Model Calibration
# Based on published clinical pharmacokinetic data
# References:
#   Verebey et al. 1976 - Single-dose PK in humans
#   Meyer et al. 1984 - Dose-ranging PK (50-200mg)
#   FDA ReVia label - Approved PK parameters
#   Krieter et al. 2019 - Intranasal PK (half-life)
#   Toljan & Vrooman 2018 - LDN review (CNS penetration)

library(deSolve)
library(dplyr)
library(ggplot2)

# ============================================================
# 1. CLINICAL PK DATA FROM LITERATURE
# ============================================================

# Verebey 1976: 100mg oral naltrexone in healthy volunteers
# Mean plasma naltrexone concentration-time profile
verebey_100mg <- data.frame(
  Time_hr = c(0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 10, 12, 24),
  Conc_nM = c(0, 60, 80, 70, 55, 35, 22, 12, 7, 4.5, 3, 0.8)
)

# Meyer 1984: 100mg oral dose
# Mean (SD) Cmax = 19.6 (17.9) ng/mL at Tmax = 1 hour
# 6-beta-naltrexol Cmax = 206.8 (78.1) ng/mL at Tmax = 1 hour
meyer_100mg <- data.frame(
  Time_hr = c(0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 24),
  # Converted from ng/mL to nM (MW naltrexone = 341.4 g/mol)
  Conc_nM = c(0, 40, 57.4, 50, 38, 25, 16, 9, 5.5, 2.5, 0.6)
)

# FDA label: 50mg oral dose
# Mean Cmax ~ 10-15 ng/mL at ~1 hour
# t1/2 naltrexone = 4 hours, t1/2 6-beta-naltrexol = 13 hours
fda_50mg <- data.frame(
  Time_hr = c(0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 10, 12, 24),
  Conc_nM = c(0, 22, 32, 28, 22, 14, 9, 5, 3, 2, 1.3, 0.3)
)

# LDN 4.5mg: No direct PK data published, but can be estimated
# by scaling from 50mg data (assuming dose-linear PK)
# 4.5mg/50mg = 0.09 scaling factor
ldn_4.5mg_est <- data.frame(
  Time_hr = fda_50mg$Time_hr,
  Conc_nM = fda_50mg$Conc_nM * (4.5/50)  # Linear scaling
)

cat("=== Naltrexone Clinical PK Data ===\n")
cat("Verebey 1976: 100mg, Cmax ~80 nM at 1 hr\n")
cat("Meyer 1984:   100mg, Cmax ~57 nM at 1 hr\n")
cat("FDA label:     50mg, Cmax ~32 nM at 1 hr\n")
cat("LDN 4.5mg:    Estimated Cmax ~2.9 nM at 1 hr\n")
cat("t1/2 naltrexone:      ~4 hours (oral)\n")
cat("t1/2 6-beta-naltrexol: ~13 hours (active metabolite)\n")

# ============================================================
# 2. PK MODEL: 2-COMPARTMENT WITH FIRST-ORDER ABSORPTION
# ============================================================
# Model structure:
#   Aa (absorption site) -> C (central) -> C_periph (peripheral)
#   C -> CNS compartment (for TLR4/OPRM1 binding)
#   C -> Metabolite (6-beta-naltrexol)
#   First-pass metabolism reduces oral bioavailability

pk_ode <- function(t, y, p) {
  with(as.list(c(y, p)), {
    # Absorption
    dAa <- -Ka * Aa

    # First-pass: fraction reaching systemic circulation
    F_sys <- Fmax * (1 - Emax_firstpass * C / (EC50_firstpass + C))

    # Central compartment
    dC <- (Ka * Aa * F_sys) / V1 -
          (CL / V1 + Q / V1) * C +
          (Q / V2) * C_periph -
          (CLmet / V1) * C

    # Peripheral compartment
    dC_periph <- (Q / V2) * C - (Q / V2) * C_periph

    # CNS compartment (brain penetration)
    # Kp = brain:plasma partition coefficient
    dC_CNS <- Kin * C * Kp - Kout * C_CNS

    # Active metabolite (6-beta-naltrexol)
    # Formed by CYP2D6/AKR1C4, has 1/12 to 1/50 potency of naltrexone
    dCmet <- (CLmet / V1) * C - (CLmet_elim / V1) * Cmet

    list(c(dAa, dC, dC_periph, dC_CNS, dCmet))
  })
}

# ============================================================
# 3. PK PARAMETERS (calibrated to literature)
# ============================================================

# Naltrexone molecular weight
MW <- 341.4  # g/mol

# Parameters from literature calibration
pk_params <- c(
  # Absorption
  Ka = 1.5,            # 1/hr, first-order absorption (fast)
  Tlag = 0.25,         # hr, lag time

  # First-pass metabolism (reduces oral bioavailability)
  Fmax = 0.50,         # Maximum bioavailability (50% - oral naltrexone has ~50% BA)
  Emax_firstpass = 0.5, # First-pass extraction (moderate)
  EC50_firstpass = 100,  # nM, concentration for half-maximal first-pass

  # Distribution
  V1 = 500,            # L, central compartment (liver-rich)
  V2 = 800,            # L, peripheral compartment
  Q = 150,             # L/hr, inter-compartmental clearance

  # Elimination
  CL = 87,             # L/hr, systemic clearance (t1/2 ~4 hr: 0.693*500/87=4)
  CLmet = 50,          # L/hr, clearance to metabolite (CYP2D6)
  CLmet_elim = 25,     # L/hr, metabolite elimination clearance (t1/2 ~13 hr)

  # CNS penetration
  Kp = 0.028,          # brain:plasma ratio (2.8% per literature)
  Kin = 0.5,           # 1/hr, CNS influx rate
  Kout = 0.5           # 1/hr, CNS efflux rate
)

# ============================================================
# 4. SIMULATE SINGLE DOSE (100mg) AND COMPARE TO VEREBEY 1976
# ============================================================

cat("\n=== PK Model Calibration: 100mg Single Dose ===\n")

y0_pk <- c(Aa=100*1000, C=0, C_periph=0, C_CNS=0, Cmet=0)
# Aa in ug (100mg = 100,000 ug)

t_seq <- seq(0, 24, by=0.1)

out_100 <- ode(y=y0_pk, times=t_seq, func=pk_ode, parms=pk_params, method="lsoda")
df_100 <- as.data.frame(out_100) %>%
  mutate(
    # Convert to nM: C is in ug/L, MW=341.4 g/mol
    # 1 ug/L = 1000 ng/L = 1000/341.4 nM = 2.93 nM
    C_nM = C * 1000 / MW,
    Cmet_nM = Cmet * 1000 / MW,
    C_CNS_nM = C_CNS * 1000 / MW
  )

# Compare to Verebey data
verebey_compare <- merge(
  df_100 %>% select(time, C_nM) %>%
    mutate(Time_hr = round(time, 1)) %>%
    filter(Time_hr %in% verebey_100mg$Time_hr),
  verebey_100mg,
  by = "Time_hr"
)

cat(sprintf("  Modeled Cmax: %.1f nM at %.1f hr\n",
            max(df_100$C_nM), df_100$time[which.max(df_100$C_nM)]))
cat(sprintf("  Verebey Cmax: 80 nM at 1.0 hr\n"))

# ============================================================
# 5. PARAMETER OPTIMIZATION (fit to Verebey 100mg data)
# ============================================================

cat("\n=== Optimizing PK Parameters ===\n")

# Objective function: sum of squared errors
pk_objective <- function(par_vec, data) {
  p <- pk_params
  p["Ka"] <- par_vec[1]
  p["Fmax"] <- par_vec[2]
  p["CL"] <- par_vec[3]
  p["V1"] <- par_vec[4]
  p["Kp"] <- par_vec[5]

  y0 <- c(Aa=100*1000, C=0, C_periph=0, C_CNS=0, Cmet=0)
  out <- ode(y=y0, times=seq(0, 24, by=0.1), func=pk_ode, parms=p, method="lsoda")
  df <- as.data.frame(out) %>%
    mutate(C_nM = C * 1000 / MW, Time_hr = round(time, 1))

  pred <- merge(df %>% select(Time_hr, C_nM), data, by = "Time_hr")
  sse <- sum((pred$C_nM.x - pred$Conc_nM)^2)
  return(sse)
}

# Initial guesses
par_init <- c(Ka=1.5, Fmax=0.50, CL=87, V1=500, Kp=0.028)

# Optimize
opt <- optim(par_init, pk_objective, data=verebey_100mg,
             method="L-BFGS-B",
             lower=c(0.5, 0.2, 50, 200, 0.01),
             upper=c(5, 0.8, 200, 1500, 0.1),
             control=list(maxit=1000))

cat(sprintf("  Optimized Ka:    %.3f 1/hr\n", opt$par[1]))
cat(sprintf("  Optimized Fmax:  %.4f (%.1f%%)\n", opt$par[2], opt$par[2]*100))
cat(sprintf("  Optimized CL:    %.1f L/hr\n", opt$par[3]))
cat(sprintf("  Optimized V1:    %.0f L\n", opt$par[4]))
cat(sprintf("  Optimized Kp:    %.4f\n", opt$par[5]))
cat(sprintf("  Final SSE:       %.2f\n", opt$value))

# Update parameters
pk_params_opt <- pk_params
pk_params_opt["Ka"] <- opt$par[1]
pk_params_opt["Fmax"] <- opt$par[2]
pk_params_opt["CL"] <- opt$par[3]
pk_params_opt["V1"] <- opt$par[4]
pk_params_opt["Kp"] <- opt$par[5]

# ============================================================
# 6. VALIDATE OPTIMIZED MODEL
# ============================================================

cat("\n=== Validation: Multiple Doses ===\n")

# Simulate 50mg single dose
y0_50 <- c(Aa=50*1000, C=0, C_periph=0, C_CNS=0, Cmet=0)
out_50 <- ode(y=y0_50, times=t_seq, func=pk_ode, parms=pk_params_opt, method="lsoda")
df_50 <- as.data.frame(out_50) %>% mutate(C_nM = C * 1000 / MW)

# Simulate 4.5mg single dose (LDN)
y0_4.5 <- c(Aa=4.5*1000, C=0, C_periph=0, C_CNS=0, Cmet=0)
out_4.5 <- ode(y=y0_4.5, times=t_seq, func=pk_ode, parms=pk_params_opt, method="lsoda")
df_4.5 <- as.data.frame(out_4.5) %>% mutate(C_nM = C * 1000 / MW, C_CNS_nM = C_CNS * 1000 / MW)

cat(sprintf("  50mg Cmax:  %.1f nM at %.1f hr (FDA: ~32 nM at 1 hr)\n",
            max(df_50$C_nM), df_50$time[which.max(df_50$C_nM)]))
cat(sprintf("  4.5mg Cmax: %.2f nM at %.1f hr\n",
            max(df_4.5$C_nM), df_4.5$time[which.max(df_4.5$C_nM)]))
cat(sprintf("  4.5mg CNS Cmax: %.4f nM at %.1f hr\n",
            max(df_4.5$C_CNS_nM), df_4.5$time[which.max(df_4.5$C_CNS_nM)]))

# Simulate 4.5mg steady-state (QD x 28 days)
cat("\n=== Steady-State Simulation: LDN 4.5mg QD x 28 days ===\n")

# Run with daily dosing events
dose_times <- seq(0, 27*24, by=24)
events_df <- data.frame(
  time = dose_times,
  var = "Aa",
  value = 4.5 * 1000,  # ug
  method = "add"
)

y0_ss <- c(Aa=0, C=0, C_periph=0, C_CNS=0, Cmet=0)
out_ss <- ode(y=y0_ss, times=seq(0, 28*24, by=0.5), func=pk_ode,
              parms=pk_params_opt, method="lsoda",
              events=list(data=events_df))

df_ss <- as.data.frame(out_ss) %>%
  mutate(
    C_nM = C * 1000 / MW,
    C_CNS_nM = C_CNS * 1000 / MW,
    Cmet_nM = Cmet * 1000 / MW,
    time_days = time / 24
  )

# Steady-state values
ss_vals <- df_ss %>% filter(time_days >= 27) %>%
  summarise(
    Cmax_nM = max(C_nM),
    Cmin_nM = min(C_nM),
    Cav_nM = mean(C_nM),
    Cmax_CNS_nM = max(C_CNS_nM),
    Cmin_CNS_nM = min(C_CNS_nM),
    Cav_CNS_nM = mean(C_CNS_nM),
    Tmax_hr = time[which.max(C_nM)] %% 24,
    Cmet_nM = mean(Cmet_nM)
  )

cat(sprintf("  Plasma Cmax (ss): %.2f nM\n", ss_vals$Cmax_nM))
cat(sprintf("  Plasma Cmin (ss): %.2f nM\n", ss_vals$Cmin_nM))
cat(sprintf("  Plasma Cav (ss):  %.2f nM\n", ss_vals$Cav_nM))
cat(sprintf("  CNS Cmax (ss):    %.4f nM = %.6f uM\n",
            ss_vals$Cmax_CNS_nM, ss_vals$Cmax_CNS_nM/1000))
cat(sprintf("  CNS Cmin (ss):    %.4f nM\n", ss_vals$Cmin_CNS_nM))
cat(sprintf("  CNS Cav (ss):     %.4f nM = %.6f uM\n",
            ss_vals$Cav_CNS_nM, ss_vals$Cav_CNS_nM/1000))
cat(sprintf("  Metabolite (ss):  %.2f nM\n", ss_vals$Cmet_nM))

# ============================================================
# 7. PLOTS
# ============================================================

# Plot 1: Single dose PK profile (100mg) vs Verebey data
p1 <- ggplot() +
  geom_line(data=df_100, aes(x=time, y=C_nM, color="Model"), linewidth=1.2) +
  geom_point(data=verebey_100mg, aes(x=Time_hr, y=Conc_nM, color="Verebey 1976"),
             size=3) +
  geom_point(data=meyer_100mg, aes(x=Time_hr, y=Conc_nM, color="Meyer 1984"),
             size=3, shape=2) +
  labs(title="Naltrexone PK: 100mg Oral (Model vs Clinical Data)",
       x="Time (hours)", y="Plasma Concentration (nM)", color="Source") +
  theme_minimal() +
  scale_color_manual(values=c("Model"="steelblue","Verebey 1976"="red","Meyer 1984"="darkgreen"))

# Plot 2: Multi-dose PK (LDN 4.5mg QD)
p2 <- ggplot(df_ss, aes(x=time_days, y=C_nM)) +
  geom_line(linewidth=0.8, color="steelblue") +
  labs(title="Naltrexone Plasma Concentration: LDN 4.5mg QD x 28 days",
       x="Time (days)", y="Plasma Concentration (nM)") +
  theme_minimal()

# Plot 3: CNS concentration (LDN 4.5mg QD)
p3 <- ggplot(df_ss, aes(x=time_days, y=C_CNS_nM)) +
  geom_line(linewidth=0.8, color="firebrick") +
  labs(title="Naltrexone CNS Concentration: LDN 4.5mg QD x 28 days",
       x="Time (days)", y="CNS Concentration (nM)") +
  theme_minimal()

# Plot 4: Dose comparison
p4 <- ggplot() +
  geom_line(data=df_100, aes(x=time, y=C_nM, color="100mg"), linewidth=1) +
  geom_line(data=df_50, aes(x=time, y=C_nM, color="50mg"), linewidth=1) +
  geom_line(data=df_4.5, aes(x=time, y=C_nM, color="4.5mg"), linewidth=1) +
  labs(title="Naltrexone PK: Dose Comparison (Single Dose)",
       x="Time (hours)", y="Plasma Concentration (nM)", color="Dose") +
  theme_minimal() +
  scale_color_manual(values=c("100mg"="firebrick","50mg"="orange","4.5mg"="steelblue"))

# Plot 5: TLR4 and OPRM1 occupancy at LDN steady-state
Kd_TLR4 <- 10  # uM
Kd_OPRM1 <- 0.001  # uM

df_occ <- df_ss %>%
  mutate(
    TLR4_occ = (C_CNS_nM/1000) / (Kd_TLR4 + C_CNS_nM/1000) * 100,
    OPRM1_occ = (C_CNS_nM/1000) / (Kd_OPRM1 + C_CNS_nM/1000) * 100
  )

p5 <- ggplot(df_occ %>% filter(time_days >= 25) %>%
               select(time_days, TLR4_occ, OPRM1_occ) %>%
               tidyr::pivot_longer(c(TLR4_occ, OPRM1_occ), names_to="R", values_to="Occ") %>%
               mutate(R=ifelse(R=="TLR4_occ","TLR4 (anti-inflammatory)","OPRM1 (mu-opioid)")),
             aes(x=time_days, y=Occ, color=R)) +
  geom_line(linewidth=1.2) +
  labs(title="Receptor Occupancy at LDN Steady-State",
       subtitle="Based on calibrated PK model",
       x="Time (days)", y="Occupancy (%)", color="Receptor") +
  theme_minimal() +
  scale_color_manual(values=c("TLR4 (anti-inflammatory)"="steelblue",
                               "OPRM1 (mu-opioid)"="firebrick"))

# Save plots
ggsave("../output/figures/pk_validation_100mg.png", p1, width=10, height=6)
ggsave("../output/figures/pk_ldn_plasma.png", p2, width=10, height=6)
ggsave("../output/figures/pk_ldn_cns.png", p3, width=10, height=6)
ggsave("../output/figures/pk_dose_comparison.png", p4, width=10, height=6)
ggsave("../output/figures/pk_receptor_occupancy.png", p5, width=10, height=6)

# ============================================================
# 8. SAVE CALIBRATED PARAMETERS
# ============================================================

calibrated_pk <- data.frame(
  Parameter = c("Ka","Fmax","CL","V1","V2","Q","Kp","CLmet","CLmet_elim",
                "t1/2_naltrexone","t1/2_metabolite","Cmax_4.5mg_nM","CNS_Cmax_4.5mg_nM",
                "TLR4_occ_LDN","OPRM1_occ_LDN"),
  Value = c(pk_params_opt["Ka"], pk_params_opt["Fmax"], pk_params_opt["CL"],
            pk_params_opt["V1"], pk_params_opt["V2"], pk_params_opt["Q"],
            pk_params_opt["Kp"], pk_params_opt["CLmet"], pk_params_opt["CLmet_elim"],
            log(2)/(pk_params_opt["CL"]/pk_params_opt["V1"]),
            log(2)/(pk_params_opt["CLmet_elim"]/pk_params_opt["V1"]),
            max(df_4.5$C_nM), max(df_4.5$C_CNS_nM),
            mean(df_occ$TLR4_occ[df_occ$time_days>=27]),
            mean(df_occ$OPRM1_occ[df_occ$time_days>=27])),
  Units = c("1/hr","fraction","L/hr","L","L","L/hr","ratio","L/hr","L/hr",
            "hr","hr","nM","nM","%","%")
)

write.csv(calibrated_pk, "../output/results/calibrated_pk_parameters.csv", row.names=FALSE)

cat("\n=== PK Calibration Complete ===\n")
cat("Figures saved to: output/figures/pk_*.png\n")
cat("Parameters saved to: output/results/calibrated_pk_parameters.csv\n")
