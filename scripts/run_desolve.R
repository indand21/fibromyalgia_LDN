# QSP Model v4 - With placebo response mechanism
# Adds expectation-mediated analgesia for all treatment arms
# Literature: Bingel et al. 2011, Colloca et al. 2008

library(deSolve)
library(dplyr)
library(ggplot2)

# ============================================================
# 1. PARAMETERS
# ============================================================

p <- c(
  # TLR4
  kon=0.1, koff=1.0, TLR4_syn=0.01, TLR4_deg=0.01,
  M_basal=0.05, M_max=1.0, EC50_M=0.5, n_M=2.0,
  T_basal=0.05, T_max=1.0, EC50_T=0.3, n_T=2.0,
  kM_NF=0.7, kT_NF=0.3,
  N_basal=0.1, N_max=1.0, EC50_N=0.3, n_N=1.5,
  # Cytokines
  kpTNF=0.5, kdTNF=0.693, kM_TNF=1.0, kT_TNF=1.5, TNF_max=1.0,
  kpIL1B=0.3, kdIL1B=0.347, kN_IL1B=3.0,
  kpIL6=0.4, kdIL6=0.116, kM_IL6=1.5, kT_IL6=1.0,
  kpIL10=0.2, kdIL10=0.116, kNR_IL10=1.5, kT_IL10=0.5,
  kpIL17=0.08, kdIL17=0.116, kT_IL17=1.0, kNF_IL17=0.5,
  kpTGF=0.06, kdTGF=0.058, kT_TGF=0.8, kNF_TGF=0.3,
  # Pain
  kpBDNF=0.1, kdBDNF=0.023, BDNF_max=1.0, kIL6_BDNF=0.5,
  kpNGF=0.08, kdNGF=0.023, NGF_max=1.0, kTNF_NGF=0.3,
  kpTAC1=0.15, kdTAC1=0.693, TAC1_max=1.0, kB_TAC1=1.0, kIL17_TAC1=0.3,
  TRPV1_b=0.1, TRPV1_m=1.0, EC50_TRPV1=0.05, n_TRPV1=2.0, kB_TRPV1=0.5,
  # HPA
  kpCRH=0.1, kdCRH=0.693, kTNF_CRH=0.2,
  kpPOMC=0.1, kdPOMC=0.347, kCRH_POMC=2.0,
  EC50_NR=0.1, n_NR=1.5, kTNF_NR=0.5,
  kpNPY=0.05, kdNPY=0.116, kCRH_NPY=0.3,
  # Monoamines
  ksHT=0.1, kdHT=0.5, VmaxHT=1.0, KmHT=0.5,
  # OPRM1 hormesis
  Kd_OPRM1=0.001, Kd_endo=0.005, endo_basal=0.002,
  OPRM1_b=1.0, k_up=0.05, k_deg_OPRM1=0.01,
  k_reb=0.02, OPRM1_amax=0.8, EC50_OPRM1=0.3, n_OPRM1=1.5, k_hyper=0.5,
  # OGF-OGFr
  Kd_OGFr=0.00001, k_OGF_prod=0.000005, k_OGF_deg=0.05,
  k_OGF_upreg=2000.0, EC50_OGF=0.00005, OGF=0.005,
  # Pain integration
  Pain_b=2.0, Em_BDNF=0.5, EC_BDNF=0.1,
  Em_TNF=0.3, EC_TNF=0.05, Em_TAC1=0.4, EC_TAC1=0.05,
  Em_IL17=0.2, EC_IL17=0.02, E_IL10=0.5, EC_IL10=0.02, nP=1.5,
  CPT_h=60.0, CPT_f=25.0, kCPT=-5.0,
  # Positive feedback loops
  k_M_act_prod=0.05, k_M_act_pain=0.02, k_M_act_deg=0.005,
  k_M_act_self=0.01, M_act_max=2.0,
  k_M_act_TNF=1.0, k_M_act_IL6=0.5, k_M_act_IL1B=0.8,
  k_pain_TNF=0.05, k_pain_IL6=0.02,
  k_central_sens=0.5, EC50_sens=3.0, n_sens=3.0,
  # Chronic pain state
  k_chronic_prod=0.08, k_chronic_deg=0.03, k_chronic_self=0.0,
  EC50_chronic=4.0, n_chronic=2.0, k_chronic_amp=0.5,
  # ============================================================
  # Placebo response mechanism
  # Literature: Bingel et al. 2011 (Lancet Neurol), Colloca et al. 2008
  # Expectation-mediated analgesia via endogenous opioid release
  # ============================================================
  k_placebo_prod=0.005,   # placebo response development rate (slower)
  k_placebo_deg=0.002,    # placebo response decay rate
  placebo_max=0.10,       # maximum placebo analgesic effect (10% pain reduction)
  EC50_placebo=1.0,       # half-maximum time (days) for placebo development
  n_placebo=2.0,          # Hill coefficient for placebo development
  k_placebo_endo=0.1,     # placebo-mediated endorphin release (reduced)
  k_placebo_IL10=0.05,    # placebo-mediated IL-10 increase (reduced)
  k_placebo_cort=0.02     # placebo-mediated cortisol reduction (reduced)
)

# ============================================================
# 2. INITIAL CONDITIONS (FMQ baseline with microglial activation)
# ============================================================

y0 <- c(
  Ccns=0, TLR4f=1.0, TLR4b=0.0,
  MyD=0.35, TRIF=0.35, NFkB=0.35, NLRP3=0.5,
  TNF=0.08, IL1B=0.05, IL6=0.10, IL10=0.015, IL17=0.03, TGF=0.04,
  BDNF=0.15, NGF=0.10, TAC1=0.08, TRPV1=0.5,
  CRH=0.05, POMC=0.04, Cort=0.15, NR3C1=0.2, NPY=0.015,
  HT=0.06, OPRM1d=1.0, Endo=0.002,
  OGF=0.005,
  M_act=0.5,           # elevated microglial activation in FMQ
  Chronic=0.5,         # elevated chronic pain state in FMQ
  PlaceboResp=0.0,     # NEW: placebo response (develops over time)
  Pain=7.0, CPT=25.0
)

# ============================================================
# 3. ODE SYSTEM WITH POSITIVE FEEDBACK LOOPS + PLACEBO RESPONSE
# ============================================================

qsp <- function(t, y, p) {
  with(as.list(c(y, p)), {
    # TLR4 binding
    TLR4b_rate <- kon*Ccns*TLR4f
    TLR4u_rate <- koff*TLR4b
    dCcns <- 0
    dTLR4f <- TLR4_syn - TLR4_deg*TLR4f - TLR4b_rate + TLR4u_rate
    dTLR4b <- TLR4b_rate - TLR4u_rate - TLR4_deg*TLR4b
    TLR4tot <- TLR4f + TLR4b
    TRIF_inhib <- TLR4b/(TLR4tot+0.001)

    # MyD88 (NOT blocked by naltrexone)
    M_t <- M_basal + (M_max-M_basal)*TLR4f^n_M/(EC50_M^n_M+TLR4f^n_M)
    dMyD <- 0.5*(M_t-MyD)

    # TRIF (BLOCKED by naltrexone via MD2)
    T_t <- T_basal + (T_max-T_basal)*TLR4f^n_T/(EC50_T^n_T+TLR4f^n_T)*(1-TRIF_inhib)
    dTRIF_ <- 0.5*(T_t-TRIF)

    # NFkB
    NF_t <- 0.05 + kM_NF*(MyD-0.05) + kT_NF*(TRIF-0.05)
    NF_t <- max(0.05, min(1.0, NF_t))
    dNFkB <- 0.5*(NF_t-NFkB)

    # NLRP3
    NLRP3_t <- N_basal + (N_max-N_basal)*NFkB^n_N/(EC50_N^n_N+NFkB^n_N)
    dNLRP3 <- 0.3*(NLRP3_t-NLRP3)

    # ============================================================
    # NEW: Microglial self-activation (persistent neuroinflammation)
    # Literature: Wakatsuki 2024, Chen 2022, Ji 2018
    # Key: M_act has SLOW degradation -> persistence
    # Key: M_act self-amplifies -> bistability
    # ============================================================
    # Activation driven by: inflammation + pain + self-amplification
    M_act_drive <- k_M_act_prod*(TNF + IL6 + IL1B) +  # inflammation-driven
                   k_M_act_pain*Pain +                   # pain-driven
                   k_M_act_self*M_act*(M_act_max - M_act)  # self-amplification
    # Deactivation (slow -> persistence)
    M_act_res <- k_M_act_deg*M_act
    dM_act <- M_act_drive - M_act_res
    # Clamp to [0, M_act_max]
    if(M_act < 0) dM_act <- max(0, dM_act)
    if(M_act > M_act_max) dM_act <- min(0, dM_act)

    # Microglial amplification factor (0 = no amplification, 1 = max)
    M_amp <- M_act / M_act_max

    # ============================================================
    # Cytokines WITH microglial amplification + pain feedback
    # NO cooperative self-amplification (removed to prevent runaway)
    # ============================================================
    # Pain-inflammation positive feedback (vicious cycle)
    pain_fb <- k_pain_TNF * max(0, Pain - 5)  # only above pain threshold

    dTNF <- kpTNF*(1 + kM_TNF*MyD + kT_TNF*TRIF +
                   k_M_act_TNF*M_amp +       # microglial amplification
                   pain_fb) -                 # pain-inflammation feedback
           kdTNF*TNF - 0.1*IL10*TNF

    dIL1B_ <- kpIL1B*(1 + kN_IL1B*NLRP3 + k_M_act_IL1B*M_amp) - kdIL1B*IL1B

    dIL6 <- kpIL6*(1 + kM_IL6*MyD + kT_IL6*TRIF +
                   k_M_act_IL6*M_amp +
                   k_pain_IL6*max(0, Pain - 5)) - kdIL6*IL6

    dIL10_ <- kpIL10*(1 + kNR_IL10*NR3C1 + kT_IL10*TRIF) - kdIL10*IL10
    dIL17_ <- kpIL17*(1 + kT_IL17*TRIF + kNF_IL17*NFkB) - kdIL17*IL17
    dTGF_ <- kpTGF*(1 + kT_TGF*TRIF + kNF_TGF*NFkB) - kdTGF*TGF

    # Pain signaling
    dBDNF <- kpBDNF*(1+kIL6_BDNF*IL6) - kdBDNF*BDNF
    if(BDNF>BDNF_max) dBDNF <- dBDNF-0.5*(BDNF-BDNF_max)
    dNGF <- kpNGF*(1+kTNF_NGF*TNF) - kdNGF*NGF
    if(NGF>NGF_max) dNGF <- dNGF-0.5*(NGF-NGF_max)
    dTAC1_ <- kpTAC1*(1+kB_TAC1*BDNF+kIL17_TAC1*IL17) - kdTAC1*TAC1
    if(TAC1>TAC1_max) dTAC1_ <- dTAC1_-0.5*(TAC1-TAC1_max)
    TRPV1_t <- TRPV1_b+(TRPV1_m-TRPV1_b)*NGF^n_TRPV1/(EC50_TRPV1^n_TRPV1+NGF^n_TRPV1)+kB_TRPV1*BDNF
    dTRPV1_ <- 0.5*(TRPV1_t-TRPV1)

    # HPA
    dCRH <- kpCRH*(1+kTNF_CRH*TNF) - kdCRH*CRH
    dPOMC_ <- kpPOMC*(1+kCRH_POMC*CRH) - kdPOMC*POMC
    dCort <- 0.5*POMC - 0.3*Cort
    NRs <- max(0.1, 1-kTNF_NR*TNF)
    NRt <- NRs*Cort^n_NR/(EC50_NR^n_NR+Cort^n_NR)
    dNR3C1_ <- 0.3*(NRt-NR3C1)
    dNPY <- kpNPY/(1+kCRH_NPY*CRH) - kdNPY*NPY

    # Monoamines
    dHT <- ksHT - VmaxHT*HT/(KmHT+HT) - kdHT*HT + 0.02*IL10

    # OGF-OGFr pathway
    OGFr_occ <- Ccns/(Kd_OGFr+Ccns)
    dOGF <- k_OGF_prod*(1 + k_OGF_upreg*OGFr_occ) - k_OGF_deg*OGF
    OGF_eff <- OGF/(OGF + EC50_OGF)

    # OPRM1 hormesis
    OPRM1_occ <- Ccns/(Kd_OPRM1+Ccns)
    dOPRM1d <- k_up*OPRM1_occ*(2-OPRM1d) - k_deg_OPRM1*(OPRM1d-OPRM1_b)
    dEndo <- k_reb*OPRM1_occ*(0.01-Endo) - 0.005*(Endo-endo_basal)
    Endo_sig <- (1-OPRM1_occ)*OPRM1d*Endo/(Kd_endo+Endo)
    OPRM1_a <- OPRM1_amax*Endo_sig^n_OPRM1/(EC50_OPRM1^n_OPRM1+Endo_sig^n_OPRM1)
    Hyper <- 1+k_hyper*OPRM1_occ^3

    anti_inflam <- OGF_eff + 0.1*OPRM1_a

    # ============================================================
    # Placebo response mechanism
    # Literature: Bingel et al. 2011, Colloca et al. 2008
    # Expectation-mediated analgesia via endogenous opioid release
    # Develops slowly over 2-4 weeks, persists for weeks
    # ============================================================
    # Placebo response development (time-dependent, cooperative)
    placebo_drive <- k_placebo_prod * t^n_placebo / (EC50_placebo^n_placebo + t^n_placebo)
    # Decay (very slow -> persistence)
    placebo_res <- k_placebo_deg * PlaceboResp
    dPlaceboResp <- placebo_drive - placebo_res
    # Clamp to [0, placebo_max]
    if(PlaceboResp > placebo_max) dPlaceboResp <- min(0, dPlaceboResp)

    # Placebo-mediated effects
    # 1. Endogenous opioid release (expectation analgesia)
    placebo_endo <- k_placebo_endo * PlaceboResp
    # 2. IL-10 increase (anti-inflammatory)
    placebo_IL10 <- k_placebo_IL10 * PlaceboResp
    # 3. Cortisol reduction (stress relief)
    placebo_cort <- k_placebo_cort * PlaceboResp

    # Placebo analgesic effect (reduces pain)
    placebo_analgesia <- PlaceboResp  # direct pain reduction

    # ============================================================
    # Chronic pain state (neural plasticity / central sensitization)
    # ============================================================
    # Production: driven by current pain level (cooperative)
    chronic_drive <- k_chronic_prod * Pain^n_chronic / (EC50_chronic^n_chronic + Pain^n_chronic)
    # Self-amplification (cooperative positive feedback)
    chronic_self <- k_chronic_self * Chronic^2 / (0.3^2 + Chronic^2)
    # Degradation (slow -> persistence)
    chronic_res <- k_chronic_deg * Chronic
    dChronic <- chronic_drive + chronic_self - chronic_res

    # Chronic state amplification factor
    chronic_amp <- 1.0 + k_chronic_amp * Chronic

    # ============================================================
    # Pain WITH central sensitization + chronic state amplification
    # ============================================================
    pB <- Em_BDNF*BDNF^nP/(EC_BDNF^nP+BDNF^nP)
    pT <- Em_TNF*TNF^nP/(EC_TNF^nP+TNF^nP)
    pTA <- Em_TAC1*TAC1^nP/(EC_TAC1^nP+TAC1^nP)
    pIL <- Em_IL17*IL17^nP/(EC_IL17^nP+IL17^nP)
    pIL10 <- E_IL10*IL10/(EC_IL10+IL10)

    # Central sensitization + chronic state amplification
    central_amp <- 1.0 + k_central_sens * M_amp *
                   Pain^n_sens / (EC50_sens^n_sens + Pain^n_sens)

    Pain_t <- Pain_b*(1+pB+pT+pTA+pIL)/(1+pIL10+placebo_IL10)*Hyper*(1-0.3*anti_inflam-placebo_analgesia)*central_amp*chronic_amp
    Pain_t <- max(0, min(10, Pain_t))
    dPain <- 0.1*(Pain_t-Pain)
    CPT_t <- max(5, min(CPT_h, CPT_f+kCPT*(Pain-7)+10*anti_inflam))
    dCPT_ <- 0.05*(CPT_t-CPT)

    list(c(dCcns,dTLR4f,dTLR4b,dMyD,dTRIF_,dNFkB,dNLRP3,
           dTNF,dIL1B_,dIL6,dIL10_,dIL17_,dTGF_,
           dBDNF,dNGF,dTAC1_,dTRPV1_,
           dCRH,dPOMC_,dCort,dNR3C1_,dNPY,dHT,
           dOPRM1d,dEndo,dOGF,
           dM_act,dChronic,dPlaceboResp,  # microglial + chronic + placebo
           dPain,dCPT_))
  })
}

# ============================================================
# 4. RUN SCENARIOS
# ============================================================

cat("=== QSP Model v4: With Placebo Response Mechanism ===\n")
cat("Adds expectation-mediated analgesia for all treatment arms\n\n")

# Calibrated CNS concentrations (from PK calibration)
scenarios <- data.frame(
  name = c("Placebo", "LDN 1.5mg", "LDN 3.0mg", "LDN 4.5mg", "LDN 6.0mg", "Standard 50mg"),
  Ccns = c(0, 0.000013, 0.000026, 0.000039, 0.000052, 0.00137)
)

t_seq <- seq(0, 28*24, by=1)

run_sim <- function(Ccns_val, name) {
  cat(sprintf("  %s (Ccns = %.6f uM)\n", name, Ccns_val))
  y0_run <- y0
  y0_run["Ccns"] <- Ccns_val
  out <- ode(y=y0_run, times=t_seq, func=qsp, parms=p, method="lsoda")
  df <- as.data.frame(out) %>%
    mutate(time_days=time/24, Scenario=name, Ccns_uM=Ccns_val) %>%
    mutate(
      TLR4_occ=TLR4b/(TLR4f+TLR4b+0.001)*100,
      OPRM1_occ=Ccns_uM/(0.001+Ccns_uM)*100,
      Pain_red=(7-Pain)/7*100,
      OPRM1_fold=OPRM1d/1.0,
      Endo_fold=Endo/0.002,
      M_act_norm=M_act/1.0,
      Chronic_norm=Chronic/1.0,
      PlaceboResp_norm=PlaceboResp/1.0
    )
  df
}

results <- list(
  run_sim(0, "Placebo"),
  run_sim(0.000013, "LDN 1.5mg"),
  run_sim(0.000026, "LDN 3.0mg"),
  run_sim(0.000039, "LDN 4.5mg"),
  run_sim(0.000052, "LDN 6.0mg"),
  run_sim(0.00137, "Standard 50mg")
)

all <- bind_rows(results)

# ============================================================
# 5. STEADY-STATE
# ============================================================

cat("\n=== Steady-State (Day 28) ===\n")

ss <- all %>%
  filter(time_days >= 27) %>%
  group_by(Scenario, Ccns_uM) %>%
  summarise(
    TLR4_occ=mean(TLR4_occ), OPRM1_occ=mean(OPRM1_occ),
    Pain_VAS=mean(Pain), Pain_red=mean(Pain_red), CPT=mean(CPT),
    TNF=mean(TNF), IL1B=mean(IL1B), IL6=mean(IL6), IL10=mean(IL10),
    IL17=mean(IL17), TGF=mean(TGF),
    BDNF=mean(BDNF), TAC1=mean(TAC1),
    OPRM1_fold=mean(OPRM1_fold), Endo_fold=mean(Endo_fold),
    M_act=mean(M_act_norm), Chronic=mean(Chronic_norm),
    PlaceboResp=mean(PlaceboResp_norm),
    .groups="drop"
  )

print(as.data.frame(ss))

# ============================================================
# 6. PLOTS
# ============================================================

cat("\nGenerating figures...\n")

p1 <- ggplot(all, aes(x=time_days, y=Pain, color=Scenario)) +
  geom_line(linewidth=1.2) +
  labs(title="Pain VAS: With Positive Feedback Loops",
       subtitle="QSP Model v3: Chronic disease state maintained",
       x="Time (days)", y="Pain VAS (0-10)") +
  theme_minimal() +
  scale_color_manual(values=c("Placebo"="gray50","LDN 1.5mg"="lightblue",
                               "LDN 3.0mg"="deepskyblue","LDN 4.5mg"="steelblue",
                               "LDN 6.0mg"="darkblue","Standard 50mg"="firebrick"))

p2 <- ggplot(all, aes(x=time_days, y=M_act_norm, color=Scenario)) +
  geom_line(linewidth=1.2) +
  labs(title="Microglial Activation: Self-Sustaining State",
       subtitle="Positive feedback maintains chronic activation",
       x="Time (days)", y="Microglial Activation (normalized)") +
  theme_minimal() +
  scale_color_manual(values=c("Placebo"="gray50","LDN 4.5mg"="steelblue","Standard 50mg"="firebrick"))

cyto <- all %>%
  filter(Scenario %in% c("Placebo","LDN 4.5mg","Standard 50mg")) %>%
  select(time_days, Scenario, TNF, IL6, IL10, IL17) %>%
  tidyr::pivot_longer(c(TNF, IL6, IL10, IL17), names_to="Cytokine", values_to="Conc")

p3 <- ggplot(cyto, aes(x=time_days, y=Conc, color=Scenario)) +
  geom_line(linewidth=0.8) + facet_wrap(~Cytokine, scales="free_y", ncol=2) +
  labs(title="Cytokines: Persistent Disease State", x="Time (days)", y="uM") +
  theme_minimal() +
  scale_color_manual(values=c("Placebo"="gray50","LDN 4.5mg"="steelblue","Standard 50mg"="firebrick"))

p4 <- ss %>%
  ggplot(aes(x=Scenario, y=Pain_red, fill=Scenario)) +
  geom_bar(stat="identity", width=0.6) + geom_hline(yintercept=0) +
  labs(title="Pain Reduction: With Positive Feedback", x="", y="Pain Reduction (%)") +
  theme_minimal() +
  scale_fill_manual(values=c("Placebo"="gray","LDN 1.5mg"="lightblue","LDN 3.0mg"="deepskyblue",
                               "LDN 4.5mg"="darkblue","LDN 6.0mg"="navy","Standard 50mg"="firebrick")) +
  theme(legend.position="none", axis.text.x=element_text(angle=30,hjust=1))

p5 <- ss %>%
  select(Scenario, TLR4_occ, OPRM1_occ) %>%
  tidyr::pivot_longer(c(TLR4_occ, OPRM1_occ), names_to="R", values_to="Occ") %>%
  mutate(R=ifelse(R=="TLR4_occ","TLR4","OPRM1")) %>%
  ggplot(aes(x=Scenario, y=Occ, fill=R)) +
  geom_bar(stat="identity", position="dodge", width=0.7) +
  labs(title="Receptor Occupancy", x="", y="%") +
  theme_minimal() +
  scale_fill_manual(values=c("TLR4"="steelblue","OPRM1"="firebrick")) +
  theme(axis.text.x=element_text(angle=30,hjust=1))

# ============================================================
# 7. SAVE
# ============================================================

ggsave("../output/figures/pain_trajectory_v4.png", p1, width=12, height=7)
ggsave("../output/figures/microglial_activation.png", p2, width=10, height=6)
ggsave("../output/figures/cytokines_v4.png", p3, width=14, height=10)
ggsave("../output/figures/dose_response_v4.png", p4, width=10, height=7)
ggsave("../output/figures/receptor_occupancy_v4.png", p5, width=10, height=7)

write.csv(all, "../output/results/full_simulation_v4.csv", row.names=FALSE)
write.csv(as.data.frame(ss), "../output/results/steady_state_v4.csv", row.names=FALSE)

# ============================================================
# 8. VALIDATION
# ============================================================

ldn <- ss %>% filter(Scenario=="LDN 4.5mg")
plac <- ss %>% filter(Scenario=="Placebo")
std <- ss %>% filter(Scenario=="Standard 50mg")

cat("\n================================================================\n")
cat("  QSP MODEL v4: WITH PLACEBO RESPONSE MECHANISM\n")
cat("================================================================\n")
cat(sprintf("\n  Disease persistence check:\n"))
cat(sprintf("    Placebo Pain VAS Day 0:   7.0\n"))
cat(sprintf("    Placebo Pain VAS Day 28:  %.2f %s\n", plac$Pain_VAS,
            ifelse(plac$Pain_VAS > 5.0, "[PERSISTENT - GOOD]", "[SELF-RESOLVING - BAD]")))
cat(sprintf("    Placebo response:         %.3f (%.1f%% pain reduction)\n",
            plac$PlaceboResp, plac$PlaceboResp*100))
cat(sprintf("    Placebo TNF Day 28:       %.4f\n", plac$TNF))
cat(sprintf("    Placebo M_act Day 28:     %.2f\n", plac$M_act))
cat(sprintf("\n  LDN 4.5mg effects:\n"))
cat(sprintf("    Pain reduction from baseline: %.1f%%\n", ldn$Pain_red))
cat(sprintf("    Pain reduction vs placebo:    %.1f%%\n",
            (plac$Pain_VAS - ldn$Pain_VAS)/plac$Pain_VAS*100))
cat(sprintf("    CPT improvement: %.1f sec\n", ldn$CPT))
cat(sprintf("    OPRM1 fold: %.2fx\n", ldn$OPRM1_fold))
cat(sprintf("    Endorphin fold: %.2fx\n", ldn$Endo_fold))
cat(sprintf("    Microglial activation: %.2f\n", ldn$M_act))
cat(sprintf("\n  Standard 50mg effects:\n"))
cat(sprintf("    Pain reduction from baseline: %.1f%%\n", std$Pain_red))
cat(sprintf("    Pain reduction vs placebo:    %.1f%%\n",
            (plac$Pain_VAS - std$Pain_VAS)/plac$Pain_VAS*100))
cat(sprintf("\n  Clinical comparison (Younger 2013):\n"))
cat(sprintf("    Model LDN vs Placebo:   %.1f%%\n",
            (plac$Pain_VAS - ldn$Pain_VAS)/plac$Pain_VAS*100))
cat(sprintf("    Clinical LDN vs Placebo: 10.8%%\n"))
cat(sprintf("    Model Placebo response:  %.1f%%\n", plac$PlaceboResp*100))
cat(sprintf("    Clinical Placebo response: ~18%%\n"))
cat("================================================================\n")
cat("Output: ../output/figures/*_v4.png, ../output/results/*_v4.csv\n")
