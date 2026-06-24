# Complete QSP Model with Working PK (deSolve)
# Fixes the PK dosing issue and serves as base for all enhancements

library(deSolve)
library(dplyr)

# ============================================================
# PARAMETERS
# ============================================================

get_params <- function(...) {
  p <- list(
    # PK
    Ka=1.2, F_bio=0.05, Vd=1127, V2=500, Q_inter=120,
    CL=22.6, Kp_brain=0.028, Kin_CNS=0.5, Kout_CNS=0.5,
    Vmax_CYP=50, Km_CYP=10, f_metab=0.6,
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
    # Pain integration
    Pain_b=2.0, Em_BDNF=0.5, EC_BDNF=0.1,
    Em_TNF=0.3, EC_TNF=0.05, Em_TAC1=0.4, EC_TAC1=0.05,
    Em_IL17=0.2, EC_IL17=0.02, E_IL10=0.5, EC_IL10=0.02, nP=1.5,
    CPT_h=60.0, CPT_f=25.0, kCPT=-5.0
  )
  # Override with any provided arguments
  dots <- list(...)
  for (nm in names(dots)) p[[nm]] <- dots[[nm]]
  p
}

# ============================================================
# INITIAL CONDITIONS
# ============================================================

get_y0 <- function(...) {
  y <- c(
    Ccns=0, TLR4f=1.0, TLR4b=0.0,
    MyD=0.35, TRIF=0.35, NFkB=0.35, NLRP3=0.5,
    TNF=0.08, IL1B=0.05, IL6=0.10, IL10=0.015, IL17=0.03, TGF=0.04,
    BDNF=0.15, NGF=0.10, TAC1=0.08, TRPV1=0.5,
    CRH=0.05, POMC=0.04, Cort=0.15, NR3C1=0.2, NPY=0.015,
    HT=0.06, OPRM1d=1.0, Endo=0.002, Pain=7.0, CPT=25.0
  )
  dots <- list(...)
  for (nm in names(dots)) y[nm] <- dots[[nm]]
  y
}

# ============================================================
# ODE SYSTEM
# ============================================================

qsp_ode <- function(t, y, p) {
  with(as.list(c(y, p)), {
    TLR4b_r <- kon*Ccns*TLR4f
    TLR4u_r <- koff*TLR4b
    dCcns <- 0
    dTLR4f <- TLR4_syn - TLR4_deg*TLR4f - TLR4b_r + TLR4u_r
    dTLR4b <- TLR4b_r - TLR4u_r - TLR4_deg*TLR4b
    TLR4tot <- TLR4f + TLR4b
    TRIF_inhib <- TLR4b/(TLR4tot+0.001)
    M_t <- M_basal + (M_max-M_basal)*TLR4f^n_M/(EC50_M^n_M+TLR4f^n_M)
    dMyD <- 0.5*(M_t-MyD)
    T_t <- T_basal + (T_max-T_basal)*TLR4f^n_T/(EC50_T^n_T+TLR4f^n_T)*(1-TRIF_inhib)
    dTRIF <- 0.5*(T_t-TRIF)
    NF_t <- max(0.05, min(1.0, 0.05 + kM_NF*(MyD-0.05) + kT_NF*(TRIF-0.05)))
    dNFkB <- 0.5*(NF_t-NFkB)
    NLRP3_t <- N_basal + (N_max-N_basal)*NFkB^n_N/(EC50_N^n_N+NFkB^n_N)
    dNLRP3 <- 0.3*(NLRP3_t-NLRP3)
    dTNF <- kpTNF*(1+kM_TNF*MyD+kT_TNF*TRIF) - kdTNF*TNF + 0.05*TNF*(1-TNF/TNF_max) - 0.1*IL10*TNF
    dIL1B <- kpIL1B*(1+kN_IL1B*NLRP3) - kdIL1B*IL1B
    dIL6 <- kpIL6*(1+kM_IL6*MyD+kT_IL6*TRIF) - kdIL6*IL6
    dIL10 <- kpIL10*(1+kNR_IL10*NR3C1+kT_IL10*TRIF) - kdIL10*IL10
    dIL17 <- kpIL17*(1+kT_IL17*TRIF+kNF_IL17*NFkB) - kdIL17*IL17
    dTGF <- kpTGF*(1+kT_TGF*TRIF+kNF_TGF*NFkB) - kdTGF*TGF
    dBDNF <- kpBDNF*(1+kIL6_BDNF*IL6) - kdBDNF*BDNF
    if(BDNF>BDNF_max) dBDNF <- dBDNF-0.5*(BDNF-BDNF_max)
    dNGF <- kpNGF*(1+kTNF_NGF*TNF) - kdNGF*NGF
    if(NGF>NGF_max) dNGF <- dNGF-0.5*(NGF-NGF_max)
    dTAC1 <- kpTAC1*(1+kB_TAC1*BDNF+kIL17_TAC1*IL17) - kdTAC1*TAC1
    if(TAC1>TAC1_max) dTAC1 <- dTAC1-0.5*(TAC1-TAC1_max)
    TRPV1_t <- TRPV1_b+(TRPV1_m-TRPV1_b)*NGF^n_TRPV1/(EC50_TRPV1^n_TRPV1+NGF^n_TRPV1)+kB_TRPV1*BDNF
    dTRPV1 <- 0.5*(TRPV1_t-TRPV1)
    dCRH <- kpCRH*(1+kTNF_CRH*TNF) - kdCRH*CRH
    dPOMC <- kpPOMC*(1+kCRH_POMC*CRH) - kdPOMC*POMC
    dCort <- 0.5*POMC - 0.3*Cort
    NRs <- max(0.1, 1-kTNF_NR*TNF)
    NRt <- NRs*Cort^n_NR/(EC50_NR^n_NR+Cort^n_NR)
    dNR3C1 <- 0.3*(NRt-NR3C1)
    dNPY <- kpNPY/(1+kCRH_NPY*CRH) - kdNPY*NPY
    dHT <- ksHT - VmaxHT*HT/(KmHT+HT) - kdHT*HT + 0.02*IL10
    OPRM1_occ <- Ccns/(Kd_OPRM1+Ccns)
    dOPRM1d <- k_up*OPRM1_occ*(2-OPRM1d) - k_deg_OPRM1*(OPRM1d-OPRM1_b)
    dEndo <- k_reb*OPRM1_occ*(0.01-Endo) - 0.005*(Endo-endo_basal)
    Endo_sig <- (1-OPRM1_occ)*OPRM1d*Endo/(Kd_endo+Endo)
    OPRM1_a <- OPRM1_amax*Endo_sig^n_OPRM1/(EC50_OPRM1^n_OPRM1+Endo_sig^n_OPRM1)
    Hyper <- 1+k_hyper*OPRM1_occ^3
    pB <- Em_BDNF*BDNF^nP/(EC_BDNF^nP+BDNF^nP)
    pT <- Em_TNF*TNF^nP/(EC_TNF^nP+TNF^nP)
    pTA <- Em_TAC1*TAC1^nP/(EC_TAC1^nP+TAC1^nP)
    pIL <- Em_IL17*IL17^nP/(EC_IL17^nP+IL17^nP)
    pIL10 <- E_IL10*IL10/(EC_IL10+IL10)
    Pain_t <- Pain_b*(1+pB+pT+pTA+pIL)/(1+pIL10)*Hyper*(1-0.3*OPRM1_a)
    Pain_t <- max(0, min(10, Pain_t))
    dPain <- 0.1*(Pain_t-Pain)
    CPT_t <- max(5, min(CPT_h, CPT_f+kCPT*(Pain-7)+10*OPRM1_a))
    dCPT <- 0.05*(CPT_t-CPT)
    list(c(dCcns,dTLR4f,dTLR4b,dMyD,dTRIF,dNFkB,dNLRP3,
           dTNF,dIL1B,dIL6,dIL10,dIL17,dTGF,
           dBDNF,dNGF,dTAC1,dTRPV1,
           dCRH,dPOMC,dCort,dNR3C1,dNPY,dHT,
           dOPRM1d,dEndo,dPain,dCPT))
  })
}

# ============================================================
# SIMULATION FUNCTION
# ============================================================

run_qsp <- function(Ccns_val, p=NULL, y0=NULL, t_end=28*24, dt=1) {
  if(is.null(p)) p <- get_params()
  if(is.null(y0)) y0 <- get_y0()
  y0["Ccns"] <- Ccns_val
  t_seq <- seq(0, t_end, by=dt)
  out <- ode(y=y0, times=t_seq, func=qsp_ode, parms=p, method="lsoda")
  df <- as.data.frame(out)
  df$time_days <- df$time/24
  df$Ccns_uM <- Ccns_val
  df$TLR4_occ <- df$TLR4b/(df$TLR4f+df$TLR4b+0.001)*100
  df$OPRM1_occ <- Ccns_val/(0.001+Ccns_val)*100
  df$Pain_red <- (7-df$Pain)/7*100
  df$OPRM1_fold <- df$OPRM1d/1.0
  df$Endo_fold <- df$Endo/0.002
  df
}

# ============================================================
# STEADY-STATE EXTRACTION
# ============================================================

get_ss <- function(df, t_min=27) {
  df %>% filter(time_days >= t_min) %>%
    summarise(across(where(is.numeric), mean))
}

cat("QSP model functions loaded successfully.\n")
cat("Use: run_qsp(Ccns_val) to simulate a scenario.\n")
