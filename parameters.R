# QSP Model Parameters: Naltrexone in Fibromyalgia
# Sources: Literature, GeneCards, STRING database, molecular docking

# ============================================================
# PHARMACOKINETIC PARAMETERS (Naltrexone)
# ============================================================

pk_params <- list(
  # Oral absorption
  Ka       = 1.2,        # 1/hr, first-order absorption rate (FDA label)
  F        = 0.05,       # fractional bioavailability at low dose (5% at 4.5mg)
  Tlag     = 0.25,       # hr, absorption lag time

  # Distribution
  Vd       = 16.1,       # L/kg * 70kg = 1127 L (total body, from literature)
  V2       = 500,        # L, peripheral compartment volume
  Q        = 120,        # L/hr, inter-compartmental clearance

  # CNS penetration
  Kp_brain = 0.028,      # brain:plasma partition coefficient (low for naltrexone)
  Kin_CNS  = 0.5,        # 1/hr, CNS influx rate
  Kout_CNS = 0.5,        # 1/hr, CNS efflux rate

  # Elimination
  CL       = 22.6,       # L/hr, systemic clearance (from FDA label, adjusted for LDN)
  CLrenal  = 2.0,        # L/hr, renal clearance component

  # Metabolism (CYP2D6 -> 6-beta-naltrexol)
  Vmax_CYP = 50,         # nmol/hr/mg protein
  Km_CYP   = 10,         # uM, CYP2D6 Km for naltrexone
  f_metabolite = 0.6     # fraction converted to active metabolite
)

# ============================================================
# TLR4 BINDING & SIGNALING PARAMETERS
# ============================================================

tlr4_params <- list(
  # Naltrexone-TLR4 binding (from molecular docking)
  # Literature: Kd ~1-10 uM (micromolar, LOW affinity compared to OPRM1)
  # This differential affinity is WHY low dose works:
  #   - At 4.5mg: C_CNS ~0.01 uM -> partial TLR4 occupancy, negligible OPRM1 blockade
  #   - At 50mg:  C_CNS ~1 uM    -> full TLR4 + full OPRM1 blockade -> opioid side effects
  kon_NTX_TLR4  = 0.1,     # 1/(uM*hr), association rate constant
  koff_NTX_TLR4 = 1.0,     # 1/hr, dissociation rate constant
  Kd_NTX_TLR4   = 10,      # uM, equilibrium dissociation constant (Kd = koff/kon)

  # TLR4 receptor dynamics
  TLR4_total    = 1.0,      # uM, total TLR4 expression (microglia, normalized)
  TLR4_syn      = 0.01,     # uM/hr, TLR4 synthesis rate
  TLR4_deg      = 0.01,     # 1/hr, TLR4 degradation rate

  # NF-kB signaling (downstream of TLR4)
  NFkB_basal    = 0.05,     # basal NF-kB activity (normalized)
  NFkB_max      = 1.0,      # maximum NF-kB activity
  EC50_TLR4_NFkB= 0.5,      # uM, EC50 for TLR4->NF-kB activation
  n_TLR4_NFkB   = 2.0,      # Hill coefficient

  # NLRP3 inflammasome priming
  NLRP3_basal   = 0.1,      # basal NLRP3 activity
  NLRP3_max     = 1.0,      # max NLRP3 activity
  EC50_NFkB_NLRP3 = 0.3,    # EC50 for NF-kB -> NLRP3 priming
  n_NLRP3       = 1.5       # Hill coefficient
)

# ============================================================
# CYTOKINE PARAMETERS
# ============================================================

cytokine_params <- list(
  # TNF (degree = 45, hub gene)
  k_prod_TNF    = 0.5,      # uM/hr, basal production rate
  k_deg_TNF     = 0.693,    # 1/hr, degradation (t1/2 ~1 hr)
  TNF_basal     = 0.01,     # uM, healthy baseline
  TNF_FMQ       = 0.08,     # uM, fibromyalgia baseline (elevated)
  k_NFkB_TNF    = 2.0,      # NF-kB driven production multiplier
  TNF_max       = 1.0,      # uM, saturation level
  k_auto_TNF    = 0.1,      # autocrine positive feedback

  # IL1B (degree = 42, hub gene)
  k_prod_IL1B   = 0.3,      # uM/hr
  k_deg_IL1B    = 0.347,    # 1/hr (t1/2 ~2 hr)
  IL1B_basal    = 0.005,    # uM, healthy baseline
  IL1B_FMQ      = 0.05,     # uM, fibromyalgia baseline
  k_NLRP3_IL1B  = 3.0,      # NLRP3-driven production multiplier
  IL1B_max      = 1.0,      # uM

  # IL6 (degree = 45, hub gene)
  k_prod_IL6    = 0.4,      # uM/hr
  k_deg_IL6     = 0.116,    # 1/hr (t1/2 ~6 hr)
  IL6_basal     = 0.01,     # uM, healthy baseline
  IL6_FMQ       = 0.1,      # uM, fibromyalgia baseline
  k_NFkB_IL6    = 2.5,      # NF-kB driven production multiplier
  IL6_max       = 1.0,      # uM

  # IL10 (degree = 36, anti-inflammatory)
  k_prod_IL10   = 0.2,      # uM/hr
  k_deg_IL10    = 0.116,    # 1/hr (t1/2 ~6 hr)
  IL10_basal    = 0.02,     # uM, healthy baseline
  IL10_FMQ      = 0.015,    # uM, fibromyalgia baseline (suppressed)
  k_NR3C1_IL10  = 1.5,      # glucocorticoid-driven production
  IL10_max      = 1.0,      # uM

  # CXCL8 (degree = 33)
  k_prod_CXCL8  = 0.2,      # uM/hr
  k_deg_CXCL8   = 0.231,    # 1/hr (t1/2 ~3 hr)
  CXCL8_basal   = 0.005,    # uM
  CXCL8_FMQ     = 0.04,     # uM
  k_NFkB_CXCL8  = 2.0,      # NF-kB driven production
  CXCL8_max     = 1.0,      # uM

  # CCL2 (monocyte chemotaxis)
  k_prod_CCL2   = 0.15,     # uM/hr
  k_deg_CCL2    = 0.116,    # 1/hr
  CCL2_basal    = 0.005,    # uM
  CCL2_FMQ      = 0.03,     # uM
  k_NFkB_CCL2   = 1.5       # NF-kB driven production
)

# ============================================================
# PAIN SIGNALING PARAMETERS
# ============================================================

pain_params <- list(
  # BDNF (degree = 41, neurotrophin)
  k_prod_BDNF   = 0.1,      # uM/hr
  k_deg_BDNF    = 0.023,    # 1/hr (t1/2 ~30 hr)
  BDNF_basal    = 0.05,     # uM, healthy baseline
  BDNF_FMQ      = 0.15,     # uM, fibromyalgia baseline (elevated)
  BDNF_max      = 1.0,      # uM
  k_IL6_BDNF    = 0.5,      # IL6-mediated BDNF upregulation

  # NGF (degree = 34)
  k_prod_NGF    = 0.08,     # uM/hr
  k_deg_NGF     = 0.023,    # 1/hr
  NGF_basal     = 0.03,     # uM
  NGF_FMQ       = 0.1,      # uM
  NGF_max       = 1.0,      # uM
  k_TNF_NGF     = 0.3,      # TNF-mediated NGF upregulation

  # TAC1 / Substance P (degree = 38)
  k_prod_TAC1   = 0.15,     # uM/hr
  k_deg_TAC1    = 0.693,    # 1/hr (t1/2 ~1 hr)
  TAC1_basal    = 0.02,     # uM
  TAC1_FMQ      = 0.08,     # uM
  TAC1_max      = 1.0,      # uM
  k_BDNF_TAC1   = 1.0,      # BDNF-driven TAC1 release

  # TRPV1 sensitization
  TRPV1_basal   = 0.1,      # basal activity (normalized)
  TRPV1_max     = 1.0,      # max activity
  EC50_NGF_TRPV1= 0.05,     # uM, EC50 for NGF -> TRPV1 sensitization
  n_NGF_TRPV1   = 2.0,      # Hill coefficient
  k_BDNF_TRPV1  = 0.5,      # BDNF contribution to TRPV1

  # Pain integration (sigmoidal)
  Pain_basal    = 5.0,      # VAS baseline (0-10 scale) in healthy
  Pain_FMQ      = 7.0,      # VAS baseline in fibromyalgia
  Emax_BDNF     = 0.5,      # max BDNF effect on pain
  EC50_BDNF     = 0.1,      # uM
  Emax_TNF      = 0.3,      # max TNF effect on pain
  EC50_TNF      = 0.05,     # uM
  Emax_TAC1     = 0.4,      # max TAC1 effect on pain
  EC50_TAC1     = 0.05,     # uM
  E_IL10        = 0.5,      # IL10 protective effect
  EC50_IL10     = 0.02,     # uM
  n_pain        = 1.5       # Hill coefficient for pain integration
)

# ============================================================
# HPA AXIS PARAMETERS
# ============================================================

hpa_params <- list(
  # CRH (corticotropin-releasing hormone)
  k_prod_CRH    = 0.1,      # uM/hr
  k_deg_CRH     = 0.693,    # 1/hr
  CRH_basal     = 0.02,     # uM
  CRH_FMQ       = 0.05,     # uM (elevated in chronic stress)
  k_TNF_CRH     = 0.2,      # TNF-driven CRH upregulation

  # POMC -> ACTH -> Cortisol
  k_prod_POMC   = 0.1,      # uM/hr
  k_deg_POMC    = 0.347,    # 1/hr
  POMC_basal    = 0.02,     # uM
  POMC_FMQ      = 0.04,     # uM
  k_CRH_POMC    = 2.0,      # CRH -> POMC stimulation

  # NR3C1 (glucocorticoid receptor)
  NR3C1_total   = 1.0,      # uM, total receptor (normalized)
  NR3C1_basal   = 0.3,      # basal activity (fraction)
  EC50_cortisol_NR3C1 = 0.1, # uM
  n_NR3C1       = 1.5,      # Hill coefficient
  k_TNF_NR3C1_resist = 0.5, # TNF-induced glucocorticoid resistance

  # NPY (neuropeptide Y, anxiolytic)
  k_prod_NPY    = 0.05,     # uM/hr
  k_deg_NPY     = 0.116,    # 1/hr
  NPY_basal     = 0.03,     # uM
  NPY_FMQ       = 0.015,    # uM (suppressed in FMQ)
  k_CRH_NPY_inhib = 0.3    # CRH-mediated NPY suppression
)

# ============================================================
# MONOAMINE NEUROTRANSMISSION PARAMETERS
# ============================================================

monoamine_params <- list(
  # SLC6A4 (serotonin transporter)
  SLC6A4_Vmax   = 1.0,      # uM/hr, max transport rate
  SLC6A4_Km     = 0.5,      # uM, Km for serotonin uptake
  SLC6A4_expression = 1.0,  # relative expression (normalized)

  # COMT (catechol-O-methyltransferase)
  COMT_Vmax     = 0.8,      # uM/hr
  COMT_Km       = 2.0,      # uM
  COMT_expression = 1.0,    # relative expression

  # OPRM1 (mu-opioid receptor) - KEY DOSE-DEPENDENT TARGET
  # Literature: Naltrexone Kd for OPRM1 = 0.1-1 nM (0.0001-0.001 uM)
  # This is ~10,000x tighter than TLR4 binding
  OPRM1_total    = 1.0,     # uM (normalized)
  OPRM1_basal    = 0.2,     # basal endogenous opioid tone
  Kd_NTX_OPRM1   = 0.001,   # uM (= 1 nM, very high affinity from literature)
  Kd_endorphin   = 0.005,   # uM, endogenous beta-endorphin Kd for OPRM1
  Endorphin_conc = 0.002,   # uM, basal endogenous beta-endorphin in CNS
  OPRM1_agonism_max = 0.8,  # max analgesic effect from OPRM1 activation
  EC50_OPRM1_pain  = 0.3,   # EC50 for OPRM1 analgesia
  n_OPRM1_pain     = 1.5,   # Hill coefficient
  # Paradoxical hyperalgesia when OPRM1 is blocked at high doses
  k_hyperalgesia   = 0.5,   # pain amplification when endorphin signaling blocked
  # OPRM1 desensitization (not relevant at LDN, relevant at high dose)
  k_beta_arrestin  = 0.5,   # desensitization rate

  # Serotonin dynamics
  HT_basal      = 0.1,      # uM, basal serotonin
  HT_FMQ        = 0.06,     # uM, depleted in FMQ
  k_syn_HT      = 0.1,      # uM/hr, synthesis rate
  k_deg_HT      = 0.5,      # 1/hr, MAO-mediated degradation
  k_COMT_HT     = 0.2       # COMT contribution
)

# ============================================================
# CROSS-TALK / FEEDBACK PARAMETERS
# ============================================================

crosstalk_params <- list(
  # Inflammatory -> Pain amplification
  k_TNF_TRPV1_sens  = 0.3,  # TNF -> TRPV1 sensitization
  k_IL1B_TAC1_rel    = 0.2,  # IL1B -> TAC1 release
  k_IL6_BDNF         = 0.3,  # IL6 -> BDNF upregulation

  # Neurotransmitter -> Inflammation
  k_HT_IL10          = 0.2,  # serotonin -> IL10 production
  k_NE_TNF_inhib     = 0.1,  # norepinephrine -> TNF suppression

  # HPA -> Immune
  k_cortisol_immune  = 0.5,  # cortisol -> broad anti-inflammatory
  k_NR3C1_TNF_inhib  = 0.3,  # NR3C1 -> TNF transcription inhibition

  # Metabolic-inflammatory link
  k_leptin_TNF       = 0.1,  # leptin -> TNF production
  k_leptin_IL6       = 0.1   # leptin -> IL6 production
)

# ============================================================
# COMPILE ALL PARAMETERS
# ============================================================

all_params <- c(
  pk_params,
  tlr4_params,
  cytokine_params,
  pain_params,
  hpa_params,
  monoamine_params,
  crosstalk_params
)
