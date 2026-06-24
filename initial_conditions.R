# Initial Conditions: Healthy vs Fibromyalgia Baseline States
# All concentrations in uM (normalized where indicated)

# ============================================================
# HEALTHY BASELINE STATE
# ============================================================

healthy_state <- list(
  # PK compartments
  A_gut       = 0,
  C_plasma    = 0,
  C_periph    = 0,
  C_CNS       = 0,
  A_metab     = 0,      # active metabolite (6-beta-naltrexol)

  # TLR4 pathway
  TLR4_free   = 1.0,    # fully available (normalized)
  TLR4_bound  = 0.0,
  NFkB        = 0.05,   # basal activity
  NLRP3       = 0.1,    # basal priming

  # Cytokines
  TNF         = 0.01,   # uM
  IL1B        = 0.005,  # uM
  IL6         = 0.01,   # uM
  IL10        = 0.02,   # uM
  CXCL8       = 0.005,  # uM
  CCL2        = 0.005,  # uM

  # Pain signaling
  BDNF        = 0.05,   # uM
  NGF         = 0.03,   # uM
  TAC1        = 0.02,   # uM (Substance P)
  TRPV1       = 0.1,    # normalized activity

  # HPA axis
  CRH         = 0.02,   # uM
  POMC        = 0.02,   # uM
  Cortisol    = 0.1,    # uM (morning level)
  NR3C1_act   = 0.3,    # fractional activity
  NPY         = 0.03,   # uM
  LEP         = 0.05,   # uM

  # Monoamines
  HT          = 0.1,    # uM (serotonin)
  DA          = 0.05,   # uM (dopamine)
  NE          = 0.03,   # uM (norepinephrine)
  OPRM1_act   = 0.2,    # fractional activity

  # Clinical outputs
  Pain_VAS    = 2.0,    # 0-10 scale
  FIQ_score   = 15.0,   # Fibromyalgia Impact Questionnaire
  Fatigue     = 2.0     # 0-10 scale
)

# ============================================================
# FIBROMYALGIA DISEASE STATE
# ============================================================

fibromyalgia_state <- list(
  # PK compartments (no drug initially)
  A_gut       = 0,
  C_plasma    = 0,
  C_periph    = 0,
  C_CNS       = 0,
  A_metab     = 0,

  # TLR4 pathway (microglial activation)
  TLR4_free   = 1.0,    # fully available
  TLR4_bound  = 0.0,
  NFkB        = 0.35,   # elevated basal activity (neuroinflammation)
  NLRP3       = 0.5,    # elevated priming

  # Cytokines (elevated pro-inflammatory, suppressed anti-inflammatory)
  TNF         = 0.08,   # 8x healthy (from network: degree=45)
  IL1B        = 0.05,   # 10x healthy (from network: degree=42)
  IL6         = 0.10,   # 10x healthy (from network: degree=45)
  IL10        = 0.015,  # 0.75x healthy (suppressed)
  CXCL8       = 0.04,   # 8x healthy (from network: degree=33)
  CCL2        = 0.03,   # 6x healthy

  # Pain signaling (sensitized)
  BDNF        = 0.15,   # 3x healthy (from network: degree=41)
  NGF         = 0.10,   # 3.3x healthy (from network: degree=34)
  TAC1        = 0.08,   # 4x healthy (from network: degree=38)
  TRPV1       = 0.5,    # 5x healthy (sensitized)

  # HPA axis (dysregulated)
  CRH         = 0.05,   # 2.5x healthy (chronic stress)
  POMC        = 0.04,   # 2x healthy
  Cortisol    = 0.15,   # elevated but blunted response
  NR3C1_act   = 0.2,    # reduced (glucocorticoid resistance from TNF)
  NPY         = 0.015,  # 0.5x healthy (depleted)
  LEP         = 0.08,   # elevated (leptin resistance)

  # Monoamines (depleted)
  HT          = 0.06,   # 0.6x healthy (serotonin deficit)
  DA          = 0.03,   # 0.6x healthy (dopamine deficit)
  NE          = 0.02,   # 0.67x healthy (norepinephrine deficit)
  OPRM1_act   = 0.15,   # reduced opioid tone

  # Clinical outputs (elevated)
  Pain_VAS    = 7.0,    # moderate-severe pain
  FIQ_score   = 65.0,   # high impact
  Fatigue     = 7.0     # severe fatigue
)

# ============================================================
# LDN-TREATED STATE (Expected after 4 weeks at 4.5 mg QD)
# ============================================================

ldn_treated_state <- list(
  A_gut       = 0,
  C_plasma    = 0,
  C_periph    = 0,
  C_CNS       = 0.05,   # steady-state CNS concentration (low)
  A_metab     = 0,

  # TLR4 pathway (partially inhibited)
  TLR4_free   = 0.7,    # 30% occupied by naltrexone
  TLR4_bound  = 0.3,
  NFkB        = 0.15,   # reduced from 0.35
  NLRP3       = 0.25,   # reduced from 0.5

  # Cytokines (normalized toward healthy)
  TNF         = 0.03,   # reduced from 0.08
  IL1B        = 0.015,  # reduced from 0.05
  IL6         = 0.04,   # reduced from 0.10
  IL10        = 0.025,  # increased from 0.015
  CXCL8       = 0.015,  # reduced from 0.04
  CCL2        = 0.012,  # reduced from 0.03

  # Pain signaling (desensitized)
  BDNF        = 0.08,   # reduced from 0.15
  NGF         = 0.05,   # reduced from 0.10
  TAC1        = 0.035,  # reduced from 0.08
  TRPV1       = 0.2,    # reduced from 0.5

  # HPA axis (partially restored)
  CRH         = 0.03,   # reduced from 0.05
  POMC        = 0.025,  # reduced from 0.04
  Cortisol    = 0.12,   # normalized
  NR3C1_act   = 0.28,   # improved from 0.2
  NPY         = 0.025,  # improved from 0.015
  LEP         = 0.06,   # improved from 0.08

  # Monoamines (partially restored)
  HT          = 0.08,   # improved from 0.06
  DA          = 0.04,   # improved from 0.03
  NE          = 0.025,  # improved from 0.02
  OPRM1_act   = 0.18,   # improved from 0.15

  # Clinical outputs (improved)
  Pain_VAS    = 4.5,    # moderate pain (36% reduction)
  FIQ_score   = 42.0,   # moderate impact (35% reduction)
  Fatigue     = 4.5     # moderate fatigue (36% reduction)
)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

get_state_vector <- function(state_list) {
  unlist(state_list)
}

get_disease_state <- function(condition = "fibromyalgia") {
  switch(condition,
    healthy      = healthy_state,
    fibromyalgia = fibromyalgia_state,
    ldn_treated  = ldn_treated_state,
    stop("Unknown condition: ", condition)
  )
}

compare_states <- function(state1, state2) {
  s1 <- get_state_vector(state1)
  s2 <- get_state_vector(state2)
  data.frame(
    Variable = names(s1),
    Healthy  = s1,
    Disease  = s2,
    Fold_Change = s2 / s1,
    Percent_Change = ((s2 - s1) / s1) * 100
  )
}
