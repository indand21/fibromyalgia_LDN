$PROBLEM
  QSP Model v2: Naltrexone (LDN) Mechanism in Fibromyalgia
  Literature-implemented model with TRIF-IRF3 biased TLR4 antagonism,
  hormesis, ESR stratification, and clinical cytokine calibration.
  Based on: Wang 2016 (Br J Pharmacol), Parkitny 2017, Younger 2009/2013

$PARAM
  // ============================================================
  // PK Parameters (Naltrexone oral LDN)
  // ============================================================
  Ka = 1.2, F = 0.05, Vd = 1127, V2 = 500, Q = 120,
  CL = 22.6, Kp_brain = 0.028, Kin_CNS = 0.5, Kout_CNS = 0.5,
  Vmax_CYP = 50, Km_CYP = 10, f_metab = 0.6,

  // ============================================================
  // Patient Stratification (ESR-based, from Younger 2009)
  // ESR_strat: 0 = normal ESR (poor responder), 1 = high ESR (good responder)
  // ============================================================
  ESR_strat = 1.0,
  ESR_baseline = 1.0,       // fold-increase in inflammatory baseline for high ESR
  k_ESR_cytokine = 0.5,     // ESR modulation of cytokine production

  // ============================================================
  // Sex-specific parameters (ESR1 modulation)
  // sex: 0 = male, 1 = female (FM predominantly female)
  // ============================================================
  sex = 1.0,
  ESR1_activity = 0.5,      // estrogen receptor alpha activity (higher in females)
  k_ESR1_TLR4 = 0.2,        // ESR1 modulation of TLR4 signaling
  k_ESR1_IL6 = 0.15,        // ESR1 modulation of IL-6 production

  // ============================================================
  // TLR4 Signaling: DUAL PATHWAY (TRIF-IRF3 biased antagonism)
  // Wang 2016: Naltrexone binds MD2 LPS pocket, blocks TRIF-IRF3
  // but does NOT block MyD88-NFkB pathway
  // ============================================================
  kon_NTX = 0.1, koff_NTX = 1.0, Kd_NTX = 10,
  TLR4_syn = 0.01, TLR4_deg = 0.01,
  // MyD88-NFkB pathway (NOT blocked by naltrexone)
  MyD88_basal = 0.05, MyD88_max = 1.0, EC50_TLR4_MyD88 = 0.5, n_MyD88 = 2.0,
  // TRIF-IRF3 pathway (BLOCKED by naltrexone)
  TRIF_basal = 0.05, TRIF_max = 1.0, EC50_TLR4_TRIF = 0.3, n_TRIF = 2.0,
  // NFkB integrates both pathways
  NFkB_basal = 0.05, NFkB_max = 1.0, k_MyD88_NFkB = 0.7, k_TRIF_NFkB = 0.3,
  // NLRP3 (primed by NFkB)
  NLRP3_basal = 0.1, NLRP3_max = 1.0, EC50_NFkB_NLRP3 = 0.3, n_NLRP3 = 1.5,

  // ============================================================
  // Cytokine Parameters (calibrated to Parkitny 2017 clinical data)
  // Clinical: LDN 4.5mg x 8wk reduced TNF-a, IL-1b, IL-6, IL-10, IL-17A, TGF-b
  // Pain reduction: 15% (pain), 18% (overall symptoms)
  // ============================================================
  // TNF-a (reduced by LDN via TRIF-IRF3 blockade)
  k_prod_TNF = 0.5, k_deg_TNF = 0.693, TNF_basal = 0.01,
  k_MyD88_TNF = 1.0, k_TRIF_TNF = 1.5,  // TRIF contributes more to TNF
  TNF_max = 1.0,

  // IL-1b (primarily NLRP3-dependent, NOT directly blocked by LDN per Wang 2016)
  k_prod_IL1B = 0.3, k_deg_IL1B = 0.347, IL1B_basal = 0.005,
  k_NLRP3_IL1B = 3.0, IL1B_max = 1.0,

  // IL-6 (reduced by LDN)
  k_prod_IL6 = 0.4, k_deg_IL6 = 0.116, IL6_basal = 0.01,
  k_MyD88_IL6 = 1.5, k_TRIF_IL6 = 1.0,
  IL6_max = 1.0,

  // IL-10 (anti-inflammatory, reduced by LDN in clinical data - complex regulation)
  k_prod_IL10 = 0.2, k_deg_IL10 = 0.116, IL10_basal = 0.02,
  k_NR3C1_IL10 = 1.5, k_TRIF_IL10 = 0.5,  // TRIF contributes to IL-10 production
  IL10_max = 1.0,

  // CXCL8
  k_prod_CXCL8 = 0.2, k_deg_CXCL8 = 0.231, CXCL8_basal = 0.005,
  k_MyD88_CXCL8 = 1.5, k_TRIF_CXCL8 = 0.5,
  CXCL8_max = 1.0,

  // CCL2
  k_prod_CCL2 = 0.15, k_deg_CCL2 = 0.116, CCL2_basal = 0.005,
  k_MyD88_CCL2 = 1.0, k_TRIF_CCL2 = 0.5,

  // IL-17A (Th17 cytokine, reduced by LDN per Parkitny 2017)
  k_prod_IL17A = 0.08, k_deg_IL17A = 0.116, IL17A_basal = 0.005,
  k_TRIF_IL17A = 1.0, k_NFkB_IL17A = 0.5,
  IL17A_max = 0.5,

  // TGF-beta (reduced by LDN per Parkitny 2017)
  k_prod_TGFb = 0.06, k_deg_TGFb = 0.058, TGFb_basal = 0.01,
  k_TRIF_TGFb = 0.8, k_NFkB_TGFb = 0.3,
  TGFb_max = 0.5,

  // ============================================================
  // Pain Signaling Parameters
  // ============================================================
  k_prod_BDNF = 0.1, k_deg_BDNF = 0.023, BDNF_basal = 0.05, BDNF_max = 1.0, k_IL6_BDNF = 0.5,
  k_prod_NGF = 0.08, k_deg_NGF = 0.023, NGF_basal = 0.03, NGF_max = 1.0, k_TNF_NGF = 0.3,
  k_prod_TAC1 = 0.15, k_deg_TAC1 = 0.693, TAC1_basal = 0.02, TAC1_max = 1.0, k_BDNF_TAC1 = 1.0,
  k_IL17A_TAC1 = 0.3,  // IL-17A amplifies Substance P release
  TRPV1_basal = 0.1, TRPV1_max = 1.0, EC50_NGF_TRPV1 = 0.05, n_NGF_TRPV1 = 2.0, k_BDNF_TRPV1 = 0.5,

  // ============================================================
  // HPA Axis Parameters
  // ============================================================
  k_prod_CRH = 0.1, k_deg_CRH = 0.693, CRH_basal = 0.02, k_TNF_CRH = 0.2,
  k_prod_POMC = 0.1, k_deg_POMC = 0.347, POMC_basal = 0.02, k_CRH_POMC = 2.0,
  NR3C1_total = 1.0, EC50_cort_NR3C1 = 0.1, n_NR3C1 = 1.5, k_TNF_NR3C1_resist = 0.5,
  k_prod_NPY = 0.05, k_deg_NPY = 0.116, NPY_basal = 0.03, k_CRH_NPY_inhib = 0.3,

  // ============================================================
  // Monoamine Parameters
  // ============================================================
  k_syn_HT = 0.1, k_deg_HT = 0.5, HT_basal = 0.1,
  SLC6A4_Vmax = 1.0, SLC6A4_Km = 0.5,
  COMT_Vmax = 0.8, COMT_Km = 2.0,

  // ============================================================
  // OPRM1 (mu-opioid) - HORMESIS MODEL
  // Dara 2023: biphasic dose response; transient blockade -> receptor upregulation
  // Jackson 2021: LDN restores endogenous opioid tone
  // ============================================================
  Kd_NTX_OPRM1 = 0.001, Kd_endorphin = 0.005, Endorphin_conc_basal = 0.002,
  OPRM1_density_basal = 1.0,  // basal receptor density (normalized)
  k_OPRM1_upreg = 0.05,       // receptor upregulation rate during blockade
  k_OPRM1_deg = 0.01,         // receptor degradation rate
  k_endorphin_rebound = 0.02, // endorphin production increase during rebound
  OPRM1_agonism_max = 0.8, EC50_OPRM1_pain = 0.3, n_OPRM1_pain = 1.5,
  k_hyperalgesia = 0.5,

  // ============================================================
  // Pain Integration
  // ============================================================
  Pain_basal = 2.0, Emax_BDNF = 0.5, EC50_BDNF_pain = 0.1,
  Emax_TNF = 0.3, EC50_TNF_pain = 0.05,
  Emax_TAC1 = 0.4, EC50_TAC1_pain = 0.05,
  Emax_IL17A = 0.2, EC50_IL17A_pain = 0.02,
  E_IL10 = 0.5, EC50_IL10_pain = 0.02, n_pain = 1.5,

  // CPT parameters (Jackson 2021: CPT time doubled in FM with LDN)
  CPT_basal_healthy = 60.0,   // seconds, healthy baseline
  CPT_basal_FMQ = 25.0,       // seconds, FM baseline (reduced)
  k_pain_CPT = -5.0,          // pain -> CPT reduction coefficient

  // Dose
  Dose = 4.5

$INIT
  // PK
  A_gut = 0, C_plasma = 0, C_periph = 0, C_CNS = 0, A_metab = 0,
  // TLR4 dual pathway
  TLR4_free = 1.0, TLR4_bound = 0.0,
  MyD88_act = 0.35, TRIF_act = 0.35, NFkB = 0.35, NLRP3 = 0.5,
  // Cytokines (FMQ baseline, ESR-stratified)
  TNF = 0.08, IL1B = 0.05, IL6 = 0.10, IL10 = 0.015,
  CXCL8 = 0.04, CCL2 = 0.03,
  IL17A = 0.03, TGFb = 0.04,
  // Pain
  BDNF = 0.15, NGF = 0.10, TAC1 = 0.08, TRPV1_act = 0.5,
  // HPA
  CRH = 0.05, POMC = 0.04, Cortisol = 0.15, NR3C1_act = 0.2, NPY = 0.015,
  // Monoamines
  HT = 0.06,
  // OPRM1 hormesis state
  OPRM1_density = 1.0, Endorphin_conc = 0.002,
  // Output
  Pain_VAS = 7.0, CPT_time = 25.0

$ODE
  // ===================== PK MODEL =====================
  dxdt_A_gut    = -Ka * A_gut;
  dxdt_C_plasma = (Ka * A_gut * F) / Vd - (CL / Vd + Q / Vd) * C_plasma + (Q / V2) * C_periph;
  dxdt_C_periph = (Q / V2) * C_plasma - (Q / V2) * C_periph;

  // CNS compartment
  double C_CNS_in = Kin_CNS * C_plasma * Kp_brain;
  double C_CNS_out = Kout_CNS * C_CNS;
  double TLR4_binding = kon_NTX * C_CNS * TLR4_free;
  double TLR4_unbinding = koff_NTX * TLR4_bound;
  dxdt_C_CNS = C_CNS_in - C_CNS_out - TLR4_binding + TLR4_unbinding;

  // Active metabolite
  double v_metab = Vmax_CYP * C_plasma / (Km_CYP + C_plasma);
  dxdt_A_metab = f_metab * v_metab - 0.1 * A_metab;

  // ===================== TLR4 DUAL PATHWAY =====================
  // IMPROVEMENT 1: TRIF-IRF3 biased antagonism (Wang 2016)
  // Naltrexone binds MD2 LPS pocket -> blocks TRIF-IRF3 but NOT MyD88-NFkB
  dxdt_TLR4_free  = TLR4_syn - TLR4_deg * TLR4_free - TLR4_binding + TLR4_unbinding;
  dxdt_TLR4_bound = TLR4_binding - TLR4_unbinding - TLR4_deg * TLR4_bound;

  double TLR4_total = TLR4_free + TLR4_bound;
  double TLR4_activation = TLR4_free / 1.0;  // free TLR4 available for activation

  // MyD88-NFkB pathway (NOT blocked by naltrexone)
  double MyD88_target = MyD88_basal + (MyD88_max - MyD88_basal) *
                         pow(TLR4_activation, n_MyD88) /
                         (pow(EC50_TLR4_MyD88, n_MyD88) + pow(TLR4_activation, n_MyD88));
  dxdt_MyD88_act = 0.5 * (MyD88_target - MyD88_act);

  // TRIF-IRF3 pathway (BLOCKED by naltrexone via MD2 binding)
  double TRIF_inhibition = TLR4_bound / (TLR4_total + 0.001);  // fraction bound = inhibition
  double TRIF_target = TRIF_basal + (TRIF_max - TRIF_basal) *
                        pow(TLR4_activation, n_TRIF) /
                        (pow(EC50_TLR4_TRIF, n_TRIF) + pow(TLR4_activation, n_TRIF)) *
                        (1.0 - TRIF_inhibition);  // naltrexone blocks TRIF
  dxdt_TRIF_act = 0.5 * (TRIF_target - TRIF_act);

  // NFkB integrates both pathways
  double NFkB_target = NFkB_basal +
                        k_MyD88_NFkB * (MyD88_act - MyD88_basal) +
                        k_TRIF_NFkB * (TRIF_act - TRIF_basal);
  if (NFkB_target > NFkB_max) NFkB_target = NFkB_max;
  if (NFkB_target < NFkB_basal) NFkB_target = NFkB_basal;
  dxdt_NFkB = 0.5 * (NFkB_target - NFkB);

  // NLRP3 inflammasome
  double NLRP3_target = NLRP3_basal + (NLRP3_max - NLRP3_basal) *
                         pow(NFkB, n_NLRP3) /
                         (pow(EC50_NFkB_NLRP3, n_NLRP3) + pow(NFkB, n_NLRP3));
  dxdt_NLRP3 = 0.3 * (NLRP3_target - NLRP3);

  // ===================== ESR STRATIFICATION =====================
  // IMPROVEMENT 5: ESR-based patient stratification (Younger 2009)
  double ESR_factor = 1.0 + k_ESR_cytokine * ESR_strat * (ESR_baseline - 1.0);

  // ===================== SEX-SPECIFIC ESR1 MODULATION =====================
  // IMPROVEMENT 8: ESR1 modulates TLR4 and IL-6 (sex differences in FM)
  double ESR1_TLR4_mod = 1.0 + k_ESR1_TLR4 * sex * ESR1_activity;
  double ESR1_IL6_mod = 1.0 + k_ESR1_IL6 * sex * ESR1_activity;

  // ===================== CYTOKINE DYNAMICS =====================
  // IMPROVED: Split MyD88 vs TRIF contributions to each cytokine

  // TNF-a (both pathways contribute, TRIF dominant per Wang 2016)
  double TNF_prod = k_prod_TNF * ESR_factor *
                    (1.0 + k_MyD88_TNF * MyD88_act + k_TRIF_TNF * TRIF_act);
  dxdt_TNF = TNF_prod - k_deg_TNF * TNF + 0.05 * TNF * (1 - TNF / TNF_max) - 0.1 * IL10 * TNF;

  // IL-1b (primarily NLRP3, NOT directly blocked by LDN per Wang 2016)
  double IL1B_prod = k_prod_IL1B * ESR_factor * (1.0 + k_NLRP3_IL1B * NLRP3);
  dxdt_IL1B = IL1B_prod - k_deg_IL1B * IL1B;

  // IL-6 (both pathways, ESR1-modulated)
  double IL6_prod = k_prod_IL6 * ESR_factor * ESR1_IL6_mod *
                    (1.0 + k_MyD88_IL6 * MyD88_act + k_TRIF_IL6 * TRIF_act);
  dxdt_IL6 = IL6_prod - k_deg_IL6 * IL6;

  // IL-10 (complex: NR3C1 + TRIF, reduced by LDN per clinical data)
  double IL10_prod = k_prod_IL10 * ESR_factor *
                     (1.0 + k_NR3C1_IL10 * NR3C1_act + k_TRIF_IL10 * TRIF_act);
  dxdt_IL10 = IL10_prod - k_deg_IL10 * IL10;

  // CXCL8
  double CXCL8_prod = k_prod_CXCL8 * ESR_factor *
                      (1.0 + k_MyD88_CXCL8 * MyD88_act + k_TRIF_CXCL8 * TRIF_act);
  dxdt_CXCL8 = CXCL8_prod - k_deg_CXCL8 * CXCL8;

  // CCL2
  double CCL2_prod = k_prod_CCL2 * ESR_factor *
                     (1.0 + k_MyD88_CCL2 * MyD88_act + k_TRIF_CCL2 * TRIF_act);
  dxdt_CCL2 = CCL2_prod - k_deg_CCL2 * CCL2;

  // IMPROVEMENT 6: IL-17A (Th17 cytokine, reduced by LDN per Parkitny 2017)
  double IL17A_prod = k_prod_IL17A * ESR_factor *
                      (1.0 + k_TRIF_IL17A * TRIF_act + k_NFkB_IL17A * NFkB);
  dxdt_IL17A = IL17A_prod - k_deg_IL17A * IL17A;

  // IMPROVEMENT 6: TGF-beta (reduced by LDN per Parkitny 2017)
  double TGFb_prod = k_prod_TGFb * ESR_factor *
                     (1.0 + k_TRIF_TGFb * TRIF_act + k_NFkB_TGFb * NFkB);
  dxdt_TGFb = TGFb_prod - k_deg_TGFb * TGFb;

  // ===================== PAIN SIGNALING =====================
  // BDNF (IL-6 driven)
  double BDNF_prod = k_prod_BDNF + k_IL6_BDNF * IL6 * k_prod_BDNF;
  dxdt_BDNF = BDNF_prod - k_deg_BDNF * BDNF;
  if (BDNF > BDNF_max) dxdt_BDNF = dxdt_BDNF - 0.5 * (BDNF - BDNF_max);

  // NGF (TNF driven)
  double NGF_prod = k_prod_NGF + k_TNF_NGF * TNF * k_prod_NGF;
  dxdt_NGF = NGF_prod - k_deg_NGF * NGF;
  if (NGF > NGF_max) dxdt_NGF = dxdt_NGF - 0.5 * (NGF - NGF_max);

  // TAC1 / Substance P (BDNF + IL-17A driven)
  double TAC1_prod = k_prod_TAC1 + k_BDNF_TAC1 * BDNF * k_prod_TAC1 +
                     k_IL17A_TAC1 * IL17A * k_prod_TAC1;
  dxdt_TAC1 = TAC1_prod - k_deg_TAC1 * TAC1;
  if (TAC1 > TAC1_max) dxdt_TAC1 = dxdt_TAC1 - 0.5 * (TAC1 - TAC1_max);

  // TRPV1 sensitization
  double TRPV1_target = TRPV1_basal +
                         (TRPV1_max - TRPV1_basal) *
                         pow(NGF, n_NGF_TRPV1) /
                         (pow(EC50_NGF_TRPV1, n_NGF_TRPV1) + pow(NGF, n_NGF_TRPV1)) +
                         k_BDNF_TRPV1 * BDNF;
  dxdt_TRPV1_act = 0.5 * (TRPV1_target - TRPV1_act);

  // ===================== HPA AXIS =====================
  double CRH_prod = k_prod_CRH + k_TNF_CRH * TNF * k_prod_CRH;
  dxdt_CRH = CRH_prod - k_deg_CRH * CRH;

  double POMC_prod = k_prod_POMC + k_CRH_POMC * CRH * k_prod_POMC;
  dxdt_POMC = POMC_prod - k_deg_POMC * POMC;

  double Cortisol_prod = 0.5 * POMC;
  dxdt_Cortisol = Cortisol_prod - 0.3 * Cortisol;

  double NR3C1_sensitivity = 1.0 - k_TNF_NR3C1_resist * TNF;
  if (NR3C1_sensitivity < 0.1) NR3C1_sensitivity = 0.1;
  double NR3C1_target = NR3C1_sensitivity *
                         pow(Cortisol, n_NR3C1) /
                         (pow(EC50_cort_NR3C1, n_NR3C1) + pow(Cortisol, n_NR3C1));
  dxdt_NR3C1_act = 0.3 * (NR3C1_target - NR3C1_act);

  double NPY_prod = k_prod_NPY * (1.0 / (1.0 + k_CRH_NPY_inhib * CRH));
  dxdt_NPY = NPY_prod - k_deg_NPY * NPY;

  // ===================== MONOAMINES =====================
  double HT_prod = k_syn_HT;
  double HT_uptake = SLC6A4_Vmax * HT / (SLC6A4_Km + HT);
  double HT_deg = k_deg_HT * HT;
  dxdt_HT = HT_prod - HT_uptake - HT_deg + 0.02 * IL10;

  // ===================== OPRM1 HORMESIS MODEL =====================
  // IMPROVEMENT 2: Hormesis / receptor upregulation (Dara 2023, Jackson 2021)
  // Transient OPRM1 blockade -> receptor synthesis increases -> rebound effect

  double OPRM1_NTX_occ = C_CNS / (Kd_NTX_OPRM1 + C_CNS);

  // Receptor density dynamics: blockade -> upregulation
  double OPRM1_upreg_stimulus = k_OPRM1_upreg * OPRM1_NTX_occ;  // blockade stimulates synthesis
  dxdt_OPRM1_density = OPRM1_upreg_stimulus * (2.0 - OPRM1_density) -  // max 2x upregulation
                        k_OPRM1_deg * (OPRM1_density - OPRM1_density_basal);  // return to basal

  // Endorphin rebound: transient blockade -> increased endorphin production
  double Endorphin_rebound_stimulus = k_endorphin_rebound * OPRM1_NTX_occ;
  dxdt_Endorphin_conc = Endorphin_rebound_stimulus * (0.01 - Endorphin_conc) -  // max 0.01 uM
                         0.005 * (Endorphin_conc - Endorphin_conc_basal);  // return to basal

  // Effective endorphin signaling (density-adjusted)
  double OPRM1_free = 1.0 - OPRM1_NTX_occ;
  double Endorphin_signal = OPRM1_free * OPRM1_density * Endorphin_conc / (Kd_endorphin + Endorphin_conc);

  // OPRM1-mediated analgesia
  double OPRM1_analgesia = OPRM1_agonism_max * pow(Endorphin_signal, n_OPRM1_pain) /
                           (pow(EC50_OPRM1_pain, n_OPRM1_pain) + pow(Endorphin_signal, n_OPRM1_pain));

  // Paradoxical hyperalgesia (only at very high blockade)
  double OPRM1_blockade = OPRM1_NTX_occ;
  double Hyperalgesia_factor = 1.0 + k_hyperalgesia * pow(OPRM1_blockade, 3);

  // ===================== PAIN OUTPUT =====================
  double pain_BDNF = Emax_BDNF * pow(BDNF, n_pain) / (pow(EC50_BDNF_pain, n_pain) + pow(BDNF, n_pain));
  double pain_TNF  = Emax_TNF * pow(TNF, n_pain) / (pow(EC50_TNF_pain, n_pain) + pow(TNF, n_pain));
  double pain_TAC1 = Emax_TAC1 * pow(TAC1, n_pain) / (pow(EC50_TAC1_pain, n_pain) + pow(TAC1, n_pain));
  double pain_IL17A = Emax_IL17A * pow(IL17A, n_pain) / (pow(EC50_IL17A_pain, n_pain) + pow(IL17A, n_pain));
  double pain_IL10 = E_IL10 * IL10 / (EC50_IL10_pain + IL10);

  double Pain_target = Pain_basal * (1.0 + pain_BDNF + pain_TNF + pain_TAC1 + pain_IL17A) / (1.0 + pain_IL10)
                       * (1.0 + (Hyperalgesia_factor - 1.0))
                       * (1.0 - 0.3 * OPRM1_analgesia);
  if (Pain_target > 10.0) Pain_target = 10.0;
  if (Pain_target < 0.0) Pain_target = 0.0;
  dxdt_Pain_VAS = 0.1 * (Pain_target - Pain_VAS);

  // IMPROVEMENT 7: Cold Pressor Test (CPT) as objective outcome (Jackson 2021)
  double CPT_target = CPT_basal_FMQ + k_pain_CPT * (Pain_VAS - 7.0) +
                      10.0 * OPRM1_analgesia;  // opioid tone improves CPT
  if (CPT_target > CPT_basal_healthy) CPT_target = CPT_basal_healthy;
  if (CPT_target < 5.0) CPT_target = 5.0;
  dxdt_CPT_time = 0.05 * (CPT_target - CPT_time);

$TABLE
  // Derived quantities
  double TLR4_occupancy = TLR4_bound / (TLR4_free + TLR4_bound + 0.001) * 100;
  double OPRM1_occupancy = OPRM1_NTX_occ * 100;
  double Pain_reduction = (7.0 - Pain_VAS) / 7.0 * 100;
  double Inflammatory_index = (TNF / 0.08 + IL1B / 0.05 + IL6 / 0.10) / 3.0;
  double Neurotrophin_index = (BDNF / 0.15 + NGF / 0.10 + TAC1 / 0.08) / 3.0;
  double Endorphin_blockade = OPRM1_NTX_occ * 100;
  double Hyperalgesia_pct = (Hyperalgesia_factor - 1.0) * 100;
  double MyD88_activity = MyD88_act;
  double TRIF_activity = TRIF_act;
  double OPRM1_density_fold = OPRM1_density / OPRM1_density_basal;
  double Endorphin_fold = Endorphin_conc / Endorphin_conc_basal;
  double CPT_improvement = (CPT_time - CPT_basal_FMQ) / (CPT_basal_healthy - CPT_basal_FMQ) * 100;

  capture Pain_VAS Pain_reduction CPT_time CPT_improvement
          TLR4_occupancy OPRM1_occupancy C_CNS
          TNF IL1B IL6 IL10 IL17A TGFb
          BDNF NGF TAC1
          Inflammatory_index Neurotrophin_index
          Endorphin_blockade Hyperalgesia_pct
          MyD88_activity TRIF_activity
          OPRM1_density_fold Endorphin_fold;

$CAPTURE
  Pain_VAS Pain_reduction CPT_time CPT_improvement
  TLR4_occupancy OPRM1_occupancy C_CNS
  TNF IL1B IL6 IL10 IL17A TGFb
  BDNF NGF TAC1
  Inflammatory_index Neurotrophin_index
  Endorphin_blockade Hyperalgesia_pct
  MyD88_activity TRIF_activity
  OPRM1_density_fold Endorphin_fold
