# QSP Model: Naltrexone Mechanism in Fibromyalgia

## 1. Project Overview

**Objective:** Build a Quantitative Systems Pharmacology (QSP) model to mechanistically characterize how Low-Dose Naltrexone (LDN) exerts its therapeutic effects in Fibromyalgia through TLR4 antagonism and downstream modulation of neuroimmune pathways.

**Rationale:** The network pharmacology analysis identified 70 overlapping targets between Fibromyalgia and potential therapeutics, with TNF, IL6, IL1B, BDNF, and TAC1 as hub genes. Naltrexone's docking to TLR4 positions it as a modulator of the neuroimmune interface central to fibromyalgia pathophysiology.

---

## 2. Pathophysiology Modules (from Network Pharmacology Data)

### Module 1: Neuroinflammation (TLR4-NLRP3 Axis)
| Target | Role | Network Degree |
|--------|------|----------------|
| TLR4 | Innate immune sensor; Naltrexone binding target | - |
| NLRP3 | Inflammasome activation | 28 |
| TNF | Pro-inflammatory cytokine (hub) | 45 |
| IL1B | Pro-inflammatory cytokine (hub) | 42 |
| IL6 | Pro-inflammatory cytokine (hub) | 45 |
| IL10 | Anti-inflammatory cytokine | 36 |
| CXCL8 | Chemokine (neutrophil recruitment) | 33 |
| CCL2 | Monocyte chemotaxis | 24 |

**ODEs:** TLR4 activation -> NLRP3 assembly -> Caspase-1 -> IL1B maturation. TNF/IL6 amplify via NF-kB feedback. IL10 provides negative regulation.

### Module 2: Central Sensitization & Pain Signaling
| Target | Role | Network Degree |
|--------|------|----------------|
| BDNF | Neurotrophin; synaptic plasticity | 41 |
| NGF | Nociceptor sensitization | 34 |
| TAC1 (Substance P) | Pain neurotransmitter | 38 |
| TACR1 | NK1 receptor for Substance P | 17 |
| TRPV1 | Nociceptive ion channel | 22 |
| SCN9A | Nav1.7 sodium channel | 15 |
| ASIC3 | Acid-sensing ion channel | 12 |

**ODEs:** NGF/BDNF -> TrkA/TrkB signaling -> TRPV1 phosphorylation -> Ca2+ influx -> TAC1 release. Positive feedback via IL1B/TNF sensitization of nociceptors.

### Module 3: HPA Axis & Stress Response
| Target | Role | Network Degree |
|--------|------|----------------|
| CRH | Corticotropin-releasing hormone | 28 |
| POMC | Proopiomelanocortin (ACTH precursor) | 35 |
| NR3C1 | Glucocorticoid receptor | 25 |
| NPY | Neuropeptide Y (anxiolytic) | 24 |
| LEP | Leptin (metabolic-inflammatory link) | 22 |

**ODEs:** CRH -> POMC -> cortisol -> NR3C1 activation -> anti-inflammatory gene expression. Chronic inflammation blunts NR3C1 sensitivity (glucocorticoid resistance).

### Module 4: Monoamine Neurotransmission
| Target | Role | Network Degree |
|--------|------|----------------|
| SLC6A4 | Serotonin transporter | 20 |
| HTR2A | 5-HT2A receptor | 18 |
| COMT | Catecholamine degradation | 22 |
| MAOA/MAOB | Monoamine degradation | 20/18 |
| OPRM1 | Mu-opioid receptor | 19 |
| DRD3 | Dopamine D3 receptor | 12 |

**ODEs:** Tonic serotonin/dopamine levels modulated by transporter kinetics and enzymatic degradation. OPRM1 signaling intersects with pain and reward circuits.

---

## 3. Drug Mechanism: Naltrexone at TLR4

### Pharmacokinetics (PK)
- **Route:** Oral
- **Dose:** 4.5 mg (Low-Dose Naltrexone)
- **Absorption:** Ka, F (oral bioavailability)
- **Distribution:** Vd (plasma + CNS compartment)
- **Metabolism:** CYP2D6 -> 6-beta-naltrexol (active metabolite)
- **Elimination:** CL (renal + hepatic)

### Pharmacodynamics (PD) at TLR4
```
Naltrexone + TLR4 <-> [Naltrexone-TLR4 complex]
    -> Inhibits TLR4/MD2 dimerization
    -> Reduces NF-kB nuclear translocation
    -> Decreases TNF, IL1B, IL6, CXCL8 transcription
    -> Increases IL10 (anti-inflammatory shift)
    -> Reduces NLRP3 inflammasome priming
```

**Binding kinetics (from docking):**
- Kd (TLR4-Naltrexone): ~10-50 uM (estimated from docking scores)
- IC50 for TLR4 inhibition: to be calibrated with clinical data

### Dose-Dependent Receptor Selectivity (WHY ONLY LOW DOSE WORKS)

This is the central mechanistic question addressed by the model.

**The Affinity Problem:**
| Receptor | Naltrexone Kd | Role in Pain |
|----------|---------------|--------------|
| OPRM1 (mu-opioid) | **0.001 uM (1 nM)** | Endogenous analgesia |
| TLR4 (innate immune) | **10 uM** | Neuroinflammation |

Naltrexone has **~10,000x higher affinity** for OPRM1 than TLR4.

**Dose-Dependent Outcomes:**

| Dose | CNS Conc | TLR4 Occ | OPRM1 Occ | Net Effect |
|------|----------|----------|-----------|------------|
| LDN 1.5mg | ~0.003 uM | <0.1% | ~75% | Minimal |
| LDN 4.5mg | ~0.01 uM | ~0.1% | ~91% | **Anti-inflammatory (beneficial)** |
| Standard 25mg | ~0.5 uM | ~5% | ~99.8% | Mixed |
| Standard 50mg | ~1 uM | ~9% | ~99.9% | **Hyperalgesia (harmful)** |

**The Paradox Explained:**
1. At LDN doses: TLR4 antagonism reduces neuroinflammation -> pain relief
2. At standard doses: OPRM1 blockade -> loss of endogenous endorphin signaling -> compensatory hyperalgesia -> pain WORSENS
3. This creates a **U-shaped dose-response curve** where 4.5mg is optimal

**Model Implementation:**
- Competitive binding: Naltrexone vs endogenous endorphins at OPRM1
- Paradoxical hyperalgesia: cubic threshold function of OPRM1 blockade
- Pain = Inflammatory_component * Opioid_modulation * Hyperalgesia_factor

---

## 4. Model Structure

```
                    +------------------+
                    |   PK Compartment  |
                    |  (Gut, Plasma,   |
                    |   CNS, Periph)   |
                    +--------+---------+
                             |
                    [Naltrexone] in CNS
                             |
                    +--------v---------+
                    |  TLR4 Antagonism  |
                    |  (Microglia)      |
                    +--------+---------+
                             |
              +--------------+--------------+
              |                             |
     +--------v--------+          +--------v--------+
     |  NF-kB Pathway   |          |  NLRP3 Inflam.  |
     |  (TNF, IL6, IL1B)|          |  (IL1B mature)   |
     +--------+--------+          +--------+--------+
              |                             |
              +--------------+--------------+
                             |
                    +--------v---------+
                    |  Neuroinflammation |
                    |  Reduction         |
                    +--------+---------+
                             |
              +--------------+--------------+
              |              |              |
     +--------v--+  +--------v--+  +--------v--+
     | Pain Signal|  | HPA Axis  |  | Monoamine |
     | Reduction  |  | Restore   |  | Balance   |
     +------------+  +-----------+  +-----------+
              |              |              |
              +--------------+--------------+
                             |
                    +--------v---------+
                    | Clinical Outcome:  |
                    | Pain Score (VAS)   |
                    | FIQ, Fatigue, Mood |
                    +------------------+
```

---

## 5. Key Equations (ODE Framework)

### 5.1 PK Model (2-compartment with first-order absorption)
```
dA_gut/dt    = -Ka * A_gut
dC_plasma/dt = Ka * A_gut / Vd - (CL/Vd + Q/Vd) * C_plasma + (Q/V2) * C_periph
dC_periph/dt = (Q/V2) * C_plasma - (Q/V2) * C_periph
dC_CNS/dt    = Kin * C_plasma - Kout * C_CNS - kon_TLR4 * C_CNS * TLR4_free + koff_TLR4 * TLR4_bound
```

### 5.2 TLR4 Signaling
```
dTLR4_free/dt   = -kon * C_CNS * TLR4_free + koff * TLR4_bound
dTLR4_bound/dt  =  kon * C_CNS * TLR4_free - koff * TLR4_bound
NFkB_activity   = f(TLR4_total - TLR4_bound)  // Hill function
```

### 5.3 Cytokine Dynamics
```
dTNF/dt  = k_prod_TNF * NFkB_activity - k_deg_TNF * TNF + k_amp * TNF * (1 - TNF/TNF_max) - k_IL10_inhib * IL10 * TNF
dIL1B/dt = k_prod_IL1B * NLRP3_activity - k_deg_IL1B * IL1B
dIL6/dt  = k_prod_IL6 * NFkB_activity - k_deg_IL6 * IL6
dIL10/dt = k_prod_IL10 * (1 + NR3C1_activity) - k_deg_IL10 * IL10
dCXCL8/dt= k_prod_CXCL8 * NFkB_activity - k_deg_CXCL8 * CXCL8
```

### 5.4 Pain Output (Sigmoidal Integration)
```
Pain_score = Pain_baseline * (1 + Emax_BDNF * BDNF^n / (EC50_BDNF^n + BDNF^n))
                            * (1 + Emax_TNF * TNF^n / (EC50_TNF^n + TNF^n))
                            * (1 + Emax_TAC1 * TAC1^n / (EC50_TAC1^n + TAC1^n))
                            / (1 + E_IL10 * IL10 / (EC50_IL10 + IL10))
```

---

## 6. Parameter Sources

| Parameter Category | Source |
|-------------------|--------|
| PK parameters | Literature (naltrexone PK studies), FDA label |
| Binding affinities | Molecular docking results (nalt_tlr4_docking.sam) |
| Cytokine kinetics | Published PBMC stimulation assays |
| Network topology | STRING PPI data (string_interactions_short.tsv) |
| Hub gene ranking | Cytoscape degree analysis (cytoscape_rank.csv) |
| Target relevance | GeneCards relevance scores |
| Clinical endpoints | Fibromyalgia RCT data (literature) |

---

## 7. Simulation Plan

1. **Baseline fibromyalgia state:** Elevate TNF, IL6, IL1B, BDNF, TAC1 above healthy controls
2. **Naltrexone PK/PD simulation:** Single dose and steady-state (4.5 mg QD x 28 days)
3. **Dose-dependent selectivity analysis (KEY):** Compare LDN (4.5mg) vs standard (50mg)
   - TLR4 occupancy across doses
   - OPRM1 occupancy across doses
   - Endorphin blockade quantification
   - Hyperalgesia threshold identification
   - U-shaped dose-response curve generation
4. **Sensitivity analysis:** Identify parameters with greatest impact on pain reduction
5. **Dose-response exploration:** 1.5 mg, 3.0 mg, 4.5 mg, 9.0 mg, 25 mg, 50 mg
6. **Virtual patient population:** Vary baseline cytokine levels and receptor densities

---

## 8. Deliverables

- [ ] R/mrgsolve QSP model code
- [ ] Parameter table with literature references
- [ ] Baseline disease state calibration
- [ ] PK/PD simulation with naltrexone
- [ ] Sensitivity analysis (local + global)
- [ ] Dose-response curves
- [ ] Clinical trial simulation framework

---

## 9. File Structure

```
QSP_Fibromyalgia_Naltrexone/
  PROJECT_PLAN.md          <- This file
  model/
    qsp_model.R            <- Main mrgsolve model file
    parameters.R           <- Parameter definitions with references
    initial_conditions.R   <- Baseline disease state
  data/
    ppi_network.csv        <- STRING interaction data (from network pharmacology)
    hub_genes.csv          <- Degree-ranked targets
    genecards_scores.csv   <- Target relevance scores
  scripts/
    run_simulation.R       <- Main simulation script
    sensitivity_analysis.R <- Local + global sensitivity
    dose_response.R        <- Dose-response exploration
    plot_results.R         <- Visualization
  output/
    figures/               <- Generated plots
    results/               <- Simulation results (CSV)
```
