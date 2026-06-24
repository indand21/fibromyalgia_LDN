# Data Preparation: Copy and format network pharmacology data for QSP model

library(dplyr)

# ============================================================
# 1. HUB GENE DATA (from Cytoscape degree analysis)
# ============================================================

hub_genes <- data.frame(
  Rank = c(1, 1, 3, 4, 4, 6, 7, 8, 9, 10),
  Gene = c("TNF", "IL6", "IL1B", "ALB", "BDNF", "TAC1", "IL10", "POMC", "NGF", "CXCL8"),
  Degree = c(45, 45, 42, 41, 41, 38, 36, 35, 34, 33),
  Pathway = c(
    "Inflammation", "Inflammation", "Inflammation", "Transport",
    "Neurotrophin", "Pain", "Inflammation", "HPA Axis",
    "Neurotrophin", "Inflammation"
  ),
  QSP_Module = c(
    "Neuroinflammation", "Neuroinflammation", "Neuroinflammation", "PK/Binding",
    "Pain Sensitization", "Pain Sensitization", "Neuroinflammation",
    "HPA Axis", "Pain Sensitization", "Neuroinflammation"
  )
)

write.csv(hub_genes, file.path("data", "hub_genes.csv"), row.names = FALSE)
cat("Hub genes saved to data/hub_genes.csv\n")

# ============================================================
# 2. FORMATTED PPI NETWORK (top interactions for QSP model)
# ============================================================

# Key interactions used in the QSP model ODEs
ppi_network <- data.frame(
  Source = c(
    "TLR4", "TLR4", "NFkB", "NFkB", "NFkB", "NFkB",
    "NLRP3", "IL10", "IL6", "TNF", "BDNF", "NGF",
    "CRH", "CRH", "TNF", "TNF", "IL1B", "IL1B",
    "NR3C1", "BDNF", "IL6", "HT"
  ),
  Target = c(
    "NFkB", "NLRP3", "TNF", "IL6", "IL1B", "CXCL8",
    "IL1B", "TNF", "BDNF", "NGF", "TAC1", "TRPV1",
    "POMC", "NPY", "CRH", "NR3C1", "TAC1", "TRPV1",
    "IL10", "TRPV1", "BDNF", "IL10"
  ),
  Effect = c(
    "activation", "activation", "activation", "activation", "activation", "activation",
    "activation", "inhibition", "activation", "activation", "activation", "activation",
    "activation", "inhibition", "activation", "inhibition", "activation", "activation",
    "activation", "activation", "activation", "activation"
  ),
  STRING_score = c(
    0.85, 0.75, 0.95, 0.95, 0.95, 0.80,
    0.90, 0.99, 0.85, 0.71, 0.84, 0.85,
    0.99, 0.95, 0.56, 0.75, 0.69, 0.62,
    0.53, 0.65, 0.85, 0.60
  ),
  QSP_equation = c(
    "NFkB = f(TLR4)", "NLRP3 = f(NFkB)",
    "dTNF/dt = k_prod * NFkB", "dIL6/dt = k_prod * NFkB",
    "dIL1B/dt = k_prod * NLRP3", "dCXCL8/dt = k_prod * NFkB",
    "dIL1B/dt = k_prod * NLRP3", "dTNF/dt -= k * IL10 * TNF",
    "dBDNF/dt += k * IL6", "dNGF/dt += k * TNF",
    "dTAC1/dt += k * BDNF", "TRPV1 = f(NGF)",
    "dPOMC/dt += k * CRH", "dNPY/dt -= k * CRH",
    "dCRH/dt += k * TNF", "NR3C1_act -= k * TNF",
    "dTAC1/dt += k * IL1B", "TRPV1 += k * IL1B",
    "dIL10/dt += k * NR3C1", "TRPV1 += k * BDNF",
    "dBDNF/dt += k * IL6", "dIL10/dt += k * HT"
  )
)

write.csv(ppi_network, file.path("data", "ppi_network.csv"), row.names = FALSE)
cat("PPI network saved to data/ppi_network.csv\n")

# ============================================================
# 3. TARGET RELEVANCE SCORES (from GeneCards)
# ============================================================

# Top targets with relevance to fibromyalgia (from GeneCards search)
target_scores <- data.frame(
  Gene = c(
    "COMT", "FMR1", "BDNF", "IGF1", "HTR3A", "SLC6A4", "HTR2A",
    "OPRM1", "IL6", "CXCL8", "TACR1", "POMC", "LEP", "CRP",
    "CCL2", "CRH", "TNF", "MAOA", "ADRB2", "NGF", "TSPO",
    "IL4", "NPY", "IL1B", "IL10", "SCN9A", "PNOC", "ESR1",
    "ACE", "MTHFR", "NOS3", "NR3C1", "TRPV1", "TRPV2", "TRPV3",
    "ASIC3", "TAC1", "NLRP3", "GSR", "TLR4"
  ),
  Relevance_score = c(
    13.34, 12.49, 11.50, 10.73, 10.28, 9.91, 9.87,
    9.53, 9.50, 9.38, 9.22, 9.07, 8.93, 8.93,
    8.60, 8.41, 8.41, 8.32, 8.24, 8.18, 8.13,
    8.11, 7.86, 7.75, 7.64, 7.37, 7.14, 1.99,
    1.99, 1.99, 2.15, 1.33, 0.71, 1.14, 1.14,
    0.71, 0.71, 1.33, 0.55, 0.55
  ),
  QSP_Module = c(
    "Monoamine", "Neurotransmission", "Pain", "Growth", "Serotonin",
    "Serotonin", "Serotonin", "Opioid", "Inflammation", "Inflammation",
    "Pain", "HPA Axis", "Metabolic", "Inflammation", "Inflammation",
    "HPA Axis", "Inflammation", "Monoamine", "Adrenergic", "Pain",
    "Neuroinflammation", "Inflammation", "HPA Axis", "Inflammation",
    "Inflammation", "Pain", "Pain", "Hormonal", "RAAS", "Metabolism",
    "NO pathway", "HPA Axis", "Pain", "Pain", "Pain",
    "Pain", "Pain", "Inflammation", "Oxidative Stress", "Innate Immune"
  ),
  In_QSP_Model = c(
    FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE,
    FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE,
    TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE,
    FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE,
    FALSE, TRUE, TRUE, FALSE, TRUE
  )
)

write.csv(target_scores, file.path("data", "genecards_scores.csv"), row.names = FALSE)
cat("Target relevance scores saved to data/genecards_scores.csv\n")

# ============================================================
# 4. COMMON TARGETS (from Venny analysis)
# ============================================================

common_targets <- readLines(file.path("..", "venny_result.txt"))
cat(sprintf("Loaded %d common targets from network pharmacology analysis\n", length(common_targets) - 1))

# ============================================================
# 5. PARAMETER REFERENCE TABLE
# ============================================================

param_references <- data.frame(
  Parameter = c(
    "Ka", "F", "Vd", "CL", "Kp_brain",
    "kon_NTX", "koff_NTX", "Kd_NTX",
    "k_deg_TNF", "k_deg_IL6", "k_deg_IL1B", "k_deg_IL10",
    "k_deg_BDNF", "k_deg_NGF", "k_deg_TAC1",
    "EC50_NGF_TRPV1", "n_pain"
  ),
  Value = c(
    "1.2 1/hr", "0.05", "1127 L", "22.6 L/hr", "0.028",
    "0.1 1/(uM*hr)", "1.0 1/hr", "10 uM",
    "0.693 1/hr", "0.116 1/hr", "0.347 1/hr", "0.116 1/hr",
    "0.023 1/hr", "0.023 1/hr", "0.693 1/hr",
    "0.05 uM", "1.5"
  ),
  Source = c(
    "FDA Label (ReVia)", "Literature (low-dose)", "FDA Label", "FDA Label", "Literature",
    "Molecular docking (estimated)", "Molecular docking (estimated)", "Kd = koff/kon",
    "Cytokine biology literature", "Cytokine biology literature",
    "Cytokine biology literature", "Cytokine biology literature",
    "Neurotrophin literature", "Neurotrophin literature", "Neuropeptide literature",
    "TRPV1 electrophysiology", "Pain modeling literature"
  ),
  Notes = c(
    "First-order oral absorption", "Reduced at low doses",
    "Total body water + tissue", "Hepatic + renal", "Low CNS penetration",
    "From AutoDock Vina score", "Estimated from binding energy",
    "Derived from kon/koff", "TNF t1/2 ~1 hr", "IL6 t1/2 ~6 hr",
    "IL1B t1/2 ~2 hr", "IL10 t1/2 ~6 hr", "BDNF t1/2 ~30 hr",
    "NGF t1/2 ~30 hr", "Substance P t1/2 ~1 hr",
    "From patch-clamp studies", "Empirical fitting"
  )
)

write.csv(param_references, file.path("data", "parameter_references.csv"), row.names = FALSE)
cat("Parameter references saved to data/parameter_references.csv\n")

cat("\n=== Data Preparation Complete ===\n")
