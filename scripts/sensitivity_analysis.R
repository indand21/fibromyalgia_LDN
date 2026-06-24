# Sensitivity Analysis: Local (OAT) + Global (PRCC)
# QSP Model v2: Naltrexone in Fibromyalgia

source("qsp_model_functions.R")
library(lhs)
library(ggplot2)

# ============================================================
# 1. LOCAL SENSITIVITY (One-at-a-Time)
# ============================================================

cat("=== Local Sensitivity Analysis ===\n")

sa_params <- list(
  kon = c(0.01, 0.05, 0.1, 0.2, 0.5),
  koff = c(0.1, 0.5, 1.0, 2.0, 5.0),
  kT_TNF = c(0.5, 1.0, 1.5, 2.5, 4.0),
  kT_IL6 = c(0.3, 0.7, 1.0, 1.5, 2.5),
  kT_IL17 = c(0.3, 0.7, 1.0, 1.5, 2.5),
  k_up = c(0.01, 0.03, 0.05, 0.08, 0.12),
  k_reb = c(0.005, 0.01, 0.02, 0.04, 0.08),
  Kd_OPRM1 = c(0.0002, 0.0005, 0.001, 0.002, 0.005),
  k_hyper = c(0.1, 0.3, 0.5, 0.8, 1.2),
  kIL6_BDNF = c(0.1, 0.3, 0.5, 0.8, 1.2),
  kB_TAC1 = c(0.3, 0.7, 1.0, 1.5, 2.5),
  Em_BDNF = c(0.2, 0.4, 0.5, 0.7, 1.0),
  Em_TNF = c(0.1, 0.2, 0.3, 0.5, 0.8),
  kdTNF = c(0.2, 0.5, 0.693, 1.0, 2.0),
  kdIL6 = c(0.05, 0.08, 0.116, 0.2, 0.4)
)

run_sa <- function(param_name, values, Ccns=0.01) {
  results <- list()
  base_p <- get_params()
  for (v in values) {
    p_mod <- base_p
    p_mod[[param_name]] <- v
    df <- run_qsp(Ccns, p=p_mod, t_end=28*24)
    ss <- get_ss(df)
    results[[as.character(v)]] <- data.frame(
      Parameter=param_name, Value=v,
      Pain_VAS=ss$Pain, Pain_red=ss$Pain_red,
      TNF=ss$TNF, IL6=ss$IL6, BDNF=ss$BDNF,
      CPT=ss$CPT, OPRM1_fold=ss$OPRM1_fold
    )
  }
  bind_rows(results)
}

local_sa <- list()
for (nm in names(sa_params)) {
  cat(sprintf("  %s\n", nm))
  local_sa[[nm]] <- run_sa(nm, sa_params[[nm]])
}
local_df <- bind_rows(local_sa)

si <- local_df %>%
  group_by(Parameter) %>%
  arrange(Value) %>%
  mutate(
    dPain = (Pain_VAS - lag(Pain_VAS))/lag(Pain_VAS),
    dParam = (Value - lag(Value))/lag(Value),
    SI = dPain/dParam
  ) %>%
  filter(!is.na(SI)) %>%
  summarise(Mean_SI=mean(abs(SI)), Max_SI=max(abs(SI)), .groups="drop") %>%
  arrange(desc(Mean_SI))

cat("\n=== Sensitivity Index Ranking (Pain VAS) ===\n")
print(as.data.frame(si))

tornado <- si %>% head(12) %>%
  mutate(Parameter=factor(Parameter, levels=rev(Parameter)))

p_tornado <- ggplot(tornado, aes(x=Parameter, y=Mean_SI)) +
  geom_bar(stat="identity", fill="steelblue", width=0.7) + coord_flip() +
  labs(title="Local Sensitivity: Pain VAS", x="Parameter", y="|Sensitivity Index|") +
  theme_minimal()
ggsave("../output/figures/tornado_plot.png", p_tornado, width=10, height=7)

spider_params <- si$Parameter[1:5]
p_spider <- local_df %>%
  filter(Parameter %in% spider_params) %>%
  ggplot(aes(x=Value, y=Pain_VAS, color=Parameter)) +
  geom_line(linewidth=1.2) + geom_point(size=2) +
  facet_wrap(~Parameter, scales="free_x", ncol=3) +
  labs(title="Spider Plots: Top 5 Parameters", x="Value", y="Pain VAS") +
  theme_minimal() + theme(legend.position="none")
ggsave("../output/figures/spider_plots.png", p_spider, width=14, height=8)

# ============================================================
# 2. GLOBAL SENSITIVITY (LHS + PRCC)
# ============================================================

cat("\n=== Global Sensitivity Analysis (LHS + PRCC) ===\n")

param_ranges <- data.frame(
  kon=c(0.01,0.5), koff=c(0.1,5.0), kT_TNF=c(0.5,4.0), kT_IL6=c(0.3,2.5),
  k_up=c(0.01,0.12), k_reb=c(0.005,0.08), Kd_OPRM1=c(0.0002,0.005),
  k_hyper=c(0.1,1.2), kIL6_BDNF=c(0.1,1.2), kB_TAC1=c(0.3,2.5),
  kdTNF=c(0.2,2.0), kdIL6=c(0.05,0.4)
)

n_samples <- 300
cat(sprintf("  Running %d LHS samples...\n", n_samples))

lhs_s <- randomLHS(n_samples, ncol(param_ranges))
colnames(lhs_s) <- names(param_ranges)
pmat <- matrix(NA, n_samples, ncol(param_ranges))
colnames(pmat) <- names(param_ranges)
for (i in 1:ncol(param_ranges))
  pmat[,i] <- param_ranges[1,i] + lhs_s[,i]*(param_ranges[2,i]-param_ranges[1,i])

gres <- list()
for (i in 1:n_samples) {
  if(i%%50==0) cat(sprintf("  %d/%d\n", i, n_samples))
  p_mod <- get_params()
  for (j in 1:ncol(param_ranges)) p_mod[[colnames(param_ranges)[j]]] <- pmat[i,j]
  df <- run_qsp(0.01, p=p_mod, t_end=28*24)
  ss <- get_ss(df)
  gres[[i]] <- cbind(data.frame(sample=i, Pain_VAS=ss$Pain, Pain_red=ss$Pain_red,
                                 TNF=ss$TNF, IL6=ss$IL6, BDNF=ss$BDNF, CPT=ss$CPT),
                      as.data.frame(t(pmat[i,])))
}
global_df <- bind_rows(gres)

compute_prcc <- function(df, out) {
  ps <- names(param_ranges)
  cors <- sapply(ps, function(p) cor.test(df[[p]], df[[out]], method="spearman")$estimate)
  data.frame(Parameter=ps, PRCC=cors, Abs_PRCC=abs(cors), Output=out) %>% arrange(desc(Abs_PRCC))
}

prcc_pain <- compute_prcc(global_df, "Pain_VAS")
cat("\n=== PRCC for Pain VAS ===\n")
print(as.data.frame(prcc_pain))

p_prcc <- prcc_pain %>%
  mutate(Parameter=factor(Parameter, levels=rev(Parameter))) %>%
  ggplot(aes(x=Parameter, y=PRCC, fill=PRCC>0)) +
  geom_bar(stat="identity", width=0.7) + coord_flip() +
  labs(title="Global Sensitivity (PRCC): Pain VAS", x="Parameter", y="PRCC") +
  theme_minimal() +
  scale_fill_manual(values=c("TRUE"="firebrick","FALSE"="steelblue"), guide="none")
ggsave("../output/figures/prcc_global_sa.png", p_prcc, width=10, height=6)

write.csv(local_df, "../output/results/local_sa_results.csv", row.names=FALSE)
write.csv(si, "../output/results/sensitivity_indices.csv", row.names=FALSE)
write.csv(global_df, "../output/results/global_sa_results.csv", row.names=FALSE)
write.csv(prcc_pain, "../output/results/prcc_pain.csv", row.names=FALSE)

cat("\nSensitivity analysis complete.\n")
