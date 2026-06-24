# Dose Titration Simulation
# Due Bruun 2021 protocol: 1.5->3->4.5->6mg weekly

source("qsp_model_functions.R")
library(ggplot2)

cat("=== Dose Titration Simulation ===\n")

# Estimated CNS concentrations for each dose level
dose_Ccns <- c("1.5mg"=0.003, "3.0mg"=0.007, "4.5mg"=0.01, "6.0mg"=0.015)

# Protocol: Titrate weekly
# Week 1: 1.5mg, Week 2: 3mg, Week 3: 4.5mg, Week 4-12: 6mg
titrate_schedule <- data.frame(
  week = 1:12,
  Ccns = c(0.003, 0.007, 0.01, rep(0.015, 9))
)

cat("Running titration protocol...\n")

# Run week by week, carrying state forward
p <- get_params()
y <- get_y0()
t_seq_day <- seq(0, 7*24, by=1)
titration_results <- list()

for (w in 1:12) {
  Ccns_w <- titrate_schedule$Ccns[w]
  y["Ccns"] <- Ccns_w
  out <- ode(y=y, times=t_seq_day, func=qsp_ode, parms=p, method="lsoda")
  df <- as.data.frame(out)
  df$week <- w
  df$Ccns <- Ccns_w
  df$time_days_total <- (w-1)*7 + df$time/24
  titration_results[[w]] <- df
  # Update initial state for next week
  state_cols <- names(y)
  last_row <- df[nrow(df), ]
  for (nm in state_cols) {
    if (nm %in% names(last_row)) y[nm] <- as.numeric(last_row[[nm]])
  }
}

titration_df <- bind_rows(titration_results) %>%
  mutate(
    Scenario="Titration 1.5->6mg",
    TLR4_occ=TLR4b/(TLR4f+TLR4b+0.001)*100,
    OPRM1_occ=Ccns/(0.001+Ccns)*100,
    Pain_red=(7-Pain)/7*100,
    OPRM1_fold=OPRM1d/1.0,
    Endo_fold=Endo/0.002
  )

# Run fixed doses for comparison
cat("Running fixed dose comparisons...\n")
fixed_results <- list()
for (nm in c("Placebo", "LDN 1.5mg", "LDN 4.5mg", "LDN 6.0mg")) {
  Ccns_val <- switch(nm, "Placebo"=0, "LDN 1.5mg"=0.003, "LDN 4.5mg"=0.01, "LDN 6.0mg"=0.015)
  df <- run_qsp(Ccns_val, t_end=12*7*24)
  df$Scenario <- nm
  fixed_results[[nm]] <- df
}
fixed_df <- bind_rows(fixed_results)

# Combine
all <- bind_rows(
  fixed_df %>% select(time_days, Scenario, Pain, Pain_red, TNF, IL6, IL10, IL17, CPT, OPRM1_fold, Endo_fold, TLR4_occ, OPRM1_occ),
  titration_df %>% select(time_days=time_days_total, Scenario, Pain, Pain_red, TNF, IL6, IL10, IL17, CPT, OPRM1_fold, Endo_fold, TLR4_occ, OPRM1_occ)
)

# Weekly averages
weekly <- all %>%
  mutate(week=ceiling(time_days/7)) %>%
  group_by(Scenario, week) %>%
  summarise(Pain=mean(Pain), Pain_red=mean(Pain_red), TNF=mean(TNF), IL6=mean(IL6),
            IL10=mean(IL10), CPT=mean(CPT), OPRM1_fold=mean(OPRM1_fold), Endo_fold=mean(Endo_fold),
            .groups="drop")

# Plots
p1 <- ggplot(weekly, aes(x=week, y=Pain, color=Scenario)) +
  geom_line(linewidth=1.2) + geom_point(size=2) +
  labs(title="Dose Titration vs Fixed Dose: Pain Trajectory",
       subtitle="Due Bruun 2021: 1.5->3->4.5->6mg weekly escalation",
       x="Week", y="Pain VAS") +
  theme_minimal() +
  scale_color_manual(values=c("Placebo"="gray50","LDN 1.5mg"="lightblue",
                               "LDN 4.5mg"="steelblue","LDN 6.0mg"="darkblue",
                               "Titration 1.5->6mg"="firebrick"))
ggsave("../output/figures/titration_pain.png", p1, width=12, height=7)

p2 <- ggplot(weekly, aes(x=week, y=CPT, color=Scenario)) +
  geom_line(linewidth=1.2) +
  labs(title="CPT: Titration vs Fixed Dose", x="Week", y="CPT (sec)") +
  theme_minimal() +
  scale_color_manual(values=c("Placebo"="gray50","LDN 4.5mg"="steelblue",
                               "LDN 6.0mg"="darkblue","Titration 1.5->6mg"="firebrick"))
ggsave("../output/figures/titration_cpt.png", p2, width=10, height=6)

p3 <- weekly %>%
  filter(Scenario %in% c("Placebo","LDN 4.5mg","Titration 1.5->6mg")) %>%
  ggplot(aes(x=week, y=OPRM1_fold, color=Scenario)) +
  geom_line(linewidth=1.2) +
  labs(title="OPRM1 Receptor Density: Titration", x="Week", y="Fold Change") +
  theme_minimal() +
  scale_color_manual(values=c("Placebo"="gray50","LDN 4.5mg"="steelblue","Titration 1.5->6mg"="firebrick"))
ggsave("../output/figures/titration_hormesis.png", p3, width=10, height=6)

write.csv(weekly, "../output/results/titration_weekly.csv", row.names=FALSE)

cat("\n=== Week 12 Summary ===\n")
w12 <- weekly %>% filter(week==12) %>% select(Scenario, Pain, Pain_red, CPT, OPRM1_fold, Endo_fold)
print(as.data.frame(w12))
cat("\nDose titration complete.\n")
