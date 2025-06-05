
# Load libraries
library(TwoSampleMR)
library(ggplot2)


### Two-step Mendelian Randomization Mediation Analysis ###
# Exposure: Gut Microbiota (GM473)
# Mediator: Protein-to-Protein Ratio (PPR)
# Outcome: Non-Alcoholic Fatty Liver Disease (NAFLD)

# Step 1: GM473 to NAFLD
harm_rt <- read.table("harmonise1.txt", sep = "\t", header = TRUE)
mr_result <- mr(harm_rt)
result_or <- generate_odds_ratios(mr_result)
write.table(result_or[, 5:ncol(result_or)], "OR1.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# Step 2: GM473 to PPR
harm_rt3 <- read.table("harmonise2.txt", sep = "\t", header = TRUE)
mr_result3 <- mr(harm_rt3)
result_or3 <- generate_odds_ratios(mr_result3)
write.table(result_or3[, 5:ncol(result_or3)], "OR3.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# Step 3: PPR to NAFLD
harm_rt4 <- read.table("harmonise3.txt", sep = "\t", header = TRUE)
mr_result4 <- mr(harm_rt4)
result_or4 <- generate_odds_ratios(mr_result4)
write.table(result_or4[, 5:ncol(result_or4)], "OR4.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# Mediation effect calculation
# Total effect: GM473 to NAFLD
beta_all <- mr_result[3, "b"]

# Indirect path:
# GM473 to PPR (beta1), and PPR to NAFLD (beta2)
beta1 <- mr_result3[1, "b"]
beta2 <- mr_result4[3, "b"]

# Indirect effect (mediation effect)
beta12 <- beta1 * beta2

# Direct effect
beta_direct <- beta_all - beta12

# Standard error for indirect effect
se <- sqrt(mr_result3[1, "b"]^2 * mr_result3[1, "se"]^2 + mr_result4[3, "b"]^2 * mr_result4[3, "se"]^2)

# Z-statistic and p-value for mediation effect
Z <- beta12 / se
P <- 2 * pnorm(q = abs(Z), lower.tail = FALSE)

# 95% Confidence Interval
lci <- beta12 - 1.96 * se
uci <- beta12 + 1.96 * se

# Proportion mediated
proportion_mediated <- beta12 / beta_all
proportion_lci <- lci / beta_all
proportion_uci <- uci / beta_all
