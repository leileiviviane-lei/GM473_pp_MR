
rm(list = ls())

library(psych)
library(doParallel)
library(TwoSampleMR)
library(ggplot2)
library(foreach)

# Set up parallel computing
numCores <- 48
cl <- makeCluster(numCores)
registerDoParallel(cl)


# Read exposure list
iddf <- read.table("id2821.txt", header = TRUE, sep = "\t")
bioid <- as.vector(iddf$id)

result <- data.frame()

# Perform MR analysis for each exposure
foreach(i = bioid, .errorhandling = "pass") %do% {

  expo_rt <- read_exposure_data(
    filename = paste0("clumped_5e-8/", i, ".csv"),
    sep = ",",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure",
    pval_col = "pval.exposure",
    samplesize_col = "samplesize.exposure"
  )

  outc_rt <- read_outcome_data(
    snps = expo_rt$SNP,
    filename = "finngen_R12_NAFLD.gz",
    sep = "\t",
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    eaf_col = "af_alt",
    pval_col = "pval"
  )

  # Add sample size for outcome data
  outc_rt$samplesize.outcome <- 500348

  # Harmonize data
  harm_rt <- harmonise_data(
    exposure_dat = expo_rt, 
    outcome_dat = outc_rt,
    action = 2
  )

  # Calculate F-statistics
  harm_rt$R2 <- (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure)) /
    (2 * (harm_rt$beta.exposure^2) * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) +
       2 * harm_rt$samplesize.exposure * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) * harm_rt$se.exposure^2)
  harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
  harm_rt$meanf <- mean(harm_rt$f)
  harm_rt <- harm_rt[harm_rt$f > 10, ]

  # Perform MR analysis
  mr_result <- mr(harm_rt)
  result_or <- generate_odds_ratios(mr_result)

  # Save MR result
  result <- rbind(result, cbind(id = i, pvalue = result_or$pval[3]))

  if (mr_result$pval[3] < 0.05) {
    filename <- paste0("result/", i)
    dir.create(filename)

    write.table(harm_rt, file = paste0(filename, "/harmonise.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
    write.table(result_or[, 5:ncol(result_or)], file = paste0(filename, "/OR.txt"), row.names = FALSE, sep = "\t", quote = FALSE)

    # Pleiotropy test
    pleiotropy <- mr_pleiotropy_test(harm_rt)
    write.table(pleiotropy, file = paste0(filename, "/pleiotropy.txt"), sep = "\t", quote = FALSE)

    # Heterogeneity test
    heterogeneity <- mr_heterogeneity(harm_rt)
    write.table(heterogeneity, file = paste0(filename, "/heterogeneity.txt"), sep = "\t", quote = FALSE)

    # Single SNP analysis
    singlesnp_res <- mr_singlesnp(harm_rt)
    singlesnpOR <- generate_odds_ratios(singlesnp_res)
    write.table(singlesnpOR, file = paste0(filename, "/singlesnpOR.txt"), row.names = FALSE, sep = "\t", quote = FALSE)

    # Leave-one-out analysis
    sen_res <- mr_leaveoneout(harm_rt)

    # MR-PRESSO
    presso <- run_mr_presso(harm_rt, NbDistribution = 1000)
    capture.output(presso, file = paste0(filename, "/presso.txt"))

    # MR Steiger directionality test
    harm_rt$r.exposure <- sqrt(harm_rt$R2)
    harm_rt$r.outcome <- with(harm_rt, sqrt(2 * (as.numeric(beta.outcome)^2) * eaf.outcome * (1 - eaf.outcome) /
      (2 * samplesize.outcome * eaf.outcome * (1 - eaf.outcome) * as.numeric(se.outcome)^2)))

    steiger_test <- directionality_test(harm_rt)
    write.table(steiger_test, file = paste0(filename, "/steiger.txt"), sep = "\t", quote = FALSE)
  }
}

# Adjust p-values
result$fdr <- p.adjust(result$pvalue, method = "BH")

# Save final result
write.table(result, "pp-finngen_R12_NAFLD-mr-all-results_FDR.txt", sep = "\t", quote = FALSE, row.names = FALSE)
