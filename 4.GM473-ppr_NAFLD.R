
library(TwoSampleMR)
library(ggplot2)
library(foreach)
library(psych)
library(doParallel)

# Set up parallel computing
numCores <- 48
cl <- makeCluster(numCores)
registerDoParallel(cl)


# Read exposure list
iddf <- read.table("473GM_ppr.txt", header = TRUE, sep = "\t")
bioid <- as.vector(iddf$id)

result <- data.frame()

# Perform MR analysis for each exposure-outcome pair
foreach(i = bioid, .errorhandling = "pass") %do% {

  expo_rt <- read_exposure_data(
    filename = paste0("clump/", i, ".txt"),
    sep = "\t",
    snp_col = "rsids2",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value"
  )

  expo_rt$samplesize.exposure <- 5959

  # Get all outcome files in the "ppr" folder
  ppr_files <- list.files("ppr", full.names = TRUE, pattern = "\\.txt\\.gz$")

  foreach(outcome_file = ppr_files, .errorhandling = "pass") %do% {

    outcome_file_name <- sub("\\.txt\\.gz$", "", basename(outcome_file))

    outc_rt <- read_outcome_data(
      snps = expo_rt$SNP,
      filename = outcome_file,
      sep = "\t",
      snp_col = "rs_id",
      beta_col = "beta",
      se_col = "standard_error",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      eaf_col = "effect_allele_frequency",
      pval_col = "p_value",
      samplesize_col = "n"
    )

    harm_rt <- harmonise_data(
      exposure_dat = expo_rt, 
      outcome_dat = outc_rt,
      action = 2
    )

    # Skip if no harmonized SNPs
    if (nrow(harm_rt) == 0) {
      next
    }

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
    result <- rbind(result, cbind(id = paste0(i, "-", outcome_file_name), result_or))

    if (mr_result$pval < 0.05) {
      filename <- paste0("result/", i, "-", outcome_file_name)
      dir.create(filename, showWarnings = FALSE)

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
}

# Save final combined result
write.table(result, "GM473-ppr_NAFLD-mr-all-result.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Stop parallel cluster
stopCluster(cl)
