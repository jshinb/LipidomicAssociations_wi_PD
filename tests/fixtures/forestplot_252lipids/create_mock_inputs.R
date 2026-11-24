# Tiny synthetic input data for smoke-testing the forest-plot workflow.
# These are fictitious lipid IDs and association results, not study observations.
lipid_classes <- data.frame(
  lipidID = c("TAG_16_0", "DAG_16_0", "CE_18_1", "PC_16_0_18_1", "PC_O-16_0_18_1", "LPC_16_0", "PE_18_0_18_1", "PE_O-18_0_18_1", "LPE_18_0", "PI_18_0_20_4", "Cer_d18_1_16_0", "SM_d18_1_16_0"),
  lipidClass = c("TAG", "DAG", "CE", "PC", "PC-O", "LPC", "PE", "PE-O", "LPE", "PI", "Cer", "SM")
)

make_mock_results <- function(offset = 0) {
  n <- nrow(lipid_classes)
  data.frame(
    lipidID = lipid_classes$lipidID,
    beta = seq(-0.18, 0.15, length.out = n) + offset,
    se = rep(0.045, n),
    n = rep(100L, n),
    pval = rep(c(0.01, 0.08, 0.03), length.out = n),
    fdr.p = rep(c(0.04, 0.12, 0.03), length.out = n)
  )
}

write.table(lipid_classes, "Class252lipids.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(make_mock_results(0), "assoc_lipid_GBA1status_any_Species-Class_adjBMI_2026-10-08.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(make_mock_results(0.025), "assoc_test_gba1eff_Species-Class_Controls_adjBMI_2026-10-15.txt", sep = "\t", quote = FALSE, row.names = FALSE)
