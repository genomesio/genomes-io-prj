pkgs <- c("openxlsx", "jsonlite")
newPkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(newPkgs)) install.packages(newPkgs)
lapply(pkgs, library, character.only = TRUE)
library(tools)

#' Arguments
#' 1: input ID
#' 2: run dir (e.g. "imputation")
#' 3: output dir
#' 4: shapeit
#' 5: plink
#' 6: gtool
#' 7: impute2
#' 8: sample_ref
#' 9: imputation reference data dir
#' 10: impute-me dir

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 10) {
    msg <- paste(
        "You need to give all these args:",
        "1: input ID",
        "2: imputation dir",
        "3: output dir",
        "4: path to shapeit",
        "5: path to plink",
        "6: path to gtool",
        "7: path to impute2",
        "8: path to sample_ref",
        "9: path to imputation reference data",
		"10: path to the impute-me tool"
    )
    stop(paste0(msg, "\n"), call. = FALSE)
}

uniqueID <- args[1]
runDir <- args[2]
destinationDir <- args[3]
shapeit <- args[4]
plink <- args[5]
gtool <- args[6]
impute2 <- args[7]
sample_ref <- args[8]
imputeDataDir <- args[9]
impute_me <- args[10]

# check imputation dir
if (!file.exists(paste0(runDir, "/", uniqueID, "_raw_data.txt"))) {
    stop(paste("No raw input file found in", runDir))
}

# MAIN
source("imputation_fn.R")

# step 1: run imputation
run_imputation(runDir, shapeit, plink, impute2, imputeDataDir, sample_ref)
# step 2: summarize imputation results
summarize_imputation(runDir, uniqueID, destinationDir, gtool, plink)
# step 3: get SNPs to analyze
crawl_for_snps_to_analyze(uniqueID, impute_me, destinationDir)
# step 4: export test results
allResults <- run_export_script(uniqueID = uniqueID, modules = NULL, impute_me = impute_me, destinationDir = destinationDir, functionFile = file_path_as_absolute("imputation_fn.R"), gtool = gtool)
