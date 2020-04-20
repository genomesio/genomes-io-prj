pkgs <- c("openxlsx", "jsonlite")
newPkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(newPkgs)) install.packages(newPkgs)
lapply(pkgs, library, character.only = TRUE)
library(tools)

#' Arguments
#' 1: full or test only
#' 2: input ID
#' 3: run dir (e.g. "imputation")
#' 4: output dir
#' 5: shapeit
#' 6: plink
#' 7: gtool
#' 8: impute2
#' 9: sample_ref
#' 10: imputation reference data dir
#' 11: impute-me dir

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 10) {
    msg <- paste(
        "You need to give all these args:",
		"1: run mode, either full or test",
        "2: input ID",
        "3: imputation dir",
        "4: output dir",
        "5: path to shapeit",
        "6: path to plink",
        "7: path to gtool",
        "8: path to impute2",
        "9: path to sample_ref",
        "10: path to imputation reference data",
		"11: path to the impute-me tool"
    )
    stop(paste0(msg, "\n"), call. = FALSE)
}

mode <- args[1]
uniqueID <- args[2]
runDir <- args[3]
destinationDir <- args[4]
shapeit <- args[5]
plink <- args[6]
gtool <- args[7]
impute2 <- args[8]
sample_ref <- args[9]
imputeDataDir <- args[10]
impute_me <- args[11]

# MAIN
source("imputation_fn.R")

if (mode == "full") {
    # check imputation dir
    if (!file.exists(paste0(runDir, "/", uniqueID, "_raw_data.txt"))) {
        stop(paste("No raw input file found in", runDir))
    }
    # step 1: run imputation
    startTime <- run_imputation(runDir, shapeit, plink, impute2, imputeDataDir, sample_ref)
    # step 2: summarize imputation results
    summarize_imputation(runDir, uniqueID, destinationDir, gtool, plink, startTime)
    # step 3: get SNPs to analyze
    crawl_for_snps_to_analyze(uniqueID, impute_me, destinationDir)
}

# step 4: export test results
allResults <- run_export_script(uniqueID = uniqueID, modules = NULL, impute_me = impute_me, destinationDir = destinationDir, gtool = gtool)
