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
#' 11: impute-functions dir

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 10) {
    msg <- paste(
        "You need to give all these args:",
        "1; trail names for report",
		"2: run mode, either full or test",
        "3: input ID",
        "4: imputation dir",
        "5: output dir",
        "6: path to shapeit",
        "7: path to plink",
        "8: path to gtool",
        "9: path to impute2",
        "10: path to sample_ref",
        "11: path to imputation reference data",
		"12: path to the impute-functions tool"
    )
    stop(paste0(msg, "\n"), call. = FALSE)
}

trails <- args[1]
mode <- args[2]
uniqueID <- args[3]
runDir <- args[4]
destinationDir <- args[5]
shapeit <- args[6]
plink <- args[7]
gtool <- args[8]
impute2 <- args[9]
sample_ref <- args[10]
imputeDataDir <- args[11]
imputeTrails <- args[12]

# MAIN
source("imputation_fn.R")

if (mode == "full" || mode == "impute") {
    # check imputation dir
    if (!file.exists(paste0(runDir, "/", uniqueID, "_raw_data.txt"))) {
        stop(paste("No raw input file found in", runDir))
    }
    # step 1: run imputation
    startTime <- run_imputation(runDir, shapeit, plink, impute2, imputeDataDir, sample_ref)
    # step 2: summarize imputation results
    summarize_imputation(runDir, uniqueID, destinationDir, gtool, plink, startTime)
    # step 3: get SNPs to analyze
    crawl_for_snps_to_analyze(uniqueID, imputeTrails, destinationDir)
}

# step 4: export test results
if (mode == "full" || mode == "test") {
    allResults <- run_export_script(
        uniqueID = uniqueID, 
        modules = trails, 
        imputeTrails = imputeTrails, 
        destinationDir = destinationDir, 
        gtool = gtool
    )
}

