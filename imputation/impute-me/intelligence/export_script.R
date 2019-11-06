export_function <- function (uniqueID, moduleDir, outputDir, functionFile, gtool) {
    if (!file.exists(outputDir)) {
        stop(paste("Did not find a output data with this id", uniqueID))
    }
	# source(functionFile)

    snps_file <- paste0(moduleDir, "/intelligence/2019-03-04_semi_curated_version_gwas_central.rdata")
    trait_file <- paste0(moduleDir, "/intelligence/2019-03-04_trait_overview.xlsx")
    all_snp_trait_file <- paste0(moduleDir, "/prs/2019-03-05_study_list.xlsx")
    all_snp_traits <- read.xlsx(all_snp_trait_file, rowNames = TRUE)

    library(openxlsx)
    load(snps_file)
    traits  <- read.xlsx(trait_file, rowNames = TRUE)
    traits <- traits[!is.na(traits[, "omit"]) & !traits[, "omit"], ]

    output <- list()
    # output[["documentation"]]  <- list()
    # output[["documentation"]][["trait_overview"]]  <- "https://github.com/lassefolkersen/impute-me/blob/589f332a148e7c0f6041637bd1c97ec0de1a14ee/intelligence/2019-03-04_trait_overview.xlsx"
    # output[["documentation"]][["snp_file"]]  <- "https://github.com/lassefolkersen/impute-me/blob/589f332a148e7c0f6041637bd1c97ec0de1a14ee/intelligence/2019-03-04_semi_curated_version_gwas_central.rdata"

    # get ethnicity parameter
    pDataFile <- paste(outputDir, "/pData.txt", sep = "")
    pData <- try(read.table(pDataFile, header = TRUE, stringsAsFactors = FALSE, sep = "\t"))
    if (class(pData) != "try-error" && "ethnicity" %in% colnames(pData)) {
        ethnicity <- pData[1, "ethnicity"]
    } else {
        ethnicity <- "global"
    }

    for (study_id in rownames(traits)) {
        SNPs_to_analyze <- data[data[, "study_id"]%in%study_id , ]
        # get genotypes
        SNPs_requested <- SNPs_to_analyze[!duplicated(SNPs_to_analyze[, "SNP"]), ]
        rownames(SNPs_requested) <- SNPs_requested[, "SNP"]
        genotypes <- get_genotypes(uniqueID = uniqueID, request = SNPs_requested, gtool = gtool, destinationDir = outputDir, namingLabel = "cached.all_gwas")
        # get correct ethnicity minor_allele_frequency
        if(ethnicity %in% c("EAS", "AMR", "AFR", "EUR", "SAS")){
            SNPs_requested[, "minor_allele_freq"] <- SNPs_requested[, paste0(ethnicity, "_AF")]
        }
        # calculate GRS
        snp_data <- SNPs_requested
        snp_data[, "genotype"]  <- genotypes[rownames(snp_data), "genotype"]
        snp_data <- get_GRS_2(snp_data, mean_scale = TRUE, unit_variance = TRUE, verbose = FALSE)
        population_sum_sd <- sqrt(sum(snp_data[, "population_score_sd"]^2, na.rm = TRUE))
        GRS_beta <- sum(snp_data[, "score_diff"], na.rm = TRUE) / population_sum_sd
        # calculate percentage
        percentage <- floor(pnorm(GRS_beta, mean = 0, sd = 1)*100)
        # calculate risk-allele
        c1 <- apply(SNPs_requested[,c("minor_allele","major_allele","effect_allele","non_effect_allele")]=="?", 1, sum)

        # gather some background info for the study
        trait <- traits[study_id, "trait"]
        pmid <- traits[study_id, "pmid"]
        link <- paste0("www.ncbi.nlm.nih.gov/pubmed/", traits[study_id, "pmid"])
        author <- traits[study_id, "first_author"]
        sampleSize <- traits[study_id, "sampleSize"]
        publication_date <- traits[study_id, "publication_date"]

        # message
        textToReturn <- paste0("Ethnicity-corrected trait Z-score is ", signif(GRS_beta, 2))
        textToReturn <- paste0(textToReturn, " This genetic risk score is higher than ", percentage, "% of the general population.")
        if (!is.na(percentage)) {
            if (percentage < 20){
                textToReturn <- paste0(textToReturn, " This is a low score.")
            } else if (percentage > 90){
                textToReturn <- paste0(textToReturn, " This is a high score. But keep in mind that additional calculation is necessary to determine a real life-time risk. For example having a very high genetic score for something that is not very heritable may make very little difference. These additional calculations typically require further studies, not always available.")
            } else {
                textToReturn <- paste0(textToReturn, " This is a fairly average score.")
            }
        }
        textToReturn <- paste0(textToReturn, " Result from the analysis of ", nrow(SNPs_to_analyze)," SNPs from <u><a target='_blank' href='http://", link, "'>", author," et al (PMID ",pmid,")</a></u>, which were reported to be associated with ", tolower(trait),".")
        textToReturn <- paste0(textToReturn, " This study reports a total sample size of ",sampleSize,", as entered on date ", publication_date,".")

        output[[study_id]] <- list()
        output[[study_id]][["GRS"]] <- GRS_beta
        output[[study_id]][["trait"]] <- tolower(trait)
        output[[study_id]][["percentage"]] <- percentage
        output[[study_id]][["message"]] <- textToReturn
    }
    return(output)
}
