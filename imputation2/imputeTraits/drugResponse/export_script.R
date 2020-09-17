export_function <- function (uniqueID, moduleDir, outputDir, gtool) {
    if (!file.exists(outputDir)) {
        stop(paste("Did not find a output data with this id", uniqueID))
    }
    suppressWarnings(dir.create(file.path(outputDir, "table_out")))
    suppressWarnings(dir.create(file.path(outputDir, "table_out/drugResponse")))

    table_file <- paste0(moduleDir, "/drugResponse/SNPs_to_analyze.txt")
    SNPs_to_analyze <- read.table(table_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    SNPs_to_analyze[, "PMID"] <- as.character(SNPs_to_analyze[, "PMID"])

    output <- list()
    # output[["documentation"]] <- list()
    # output[["documentation"]][["source"]] <- "https://github.com/lassefolkersen/impute-me/blob/fbbd86ac652c74d59d94b5e01ab5c5600bd8fea7/drugResponse/SNPs_to_analyze.txt"

    # retrieving SNPs
    SNPs_to_retrieve <- SNPs_to_analyze
    SNPs_to_retrieve <- SNPs_to_retrieve[!duplicated(SNPs_to_retrieve[, "SNP"]), ]
    rownames(SNPs_to_retrieve) <- SNPs_to_retrieve[, "SNP"]
    SNPs_to_retrieve <- get_genotypes(uniqueID = uniqueID, request = SNPs_to_retrieve, gtool = gtool, destinationDir = outputDir)
    # inserting SNPs and calculating GRS
    SNPs_to_analyze[, "genotype"] <- SNPs_to_retrieve[SNPs_to_analyze[, "SNP"], "genotype"]
    for (study in unique(SNPs_to_analyze[, "PMID"])){
        d2 <- SNPs_to_analyze[SNPs_to_analyze[, "PMID"] %in% study, ]
        disease <- unique(d2[, "disease"])
        drug <- unique(d2[, "drug"])

        # handling duplicated SNPs
        if(any(duplicated(d2[, "SNP"]))){
            d2 <- d2[!duplicated(d2[, "SNP"]), ]
        }

        rownames(d2) <- d2[, "SNP"]
        d3 <- try(get_GRS_2(d2, mean_scale = TRUE, unit_variance = TRUE, verbose = FALSE))
        if (class(d3) == "try-error") {
            z_score <- "Not calculated"
            percentage <- "Not calculated"

        } else {
            population_sum_sd <- sqrt(sum(d3[, "population_score_sd"]^2, na.rm = TRUE))
            if (population_sum_sd == 0) {
                z_score <- "Not calculated"
                percentage <- "Not calculated"
            } else {
                GRS <- sum(d3[, "score_diff"], na.rm = TRUE) / population_sum_sd
                percentage <- floor(pnorm(GRS, mean = 0, sd = 1) * 100)
                z_score <- signif(GRS, 2)
                percentage <- signif(percentage, 2)
            }
        }

        output[[study]] <- list(
            disease = disease,
            drug = drug,
            GRS = z_score,
            percentage = percentage,
            pop_sd = population_sum_sd
        )
        
        ### save snp table
        table <- SNPs_to_analyze[SNPs_to_analyze[, "PMID"] %in% study, ]
        # join with score_diff from d3
        d3sub <- d3[,c("SNP","score_diff")]
        table <- merge(table, d3sub, by = "SNP", all.x = TRUE)
        table[, "minor/major allele"] <- apply(table[,c("minor_allele", "major_allele")], 1, paste, collapse = "/")
        table[, "effect/alternative allele"] <- apply(table[,c("effect_allele", "non_effect_allele")], 1, paste, collapse = "/")
        order <- c("SNP", "genotype" , "effect/alternative allele", "effect_size", "effect_direction", "effect_measure", "score_diff", "minor/major allele", "minor_allele_freq", "gene", "PMID")
        missing <- order[!order %in% colnames(table)]
        if (length(missing) > 0) stop(safeError("Missing some columns"))
        table[,"minor_allele_freq"] <- signif(table[,"minor_allele_freq"], 2)
        table[,"effect_size"] <- signif(table[,"effect_size"], 2)
        table <- table[,order]
        colnames(table) <- c(
            "SNP","Your Genotype","Risk/ non-risk Allele","Effect Size","Effect Direction", "Effect Measure",
            "SNP-score (population normalized)","Minor/ major Allele","Minor Allele Frequency","Reported Gene","PMID"
        )
        write.table(table, paste0(outputDir, "/table_out/drugResponse/", study, ".txt"), 
                    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
    return(output)
}
