export_function <- function (uniqueID, moduleDir, outputDir, functionFile, gtool) {
    if (!file.exists(outputDir)) {
        stop(paste("Did not find a output data with this id", uniqueID))
    }

    output <- list()
    #############
    # First get height
    #############

    # get gender
    pDataFile <- paste(outputDir, "/pData.txt", sep = "")
    gender <- read.table(pDataFile, header = TRUE, stringsAsFactors = FALSE, sep = "\t")[1, "gender"]

    # get the current best predictor SNPs
    giant_sup_path <- paste0(moduleDir, "/guessMyHeight/SNPs_to_analyze.txt")
    giant_sup <- read.table(giant_sup_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

    # get genotypes and calculate gheight
    genotypes <- get_genotypes(uniqueID = uniqueID, request = giant_sup, gtool = gtool, destinationDir = outputDir)
    gheight <- get_GRS(genotypes = genotypes, betas = giant_sup)

    # calculate auxilary info
    gheight_snp_count <- paste(sum(!is.na(genotypes[, "genotype"])), "of", nrow(genotypes))

    # calculate estimated_real height (note this one may need a lot of fine tuning later)
    women_mean <- 1.62
    women_sd <- 0.0796
    men_mean <- 1.76
    men_sd <- 0.0802
    if (gender == 1) {
        mean_here <- men_mean
        sd_here <- men_sd
    } else {
        mean_here <- women_mean
        sd_here <- women_sd
    }
    gheight_m_estimate <- mean_here + gheight*sd_here
    gheight_m_estimate_cm <- gheight_m_estimate * 100
    output[["gheight_Z_score"]] <- gheight
    output[["gheight_snp_count"]] <- gheight_snp_count
    output[["gheight_m_estimate"]] <- gheight_m_estimate_cm

    #############
    # Then get colour
    #############

    #get the gColour
    GRS_file_name <- paste0(moduleDir, "/hairColour/SNPs_to_analyze.txt")
    GRS_file <- read.table(GRS_file_name, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    for (component in c("blonde", "red")) {
        s1 <- GRS_file[GRS_file[, "Category"] %in% component, ]
        rownames(s1) <- s1[, "SNP"]
        # get genotypes and calculate gHairColour
        s1[, "genotype"] <- get_genotypes(uniqueID = uniqueID, request = s1, gtool = gtool, destinationDir = outputDir)
        s1 <- get_GRS_2(s1, mean_scale = TRUE, unit_variance = TRUE)
        population_sum_sd <- sqrt(sum(s1[, "population_score_sd"]^2, na.rm = TRUE))
        GRS <- sum(s1[, "score_diff"], na.rm = TRUE) / population_sum_sd
        assign(paste("gColour", component, sep = "_"), GRS)
    }

    blond_calibrate <- function(x) {max(c(0, min(c(1, (x+1)/6))))}
    red_calibrate <- function(x) {max(c(0, min(c(1, (x+1)/5))))}

    blondeness <- blond_calibrate(gColour_blonde)
    redheadness <- red_calibrate(gColour_red)

    # Calculate colour #with the line from the image map
    colour <- hsv(h = 0.1 - (redheadness/10), s = min(c(1, 1 - blondeness + (redheadness/2))), v = blondeness)

    output[["hair_colour"]] <- colour
    return(output)
}
