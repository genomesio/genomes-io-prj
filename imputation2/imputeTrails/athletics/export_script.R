export_function <- function(uniqueID, moduleDir, outputDir, gtool) {
  if (!file.exists(outputDir)) {
    stop(paste("Did not find a output data with this id", uniqueID))
  }
  
  table_file <- paste0(moduleDir, "/athletics/SNPs_to_analyze.txt")
  request <- table <- read.table(table_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "", quote = "")
  
  # get genotypes
  request <- request[!duplicated(request[,"SNP"]),]
  rownames(request) <- request[,"SNP"]
  genotypes <- get_genotypes(uniqueID = uniqueID, request = request, gtool = gtool, destinationDir = outputDir)
  table[,"genotype"] <- genotypes[table[,"SNP"],]
  
  if ("Beta" %in% colnames(table)) {
    colnames(table)[colnames(table) %in% "Beta"] <- "effect_size"
  }
  
  output <- list()
  
  ### athletics snp list
  output[["athletics_snp"]] < list()
  athleticSnp <- table[table[,"Domain"] %in% "Table1",]
  athleticSnp <- athleticSnp[, c("SNP", "genotype", "Comment")]
  colnames(athleticSnp)<-c("SNP", "Your genotype", "Description")

  for (i in athleticSnp$SNP) {
    output[["athletics_snp"]][[i]] <- as.list(athleticSnp[athleticSnp$SNP == i,])
  }
  
  
  ### injuries and dietary
  output[["injuries_and_dietary"]] <- list()
  output[["injuries_and_dietary"]][["note"]] <- paste0(
    "This table calculates genetic risk scores for all domains covered in Goodlin et al (https://www.ncbi.nlm.nih.gov/pubmed/25919592, ",
    "which covers a number of often encountered injuries and dietary needs in athletics. The risk score is indicated as percentile, ",
    "i.e. 'how many percent of people have a lower score'. So it should not be translated as the direkt risk probability. ",
    "It is just a measure of how you scale relative to other people, based on a measurement of the known genetic component."
  )
  
  domains <- c(
    'ACL rupture', 'Achilles tendon', 'Stress fracture', 'Osteoarthritis', 'Iron Biomarker', 
    'Vitamin E', 'Vitamin D', 'Magnesium', 'Vitamin B', 'Phytosterols', 'Bone mineral density'
  )
  domains <- data.frame(row.names = domains, Domain = domains)
  
  for (i in 1:nrow(domains)) {
    d <- table[table[,"Domain"] %in% rownames(domains)[i],]
    rownames(d) <- d[, "SNP"]

    # try to see if it is numeric-ok
    ef <- suppressWarnings(as.numeric(d[, "effect_size"]))
    if (sum(is.na(ef)) == 0) {
      d[,"effect_size"] <- ef
    } else { #else just insert 1 because it is "effect allele increase risk"
      d[,"effect_size"] <- 1
    }

    d <- suppressWarnings(get_GRS_2(d, mean_scale = TRUE, unit_variance = TRUE, verbose = TRUE))
    population_sum_sd <- sqrt(sum(d[, "population_score_sd"]^2, na.rm = TRUE))
    GRS_beta <- sum(d[,"score_diff"], na.rm = T) / population_sum_sd
    domains[i,"Number of SNPs"] <- sum(!is.na(d[, "score_diff"]))
    domains[i,"Level-score"] <- paste(signif(pnorm(GRS_beta, mean = 0, sd = 1)*100, 3), "%")

    d[,"snps_line"] <- rownames(d)
    d[,"duplicated"] <- duplicated(d[, "Comment"])
    duplicated_snps <- unique(rownames(d)[d[, "duplicated"]])
    for (duplicated_snp in duplicated_snps) {
      w <- which(d[, "Comment"]  %in% d[duplicated_snp, "Comment"])
      d[!d[w,"duplicated"], "snps_line"] <- paste(sort(rownames(d)[w]), collapse = ", ")
    }
    d <- d[!d[,"duplicated"],]
    domains[i, "Source notes"] <- paste(paste0(d[, "snps_line"], ": ", d[,"Comment"]), collapse = "; ")
  }
  
  for(i in rownames(domains)) {
    output[["injuries_and_dietary"]][[i]] <- as.list(domains[i,])
  }

  return(output)
}


