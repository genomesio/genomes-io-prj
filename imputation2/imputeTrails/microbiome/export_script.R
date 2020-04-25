export_function <- function(uniqueID, moduleDir, outputDir, gtool) {
    if (!file.exists(outputDir)) {
        stop(paste("Did not find a output data with this id", uniqueID))
    }
    
	table_file <- paste0(moduleDir, "/microbiome/SNPs_to_analyze.txt")
	table <- read.table(table_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
	rownames(table) <- table[,"SNP"]
	genotypes <- get_genotypes(uniqueID = uniqueID, request = table, gtool = gtool, destinationDir = outputDir)
	table[,"Your genotype"] <- genotypes[rownames(table),]
	table <- table[,c("SNP", "Your genotype", "Increasing_allele", "Bacteria")]

	# return
	output <- list()
	output[["note"]] <- paste(
	    "For each of the strains of bacteria you can read your own genotype and see",
	    "if your are host-genetically disposed to have an increased proportion of this particular strain.",
	    "The consequence of being particularly disposed are fairly unclear, however,",
	    "but refer to the primary literature on gut microbiome for more information.",
	    "The findings were taken from the study by Blekhman et al (2015) (http://www.ncbi.nlm.nih.gov/pubmed/?term=26374288,",
	    "Supplementary material table S5.)"
	)
	for (i in table$SNP) {
	    output[[i]] <- as.list(table[table$SNP == i,])
	}
	return(output)
}