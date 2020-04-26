export_function <- function(uniqueID, moduleDir, outputDir, gtool) {
    if (!file.exists(outputDir)) {
        stop(paste("Did not find a output data with this id", uniqueID))
    }
    
    load(paste0(moduleDir, "/nonsenser/2017-04-05_all_coding_SNPs.rdata"))
    coding_snps[,"SNP"] <- rownames(coding_snps)
    # get genotypes
    genotypes <- get_genotypes(uniqueID = uniqueID, request = coding_snps, gtool = gtool, destinationDir = outputDir, namingLabel = "cached.nonsenser")
    coding_snps[,"Your genotype"] <- genotypes[rownames(coding_snps),]
    
    # get ethnicity parameter
    pDataFile <- paste(outputDir, "/pData.txt", sep = "")
    pData <- try(read.table(pDataFile, header = TRUE, stringsAsFactors = FALSE, sep = "\t"))
    if (class(pData) != "try-error" && "ethnicity" %in% colnames(pData)) {
        ethnicity <- pData[1, "ethnicity"]
    } else {
        ethnicity <- "global"
    }
    
    # get Frequency and 'Common allele' based on 1kg data
    coding_snps[,"Common allele"] <- coding_snps[,"Frequency"] <- NULL #remove the old annovar derived frequency
    if (ethnicity == "global"){
        coding_snps[,"new_freq"] <- coding_snps[, paste0("AF")]
    } else {
        coding_snps[,"new_freq"] <- coding_snps[, paste0(ethnicity,"_AF")]
    }
    
    flips <- which(coding_snps[,"new_freq"] > 0.5)
    coding_snps[, "Common allele"] <- coding_snps[, "REF"]
    coding_snps[, "Minor allele"] <- coding_snps[, "ALT"]
    coding_snps[, "Frequency"] <- coding_snps[, "new_freq"] 
    coding_snps[flips, "Common allele"] <- coding_snps[flips, "ALT"]
    coding_snps[flips, "Minor allele"] <- coding_snps[flips, "REF"]
    coding_snps[flips, "Frequency"] <- 1 - coding_snps[flips, "new_freq"] 
    
    out <- coding_snps[, c("SNP", "Your genotype", "Common allele", "Frequency", "SIFT_pred", "MutationTaster_pred", "Polyphen2_HDIV_pred")]
    colnames(out)[5] <- "SIFT"
    colnames(out)[6] <- "MutationTaster"
    colnames(out)[7] <- "PolyPhen2"
    
    out[,"Heterozygote"] <- sapply(strsplit(out[,"Your genotype"],"/"), function(x) {x[1] != x[2]})
    out[,"Unmeasured"] <- is.na(out[,"Your genotype"])
    out[,"Homozygote allele"] <- sub("/.$","", out[,"Your genotype"])
    out[out[,"Heterozygote"] & !out[,"Unmeasured"], "Homozygote allele"] <- NA
    out[,"Homozygote minor"] <- out[,"Homozygote allele"] != out[,"Common allele"] & out[,"Common allele"] != ""
    out[is.na(out[,"Homozygote minor"]),"Homozygote minor"] <- FALSE
    type <- rep("Homozygote major?", nrow(out))
    type[out[,"Homozygote minor"]] <- "Homozygote minor?"
    type[out[,"Heterozygote"]] <- "Heterozygote"
    type[out[,"Unmeasured"]] <- "Unmeasured"
    out[,"Type"] <- factor(type, levels = c("Homozygote minor?", "Heterozygote", "Homozygote major?", "Unmeasured"))
    out <- out[order(out[, "Type"], out[, "Frequency"]),]
    out[,"Heterozygote"] <- out[,"Unmeasured"] <- out[,"Homozygote allele"] <- out[,"Homozygote minor"]<-NULL
    
    # return
    output <- list()
    output[["note"]] <- paste(
        "Most SNPs in the genome are not actually found within a gene: They are 'intergenic'. ",
        "When talking about a gene-mutation however, as is done in popular media, most often the meaning is a SNP that alters the sequence of a gene.",
        "Because of selection pressure throughout our evolution, these are rare.",
        "Also, they are often the focus of scientific studies using DNA-sequencing technology to discover causes of rare diseases.",
        "However, interestingly many of us actually have these 'gene-breaking' SNPs while nonetheless being perfectly healthy.",
        "The imputation technology used with this site, gives the opportunity to identify a number of these based on just on genotyping",
        "microarray results. If you give your ID-code to this module a table of all measured missense and nonsense mutations will be presented.",
        "Interpretation of the table can be done in many ways and unlike other modules, this does not give 'one true answer'.",
        "One method is to search for SNPs where you have one or two copies of the non-common allele and then investigate the consequence",
        "using other resources such as http://www.ncbi.nlm.nih.gov/SNP/ or http://exac.broadinstitute.org/.",
        "Note however that the definition of 'common' is very dependent on ethnicity: in this browser common just means the allele most often found in impute.me-users.",
        "However, it is recommended to check the ethnical distribution in e.g. the http://www.1000genomes.org/.",
        "Another help provided is the polyphen-scores (http://genetics.bwh.harvard.edu/pph2/) and SIFT-scores (http://sift.jcvi.org/),",
        "which can give an indication of the consequence. Ultimately the goal of this is to satisfy ones curiousity about the state of your functional genes.",
        "By being healthy, in spite of a specific broken gene, you'd be contributing to complete our view of genes and how they work."
    )
    for (i in out$SNP) {
        output[[i]] <- as.list(out[out$SNP == i,])
    }
    return(output)
}
