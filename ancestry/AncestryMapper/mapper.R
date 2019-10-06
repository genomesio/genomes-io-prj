pkgs <- c("reshape2", "AncestryMapper")
newPkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(newPkgs)) install.packages(newPkgs)
lapply(pkgs, library, character.only = TRUE)

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
    stop("At least one argument must be given (input folder)\n", call. = FALSE)
} else if (length(args) == 1) {
    # default output file
    args[2] = "mapperOut"
}

# Path to folder containing population references
Refs <- system.file('data', package = 'AncestryMapper')

# Path to CorPheno file
Corpheno <- system.file('extdata', 'CorPheno', package = 'AncestryMapper')

# Path to dbSNP allele data file
# NOTE: changed MinMaxFreq.rda to modified MinMaxFreq.rds
All00Frq <- system.file('data', 'MinMaxFreq.rds', package = 'AncestryMapper')

# function to calculate mean region distances
regionDist <- function(distanceDf, cophenoFile) {
	# convert to long format
	distanceDfMelted <- melt(distanceDf, id.vars = "UNIQID")
	colnames(distanceDfMelted) <- c("id", "Pheno_Pop", "value")
	cList <- distanceDfMelted[grepl("C_", distanceDfMelted$Pheno_Pop),]
	# remove "C_" out of population names
	cList$Pheno_Pop <- gsub("C_", "", cList$Pheno_Pop)
	# read corpheno file to get population info
	popDf <- read.csv(
		cophenoFile, header = TRUE, stringsAsFactors = FALSE, sep = "\t"
	)
	# get populations and their regions
	filteredPopDf <- unique(popDf[,c("Pheno_Pop", "Pheno_Region")])
	# add regions to distanceDf
	mergeDf <- merge(cList, filteredPopDf, by = "Pheno_Pop", all.x = TRUE)
	# rename NA region if present
	if(NA %in% mergeDf$Pheno_Region)
		mergeDf[is.na(mergeDf$Pheno_Region),]$Pheno_Region <- "NaN"
	# calculate region percentage
	regionDf <- aggregate(
		mergeDf$value,
		by = list("region" = mergeDf$Pheno_Region),
		mean # consider using median
	)
	colnames(regionDf) <- c("region", "dist")
	# return
	return(regionDf)
}

############### MAIN ###############

# calculate genetic distance
genetic.distance <- calculateAMidsArith(
	pathTotpeds = args[1],
	NameOut = args[2],
	pathToAriMedoids = Refs,
	pathAll00 = All00Frq
)
# # or read existing distance file
# genetic.distance <- read.csv(
# 	"AMidtped_ref160_inds1_SNPs230308.amid",
# 	header = TRUE,
# 	sep = " ",
# 	stringsAsFactors = FALSE
# )

# plot
pdf(paste0(args[2], ".pdf"), 8, 3)
plotAMids(AMids = genetic.distance, phenoFile = Corpheno, columnPlot = 'I')
dev.off()

# get region distribution
region <- regionDist(genetic.distance, Corpheno)
outDf <- region[order(region$dist),]
write.table(
    region[order(region$dist),],
    file = paste0(args[2], ".txt"),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
)

msg <- paste0("NOTE: remember to change MinMaxFreq.rds back to MinMaxFreq.rda ",
              "when AncestryMapper has a new release")
message(msg)