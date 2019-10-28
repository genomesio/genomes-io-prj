pkgs <- c("reshape2", "stringr", "AncestryMapper")
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
All00Frq <- system.file('data', 'MinMaxFreq.rda', package = 'AncestryMapper')

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
	colnames(regionDf) <- c("region", "distance")
	# return
	regionDf <- regionDf[order(regionDf$distance),]
	return(regionDf)
}

popDist <- function(distanceDf) {
    # convert to long format
    distanceDfMelted <- melt(distanceDf, id.vars = "UNIQID")
    colnames(distanceDfMelted) <- c("id", "Pheno_Pop", "value")
    cList <- distanceDfMelted[grepl("C_", distanceDfMelted$Pheno_Pop),]
    # remove "C_" out of population names
    cList$Pheno_Pop <- gsub("C_", "", cList$Pheno_Pop)
    # split pop name and source
    cList$Pheno_Pop <- gsub("1000.Genomes", "1KG", cList$Pheno_Pop)
    cList$popName <- strReverse(str_split_fixed(strReverse(cList$Pheno_Pop), "\\.", 2)[,2])
    cList$popSource <- strReverse(str_split_fixed(strReverse(cList$Pheno_Pop), "\\.", 2)[,1])
    # return
    colnames(cList) <- c("id", "Pheno_Pop", "distance", "ethnic_group", "source")
    cList <- cList[order(cList$distance),]
    return(cList[,c("ethnic_group", "distance", "source")])
}

strReverse <- function(x) {
    sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")
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
# 	"/Volumes/External/work/genomesio/ancestry/AMidExample_ref143_inds567_SNPs1000.amid",
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
write.table(
    region,
    file = paste0(args[2], "_region_dist.txt"),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
)

# get population distribution
pop <- popDist(genetic.distance)
write.table(
    pop,
    file = paste0(args[2], "_ethnic_dist.txt"),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
)

msg <- paste0("NOTE: remember to change MinMaxFreq.rds back to MinMaxFreq.rda ",
              "when AncestryMapper has a new release")
message(msg)
