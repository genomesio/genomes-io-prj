#' Do imputation
run_imputation <- function (
	runDir,
	shapeit = NULL,
	plink = NULL,
	impute2 = NULL,
	imputeDataDir = NULL,
	sample_ref = NULL
) {
	# check working directory
	if (class(runDir) != "character") stop(paste("runDir must be character, not", class(runDir)))
	if (length(runDir) != 1) stop(paste("runDir must be lengh 1, not", length(runDir)))
	if (!file.exists(runDir)) stop(paste("Did not find runDir at path:", runDir))
	if (length(grep("/$", runDir)) != 0) stop("Please don't use a trailing slash in the runDir")

	# load and check variables.rdata and raw_data
	# load(paste(runDir, "/variables.rdata", sep=""))
	rawdata <- paste(runDir, "/", uniqueID, "_raw_data.txt", sep="")

	if (class(rawdata) != "character") stop(paste("rawdata must be character, not", class(rawdata)))
	if (length(rawdata) != 1) stop(paste("rawdata must be lengh 1, not", length(rawdata)))
	if (!file.exists(rawdata)) stop(paste("Did not find rawdata at path:", rawdata))

	# check paths to tools
	if (class(shapeit) != "character") stop(paste("shapeit must be character, not", class(shapeit)))
	if (length(shapeit) != 1) stop(paste("shapeit must be lengh 1, not", length(shapeit)))
	if (!file.exists(shapeit)) stop(paste("Did not find shapeit at path:", shapeit))

	if (class(plink) != "character") stop(paste("plink must be character, not", class(plink)))
	if (length(plink) != 1) stop(paste("plink must be lengh 1, not", length(plink)))
	if (!file.exists(plink)) stop(paste("Did not find plink at path:", plink))

	if (class(impute2) != "character") stop(paste("impute2 must be character, not", class(impute2)))
	if (length(impute2) != 1) stop(paste("impute2 must be lengh 1, not", length(impute2)))
	if (!file.exists(impute2)) stop(paste("Did not find impute2 at path:", impute2))

	if (class(imputeDataDir) != "character") stop(paste("imputeDataDir must be character, not", class(imputeDataDir)))
	if (length(imputeDataDir) != 1) stop(paste("imputeDataDir must be lengh 1, not", length(imputeDataDir)))
	if (!file.exists(imputeDataDir)) stop(paste("Did not find imputeDataDir at path:", imputeDataDir))

	if (class(sample_ref) != "character") stop(paste("sample_ref must be character, not", class(sample_ref)))
	if (length(sample_ref) != 1) stop(paste("sample_ref must be lengh 1, not", length(sample_ref)))
	if (!file.exists(sample_ref)) stop(paste("Did not find sample_ref at path:", sample_ref))

	# need to always check if the genes_for_good_cleaner should be run
	if (length(grep("genes for good", tolower(readLines(rawdata, n = 5))) > 0)) {
	    genes_for_good_cleaner(rawdata)
	}

	# print start message
	cat(paste0(Sys.time(), "\nStarting imputation\n"))
	setwd(runDir)

	# load data using plink and convert into plink format
	cmd1 <- paste0(plink, " --23file ", rawdata, " ", uniqueID, " ", uniqueID, " --recode --out step_1")
	out1 <- system(cmd1)

	# omit duplicates
	map <- read.table('step_1.map', sep = '\t', stringsAsFactors = FALSE, comment.char = "")
	exclude <- map[duplicated(map[, 4]), 2]
	print(paste('Removed', length(exclude), 'SNPs that were duplicated'))
	write.table(exclude, file = 'step_2_exclusions', sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

	# loop over chromosomes
	for (chr in c("X", as.character(1:22))) {
		# extract only one specific chromosome
		cmd2 <- paste(plink, " --file step_1 --chr ", chr, " --recode --out step_2_chr", chr, " --exclude step_2_exclusions", sep = "")
		out2 <- system(cmd2)

		# if X chromosome is missing it is allowed to skip forward
		if (out2 == 13 & chr == "X") {
			print("Didn't find X-chr data, so skipping that")
			next
		}

		# then check for strand flips etc.
		cmd3 <- paste(
			shapeit,
			" -check --input-ped step_2_chr", chr, ".ped step_2_chr", chr, ".map",
			" -M ", imputeDataDir, "/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr", chr, "_combined_b37.txt",
			" --input-ref ", imputeDataDir, "/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr", chr, "_impute.hap.gz ", imputeDataDir, "/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr", chr, "_impute.legend.gz ", sample_ref,
			" --output-log step_2_chr", chr, "_shapeit_log",
			sep = ""
		)
		system(cmd3)

		# many homozygote SNPs will fail the check, because they don't have the ref-allele. So we need to sort them
		logFile <- read.table(paste("step_2_chr", chr, "_shapeit_log.snp.strand", sep = ""), sep = '\t', stringsAsFactors = FALSE, header = FALSE, skip = 1)
		omitMissing <- logFile[logFile[, 1] %in% 'Missing', 3]
		logStrand <- logFile[logFile[, 1] %in% 'Strand', ]
		omitNonIdentical <- logStrand[logStrand[, 5] != logStrand[, 6], 3]
		omitBlank <- logStrand[logStrand[, 5]%in%'', 3]

		# create another (fake) person with the alternative allele just for their sake.
		# this next command takes all the homozygotes, minus the indels (which are too complicated to lift out from 23andme)
		forceHomozygoteTable <- logStrand[
			logStrand[, 5] == logStrand[, 6] &
				nchar(logStrand[, 9]) == 1 &
				nchar(logStrand[, 10]) == 1 &
				!logStrand[, 5] %in% c("D", "I") &
				!logStrand[, 6] %in% c("D", "I"),
			]

		# remove any cases where there are more than two alleles involved
		forceHomozygoteTable <- forceHomozygoteTable[sapply(apply(forceHomozygoteTable[, c(5, 6, 9, 10)], 1, unique), length) == 2, ]

		# remove any duplicates there might be
		forceHomozygoteTable <- forceHomozygoteTable[!duplicated(forceHomozygoteTable[, 4]), ]
		map <- read.table(paste("step_2_chr", chr, ".map", sep = ""), sep = "\t", stringsAsFactors = FALSE, comment.char = "")

		# load the ped file, and doubles it
		ped2 <- ped1 <- strsplit(readLines(paste("step_2_chr", chr, ".ped", sep="")), " ")[[1]]
		ped2[1] <- "Temporary"
		ped2[2] <- "Non_person"
		if ((length(ped1)-6) / 2 != nrow(map)) stop("mismatch between map and ped")
		replacementPos <- which(map[, 2] %in% forceHomozygoteTable[, 4])
		A1_pos <- 7 + 2 * (replacementPos - 1)
		A2_pos <- 8 + 2 * (replacementPos - 1)
		ped2[A1_pos] <- forceHomozygoteTable[, 9]
		ped2[A2_pos] <- forceHomozygoteTable[, 10]
		ped <- rbind(ped1, ped2)
		write.table(ped, paste("step_3_chr", chr, ".ped", sep=""), sep=" ", col.names=FALSE, row.names=FALSE, quote=F)

		omitRemaining <- logStrand[!logStrand[, 4] %in% forceHomozygoteTable[, 4], 3]
		print(paste('Omitting', length(omitMissing), 'because of missing', length(omitBlank), 'because they are blank, and', length(omitNonIdentical), 'true strand flips'))
		write.table(
			c(omitNonIdentical, omitBlank, omitMissing, omitRemaining),
			file = paste("step_3_chr", chr, "_exclusions", sep = ""),
			sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE
		)

		# run the shapeit command (with two people, the right one and a placeholder heterozygote)
		cmd4 <- paste(
			shapeit,
			" --input-ped step_3_chr", chr, ".ped step_2_chr", chr, ".map",
			" -M ", imputeDataDir, "/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr", chr, "_combined_b37.txt",
			" --input-ref ", imputeDataDir, "/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr", chr, "_impute.hap.gz ", imputeDataDir, "/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr", chr, "_impute.legend.gz ", sample_ref,
			" --output-log step_4_chr", chr, "_shapeit_log",
			" --exclude-snp step_3_chr", chr, "_exclusions -O step_4_chr", chr,
			sep = ""
		)
		system(cmd4)

		# check for errors and stopping if there are any. No point to continue otherwise
		log <- readLines(paste("step_4_chr", chr, "_shapeit_log.log", sep=""))
		if (substr(log[length(log)], 1, 5) == "ERROR") {
			stop(paste("At chr", chr, " the shapeit failed. Check this file for explanation: step_4_chr", chr, "_shapeit.log", sep = ""))
		}

		# remove the placeholder person
		cmd5_1 <- paste("cut --delimiter=\" \" -f 1-7 step_4_chr", chr, ".haps > step_5_chr", chr, ".haps", sep = "")
		system(cmd5_1)
		cmd5_2 <- paste("head -n 3 step_4_chr", chr, ".sample > step_5_chr", chr, ".sample", sep = "")
		system(cmd5_2)

		# detect max length of each chromosome
		cmd6 <- paste("zcat ", imputeDataDir, "/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr", chr, "_impute.legend.gz | tail -n 1 | cut --delimiter=\" \" -f 2", sep = "")
		print(cmd6)
		maxPos <- as.numeric(system(cmd6, intern = TRUE))

		# iterate over 5e6 chunks
		# chr[X]_nonPAR_combined_b37.txt, _nonPAR_impute.hap.gz and _nonPAR_impute.legend.gz
		starts <- seq(0, maxPos, 5e6)
		for (i in 1:length(starts)) {
			start <- starts[i]
			end <- start + 5e6

			cmd7 <- paste(
				impute2,
				" -m ", imputeDataDir, "/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr", chr, "_combined_b37.txt",
				" -h ", imputeDataDir, "/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr", chr, "_impute.hap.gz",
				" -l ", imputeDataDir, "/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr", chr, "_impute.legend.gz -known_haps_g step_5_chr", chr, ".haps",
				" -int ", start, " ", end, " -Ne 20000 -o step_7_chr", chr, "_", i,
				sep = ""
			)
			step_7_log <- system(cmd7)

			# test for memory-lack bug (step_7_log will be 137 if killed, otherwise 0)
			if (step_7_log == 137) {
				# we divide the job in smaller bits
				divisions <- 3
				for (j in 1:divisions) {
					start_2 <- floor(starts[i] + (j-1)*(5e6/ divisions))
					end_2 <- floor(starts[i]+ (j)*(5e6/ divisions))
					print(paste("restart imputation with new subset to avoid memory-lack bug:", start_2, "to", end_2)  )

					cmd7<-paste(
						impute2,
						" -m ", imputeDataDir, "/ALL_1000G_phase1integrated_v3_impute/genetic_map_chr", chr, "_combined_b37.txt",
						" -h ", imputeDataDir, "/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr", chr, "_impute.hap.gz",
						" -l ", imputeDataDir, "/ALL_1000G_phase1integrated_v3_impute/ALL_1000G_phase1integrated_v3_chr", chr, "_impute.legend.gz -known_haps_g step_5_chr", chr, ".haps",
						" -int ", start_2, " ", end_2, " -Ne 20000 -o step_7_chr", chr, "_", i, "-", j,
						sep=""
					)
					step_7_log_2<-system(cmd7)
					if (step_7_log_2 == 137) print("the memory problem was still active after second round")
				}
			}
		}
	}
}

#' Summarize imputation
summarize_imputation<-function(
	runDir = NULL,
	uniqueID = NULL,
	destinationDir = NULL,
	gtool = NULL,
	plink = NULL
) {
	# check paths
	if (class(runDir) != "character") stop(paste("runDir must be character, not", class(runDir)))
	if (length(runDir) != 1) stop(paste("runDir must be lengh 1, not", length(runDir)))
	if (!file.exists(runDir)) stop(paste("Did not find runDir at path:", runDir))
	if (length(grep("/$", runDir)) != 0) stop("Please don't use a trailing slash in the runDir")
	setwd(runDir)

	if (class(uniqueID) != "character") stop(paste("uniqueID must be character, not", class(uniqueID)))
	if (length(uniqueID) != 1) stop(paste("uniqueID must be lengh 1, not", length(uniqueID)))

	if (class(destinationDir) != "character") stop(paste("destinationDir must be character, not", class(destinationDir)))
	if (length(destinationDir) != 1) stop(paste("destinationDir must be lengh 1, not", length(destinationDir)))
	if (!file.exists(destinationDir)) stop(paste("Did not find destinationDir at path:", destinationDir))
	if (length(grep("/$", destinationDir)) != 0) stop("Please don't use a trailing slash in the destinationDir")

	if (class(gtool) != "character") stop(paste("gtool must be character, not", class(gtool)))
	if (length(gtool) != 1) stop(paste("gtool must be lengh 1, not", length(gtool)))
	if (!file.exists(gtool)) stop(paste("Did not find gtool at path:", gtool))

	if (class(plink) != "character") stop(paste("plink must be character, not", class(plink)))
	if (length(plink) != 1) stop(paste("plink must be lengh 1, not", length(plink)))
	if (!file.exists(plink)) stop(paste("Did not find plink at path:", plink))

	if (file.exists(destinationDir)) {
		if (length(list.files(destinationDir)) > 0) {
			stop(paste0("The destinationDir '", destinationDir, "' already exists and has files in it. This is a major unforeseen error") )
		} else {
			dir.create(paste0(destinationDir))
		}
	}

	# print start message
	cat(paste0(Sys.time(), "\nStarting summarize imputation results for: ", uniqueID, "\n"))

	# get imputation files
	allFiles<-list.files(runDir)
	step7Files<-grep("^step_7_chr", allFiles, value = TRUE)
	step7ResultsFiles<-grep("[0-9]$", step7Files, value = TRUE)
	chromosomes<-unique(sub("_[0-9-]+$", "", sub("^step_7_chr", "", step7ResultsFiles)))
	chromosomes<-chromosomes[order(suppressWarnings(as.numeric(chromosomes)))]

	# get the right order, not this is somewhat complicated by the fact that most chunks are marked numerically, then some are marked as e.g. 10-1, 10-2 (see the memory limit issue)
	for (chr in chromosomes) {
		print(paste("Merging chunks in chromosome", chr))
		s <-grep(paste("^step_7_chr", chr, "_", sep = ""), step7ResultsFiles, value = TRUE)
		main_number <- as.numeric(sub("-[0-9]", "", sub("^.+_", "", s)))
		suffix_number <- suppressWarnings(as.numeric(sub("^-", "", sub("^[0-9]+", "", sub("^.+_", "", s)))))
		s <- s[order(main_number, suffix_number)]
		print(paste("For chr", chr, "these were the files to merge:", paste(s, collapse = ", ")))
		cmd1 <- paste("cat ", paste(s, collapse=" "), " > ", uniqueID, "_chr", chr, ".gen", sep = "")
		system(cmd1)
		unlink(s)
	}

	genFiles<-paste(uniqueID, "_chr", chromosomes, ".gen", sep = "")
	if (length(genFiles)==0) stop("Didn't find a single gen-file")

	# run a file conversion, first to plink then to 23andme/simple-format
	for (genFile in genFiles) {
		chr <- sub("\\.gen$", "", sub("^.+_chr", "", genFile))
		print(paste("Simplifying in chromosome", chr))
		sampleFile <- paste("step_4_chr", chr, ".sample", sep = "")

		# make list of indels
		cmd2<-paste(
			"awk -F' ' '{ if ((length($4) > 1 ) || (length($5) > 1 ) || $4 == \"-\" || $5 == \"-\") print $2 }'",
			genFile, ">",
			paste("step_8_chr", chr, "_snps_to_exclude", sep = "")
		)
		system(cmd2)

		# exclude indels
		cmd3 <- paste(
			gtool,
			" -S --g ", genFile,
			" --s ", sampleFile,
			" --exclusion step_8_chr", chr, "_snps_to_exclude",
			" --og step_8_chr", chr, ".gen",
			sep = ""
		)
		system(cmd3)

		# convert to ped format
		cmd4 <- paste(
			gtool,
			" -G --g step_8_chr", chr, ".gen",
			" --s ", sampleFile,
			" --chr ", chr,
			" --snp",
			sep = ""
		)
		system(cmd4)

		# convert to 23andme format
		cmd5 <- paste(plink, " --file step_8_chr", chr, ".gen --recode 23 --out step_9_chr", chr, sep = "")
		system(cmd5)

		# remove comment lines
		cmd6 <- paste("less step_9_chr", chr, ".txt | grep -v \"#\" > step_10_chr", chr, ".txt", sep = "")
		system(cmd6)

		# the next step 8 and 9 sometime fails for no apparent reason. Probably memory. We therefore make a checkup, where
		# it is checked if the file actually exists and if not - a more complicated step splits it up in chunks.
		fileExists <- file.exists(paste("step_10_chr", chr, ".txt", sep = ""))
		if (fileExists) {
			size <- file.info(paste("step_10_chr", chr, ".txt", sep = ""))["size"]
		} else {
			size <- 0
		}

		# re-run if it's less than 100 bytes (fair to assume something was wrong then)
		if (size < 100) {
			print(paste("Retrying step 8-9 command for chr", chr, ". Trying to split it in pieces (non-normal low memory running)"))
			cmd7 <- paste("split --verbose --lines 5000000 step_8_chr", chr, ".gen step_8_extra_chr", chr, ".gen", sep = "")
			system(cmd7)
			chunks<-grep(paste("step_8_extra_chr", chr, "\\.gena[a-z]$", sep = ""), list.files(runDir), value = TRUE)
			for (chunk in chunks) {
				ch <- sub("^.+\\.", "", chunk)
				cmd8 <- paste(gtool, " -G --g ", chunk, " --s ", sampleFile, " --chr ", chr, " --snp", sep = "")
				system(cmd8)
				# reform to plink fam/bim/bed file
				cmd9 <- paste(plink, " --file ", chunk, " --recode --transpose --noweb --out step_9_chr", chr, "_", ch, sep = "")
				system(cmd9)
				#re-order to 23andme format
				cmd10 <- paste("awk '{ print $2 \"\t\" $1 \"\t\"$4\"\t\" $5 $6}' step_9_chr", chr, "_", ch, ".tped > step_9_chr", chr, "_split_", ch, ".txt", sep = "")
				system(cmd10)
			}
			cmd11 <- paste("cat ", paste(paste("step_9_chr", chr, "_split_", sub("^.+\\.", "", chunks), ".txt", sep = ""), collapse=" "), " > step_10_chr", chr, ".txt", sep = "")
			system(cmd11)
		}

		# remove NN
		cmd12 <- paste("awk '{ if ($4 != \"NN\") print}' step_10_chr", chr, ".txt > step_11_chr", chr, ".txt", sep = "")
		system(cmd12)

		# remove duplicates
		cmd13 <- paste("awk -F', ' '!seen[$1]++' step_11_chr", chr, ".txt > ", sub("\\.gen$", "", genFile), ".simple_format.txt", sep = "")
		system(cmd13)

		# removing temporary files from step 8, 9 and 10
		unlink(list.files(runDir, pattern=paste0("^step_8_chr", chr), full.names=T))
		unlink(list.files(runDir, pattern=paste0("^step_9_chr", chr), full.names=T))
		unlink(list.files(runDir, pattern=paste0("^step_10_chr", chr), full.names=T))
	}

	##### move summarized files into destination folder
	# prepDestinationDir <- paste(destinationDir, "/", uniqueID, sep = "")
	prepDestinationDir <- destinationDir
	if (!file.exists(prepDestinationDir)) dir.create(prepDestinationDir)

	# zip and move simple_format files to destinationDir
	zipFile_simpleformat <- paste(runDir, paste(uniqueID, ".simple_format.zip", sep = ""), sep = "/")
	twentythreeandmeFiles <- paste(uniqueID, "_chr", chromosomes, ".simple_format.txt", sep = "")
	zip(zipFile_simpleformat, twentythreeandmeFiles, flags = "-r9X", extras = "", zip = Sys.getenv("R_ZIPCMD", "zip"))
	zipFile_simpleformat_destination <- paste(prepDestinationDir, basename(zipFile_simpleformat), sep = "/")
	move_result <- file.rename(zipFile_simpleformat, zipFile_simpleformat_destination)
	if (!move_result) { # this would fail, for example if data is on another volume. But we can still copy
		file.copy(zipFile_simpleformat, zipFile_simpleformat_destination)
		unlink(zipFile_simpleformat)
	}
	unlink(list.files(runDir, pattern = "23andme", full.names = TRUE)) ### VINH: NO EFFECT? SHOULD USE pattern="simple_format" INSTEAD?

	# zip and move gen files
	zipFileGen <- paste(runDir, paste(uniqueID, ".gen.zip", sep = ""), sep = "/")
	zip(zipFileGen, genFiles, flags = "-r9X", extras = "", zip = Sys.getenv("R_ZIPCMD", "zip"))
	zipFileGen_destination <- paste(prepDestinationDir, basename(zipFileGen), sep = "/")
	move_result <- file.rename(zipFileGen, zipFileGen_destination)
	if (!move_result) { #this would fail, for example if data is on another volume. But we can still copy
		file.copy(zipFileGen, zipFileGen_destination)
		unlink(zipFileGen)
	}
	unlink(genFiles)

	# move the original input file
	zipFileOriginal <- paste(runDir, paste(uniqueID, ".input_data.zip", sep = ""), sep = "/")
	zip(zipFileOriginal, paste(uniqueID, "_raw_data.txt", sep = ""), flags = "-r9X", extras = "", zip = Sys.getenv("R_ZIPCMD", "zip"))
	zipFileOriginal_destination <- paste(prepDestinationDir, basename(zipFileOriginal), sep = "/")
	move_result <- file.rename(zipFileOriginal, zipFileOriginal_destination)
	if (!move_result) { #this would fail, for example if data is on another volume. But we can still copy
		file.copy(zipFileOriginal, zipFileOriginal_destination)
		unlink(zipFileOriginal)
	}

	# create the pData file
	# load(paste0(runDir, "/variables.rdata"))
	# if (!exists("should_be_imputed")) should_be_imputed  <- NA
	# if (!exists("imputemany_upload")) imputemany_upload  <- NA
	# if (!exists("upload_time")) upload_time <- NA

	timeStamp <- format(Sys.time(), "%Y-%m-%d-%H-%M")
	# md5sum <- md5sum(paste(uniqueID, "_raw_data.txt", sep = ""))
	gender <- system(paste("cut --delimiter=\" \" -f 6 ", runDir, "/step_4_chr22.sample", sep = ""), intern = TRUE)[3]

	f <- file(paste0(prepDestinationDir, "/pData.txt"), "w")
	# writeLines(paste(c("uniqueID", "filename", "email", "first_timeStamp", "md5sum", "gender", "protect_from_deletion", "should_be_imputed", "imputemany_upload", "upload_time"), collapse="\t"), f)
	# writeLines(paste(c(uniqueID, file, email, timeStamp, md5sum, gender, protect_from_deletion, should_be_imputed, imputemany_upload, upload_time), collapse="\t"), f)
	writeLines(paste(c("uniqueID", "first_timeStamp", "gender"), collapse = "\t"), f)
	writeLines(paste(c(uniqueID, timeStamp, gender), collapse = "\t"), f)
	close(f)
}

#' get genotypes
get_genotypes <- function (
	uniqueID = NULL,
	request = NULL,
	gtool = NULL,
	destinationDir = NULL,
	namingLabel = "cached",
	call_threshold = 0.8 # threshold for calling SNP. Ok with 0.8 for multi-SNP signatures, but should definetly be increased in special high-importance SNPs. Default from gtool is suggested at 0.9.
) {
	# check namingLabel
	if (class(namingLabel) != "character") stop(paste("namingLabel must be character, not", class(namingLabel)))
	if (length(namingLabel) != 1) stop(paste("namingLabel must be lengh 1, not", length(namingLabel)))

	# check data in uniqueID's data folder
	if (class(uniqueID) != "character") stop(paste("uniqueID must be character, not", class(uniqueID)))
	if (length(uniqueID) != 1) stop(paste("uniqueID must be lengh 1, not", length(uniqueID)))
	idFolder <- destinationDir
	if (!file.exists(idFolder)) stop(paste("Did not find an output folder at", idFolder))
	genZipFile <- paste(idFolder, "/", uniqueID, ".gen.zip", sep = "")
	inputZipFile <- paste(idFolder, "/", uniqueID, ".input_data.zip", sep = "")
	cachedGenotypeFile <- paste(idFolder, "/", uniqueID, ".", namingLabel, ".gz", sep = "")
	if (!file.exists(cachedGenotypeFile)) print(paste0("Did not find a '", namingLabel, "' chachedGenotypeFile file in idFolder at '", idFolder, "' but that's no problem"))

	# create a temp folder to use
	idTempFolder <- paste(destinationDir, "temp", sep = "/")
	if (file.exists(idTempFolder)) stop(safeError(paste("Temp folder exists, this could indicate that", uniqueID, "is already worked on!")))

	# check other variables
	if (class(gtool) != "character") stop(paste("gtool must be character, not", class(gtool)))
	if (length(gtool) != 1) stop(paste("gtool must be lengh 1, not", length(gtool)))
	if (!file.exists(gtool)) stop(paste("Did not find gtool at path:", gtool))

	if (class(request) != "data.frame") stop(paste("request must be data.frame, not", class(request)))
	if (!"chr_name" %in% colnames(request)) stop("request object must have a column 'chr_name'")
	if ("already_exists" %in% colnames(request)) print("request object had a column 'already_exists', this will be overwritten")
	if (!any(substr(rownames(request), 1, 2) %in% "rs")) {
		if (!any(substr(rownames(request), 1, 1) %in% "i")) {
			stop("Not a single rs id was found among the rownames of the request. Really?")
		}
	}

	# check existence of already cached genotypes
	if (file.exists(cachedGenotypeFile)) {
		cachedGenotypes <- read.table(cachedGenotypeFile, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
		snpsAlreadyCached <- rownames(cachedGenotypes)
		requestDeNovo <- request[!rownames(request)%in%snpsAlreadyCached, , drop = FALSE]
	} else {
		requestDeNovo <- request
	}

	# if there are anything novel, extract it from zip (takes a long time)
	if (nrow(requestDeNovo) > 0) {
		# check if raw data is there
		if (!file.exists(genZipFile) | !file.exists(inputZipFile)) {
			stop(safeError(paste("Did not find genome-wide data for this uniqueID. This probably means that raw data has been deleted and no furter SNP-data can be retrieved. Since our policy is to delete personally traceable data, such as genome-wide SNP data, after 14 days, this can often affect users that submitted their data before an update. There is no solution to this, other than re-upload of data.")))
		}

		dir.create(idTempFolder)
		chromosomes <- unique(requestDeNovo[, "chr_name"])
		contents <- unzip(genZipFile, list = TRUE)

		# create a blank genotypes object
		genotypes <- data.frame(genotype = vector(), stringsAsFactors = FALSE)

		# if input is in as a chromosome, use the 23andmefile as input
		if ("input" %in% chromosomes) {
			snpsFromInput <- requestDeNovo[requestDeNovo[, "chr_name"] %in% "input", "SNP"]
			outZip <- unzip(inputZipFile, overwrite = TRUE, exdir = idTempFolder, unzip = "internal")
			cmd0 <- paste("grep -E '", paste(paste(snpsFromInput, "\t", sep = ""), collapse = "|"), "' ", outZip, sep = "")
			input_genotypes <- system(cmd0, intern = TRUE)
			if (length(input_genotypes) > 0) {
				input_genotypes <- do.call(rbind, strsplit(input_genotypes, "\t"))
				input_genotypes[, 4] <- sub("\r$", "", input_genotypes[, 4])
				if (any(nchar(input_genotypes[, 4]) != 2)) {
					print("WARNING: input data should have length 2 genotypes. We try to clean the \r ending once more")
					input_genotypes[, 4] <- sub("\r$", "", input_genotypes[, 4])
				}
				if (any(nchar(input_genotypes[, 4])!=2)) stop("Input data must have length 2 genotypes")

				input_genotypes[, 4] <- paste(substr(input_genotypes[, 4], 1, 1), substr(input_genotypes[, 4], 2, 2), sep = "/")
				genotypes <- data.frame(rsids = input_genotypes[, 1], genotype = input_genotypes[, 4], stringsAsFactors = FALSE)
				genotypes <- genotypes[!duplicated(genotypes[, "rsids"]), ]
				rownames(genotypes) <- genotypes[, "rsids"]
				genotypes[, "rsids"] <- NULL
			}
		}

		# if any normal style chromosome names are in use the gen files
		if (any(c(as.character(1:22), "X") %in% chromosomes)) {
			chromosomes <- chromosomes[chromosomes %in% c(as.character(1:22), "X")]
			chromosomes <- chromosomes[order(suppressWarnings(as.numeric(chromosomes)))]

			gensToExtract <- paste(uniqueID, "_chr", chromosomes, ".gen", sep = "")
			if (!all(gensToExtract %in% contents[, "Name"])) {
				missing <- gensToExtract[!gensToExtract %in% contents[, "Name"]]
				gensToExtract <- gensToExtract[!gensToExtract %in% missing]
				chromosomes <- chromosomes[!chromosomes %in% sub("\\.gen$", "", sub("^.+_chr", "", missing))]
				warning(paste("These were missing in the zip-gen file:", paste(missing, collapse = ", ")))
			}
			outZip <- unzip(genZipFile, files = gensToExtract, overwrite = TRUE, exdir = idTempFolder, unzip = "internal")

			f <- file(paste(idTempFolder, "/samples.txt", sep = ""), "w")
			writeLines("ID_1 ID_2 missing sex", f)
			writeLines("0 0 0 D", f)
			writeLines(paste(uniqueID, uniqueID, "0.0 2 "), f)#gender probably doesn't matter here
			close(f)

			#loop over all chromosomes and extract the relevant genotypes in each using gtool
			for (chr in chromosomes) {
				genotypes_here <- data.frame(row.names = vector(), genotype = vector(), stringsAsFactors = FALSE)

				# this is wrapped in a try block, because it has previously failed from unpredictable memory issues, so it's better to give a few tries
				for (tryCount in 1:3) {
					print(paste("Getting ped and map file at chr", chr, " - try", tryCount))
					gen <- paste(idTempFolder, "/", uniqueID, "_chr", chr, ".gen", sep = "")
					snpsHere <- rownames(requestDeNovo)[requestDeNovo[, "chr_name"] %in% chr]
					write.table(snpsHere, file=paste(idTempFolder, "/snps_in_chr", chr, ".txt", sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE)
					cmd1 <- paste(
						gtool,
						" -S --g " , gen,
						" --s ", idTempFolder, "/samples.txt",
						" --inclusion ", idTempFolder, "/snps_in_chr", chr, ".txt",
						sep=""
					)
					system(cmd1)
					subsetFile <- paste(idTempFolder, "/", uniqueID, "_chr", chr, ".gen.subset", sep = "")
					if (!file.exists(subsetFile)) {
						print(paste("Did not find any of the SNPs on chr", chr))
						next
					}
					cmd2 <- paste(
						gtool, " -G --g " , subsetFile, " --s ", idTempFolder, "/samples.txt --snp --threshold ", call_threshold, sep = ""
					)
					system(cmd2)

					ped <- try(strsplit(readLines(paste(idTempFolder, "/", uniqueID, "_chr", chr, ".gen.subset.ped", sep = "")), "\t")[[1]], silent = TRUE)
					map <- try(read.table(paste(idTempFolder, "/", uniqueID, "_chr", chr, ".gen.subset.map", sep = ""), stringsAsFactors = FALSE), silent = TRUE)

					if (class(ped) != "try-error" & class(map) != "try-error") {
						ped <- ped[7:length(ped)]
						genotypes_here <- try(data.frame(row.names=map[, 2], genotype=sub(" ", "/", ped), stringsAsFactors = FALSE))
						# error where there's duplicate row names
						if (class(genotypes_here) == "try-error") {
							print("Found a duplicate row names error")
							include_these <- which(!duplicated(map[, 2]))
							genotypes_here <- try(data.frame(row.names = map[include_these, 2], genotype = sub(" ", "/", ped[include_these]), stringsAsFactors = FALSE))
						}
						break
					} else {
						genotypes_here <- data.frame(row.names = vector(), genotype = vector(), stringsAsFactors = FALSE)
					}
				}
				genotypes <- rbind(genotypes, genotypes_here)
			}
		}

		if ("N/N" %in% genotypes[, "genotype"]) {
			genotypes[genotypes[, "genotype"] %in% "N/N", "genotype"] <- NA
		}

		stillMissing <- rownames(requestDeNovo)[!rownames(requestDeNovo) %in% rownames(genotypes)]
		genotypes <- rbind(genotypes, data.frame(row.names = stillMissing, genotype = rep(NA, length(stillMissing), stringsAsFactors = FALSE)))

		# remove temporary folder
		unlink(idTempFolder, recursive = TRUE)
	}

	# merge with cachedGenotypes
	if (nrow(requestDeNovo) > 0) {
		if (file.exists(cachedGenotypeFile)) {
			genotypes <- rbind(cachedGenotypes, genotypes)
			unlink(cachedGenotypeFile)
		}
		f <- gzfile(cachedGenotypeFile, "w")
		write.table(genotypes, file = f, sep = "\t", col.names = NA)
		close(f)
	} else {
		genotypes <- cachedGenotypes
	}

	return(genotypes[rownames(request), , drop = FALSE])
}

#' Get SNPs for analyzing
crawl_for_snps_to_analyze <- function (uniqueID = NULL, impute_me = NULL, destinationDir = NULL) {
	if (class(uniqueID) != "character") stop("uniqueID must be of class character")

	# print start message
	cat(paste0(Sys.time(), "\nGetting SNPs for analyzing for: ", uniqueID, "\n"))

	# get list of SNPs to analyze
	all_SNPs <- data.frame(SNP = vector(), chr_name = vector(), stringsAsFactors = FALSE)
	for (module in list.files(impute_me, full.names = TRUE)) {
		if (!file.info(module)["isdir"]) next
		if ("SNPs_to_analyze.txt" %in% list.files(module)) {
			SNPs_to_analyze <- read.table(
				paste(module, "/SNPs_to_analyze.txt", sep = ""), sep = "\t", stringsAsFactors = FALSE, header = TRUE, quote = "", comment = ""
			)
			if (!all(c("SNP", "chr_name") %in% colnames(SNPs_to_analyze))) stop(paste("In", module, "a SNPs_to_analyze file was found that lacked the SNP and chr_name column"))
			SNPs_to_analyze[, "chr_name"] <- as.character(SNPs_to_analyze[, "chr_name"])
			if (!all(SNPs_to_analyze[, "chr_name"] %in% c(1:22, "X", "input"))) stop(paste("In", module, "a SNPs_to_analyze had a chr_name column that contained something else than 1:22 and X"))
			all_SNPs <- rbind(all_SNPs, SNPs_to_analyze[, c("SNP", "chr_name")])
		}
	}

	# extra check for non-discrepant chr info
	if (any(duplicated(all_SNPs[, "SNP"]))) {
		duplicates <- all_SNPs[duplicated(all_SNPs[, "SNP"]), "SNP"]
		for (duplicate in duplicates) {
			if (length(unique(all_SNPs[all_SNPs[, "SNP"] %in% duplicate, "chr_name"])) != 1) stop(paste("Found a multi-entry SNP", duplicate, "with discrepant chr info"))
		}
		all_SNPs <- all_SNPs[!duplicated(all_SNPs[, "SNP"]), ]
	}
	rownames(all_SNPs) <- all_SNPs[, "SNP"]

	print(paste("Checking all requested SNPs from", uniqueID))
	genotypes <- try(get_genotypes(uniqueID = uniqueID, request = all_SNPs, gtool = gtool, destinationDir = destinationDir))
	if (class(genotypes) == "try-error") {
		if (file.exists(paste("/share/project/vinh/genomes/data/", uniqueID, "/temp", sep = ""))) {
			next
		} else {
			print("Some other error happened in the extraction crawler, but probably no cause for alarm.")
			print(genotypes)
		}
	}
	genotypes <- try(get_genotypes(uniqueID = uniqueID, request = all_SNPs, gtool = gtool, destinationDir = destinationDir))

	# get the nonsenser SNPs if possible
	e <- try(load(paste0(impute_me, "/nonsenser/2015-12-16_all_coding_SNPs.rdata")))
	if (class(e) != "try-error") {
		genotypes <- try(get_genotypes(uniqueID, coding_snps, gtool = gtool, destinationDir = destinationDir, namingLabel = "cached.nonsenser"))
	}

	# get the AllDiseases + ukbiobank SNPs if possible
	load(paste0(impute_me, "/AllDiseases/2019-03-04_all_gwas_snps.rdata"))
	e1 <- gwas_snps
	load(paste0(impute_me, "/ukbiobank/2017-09-28_all_ukbiobank_snps.rdata"))
	e2 <- gwas_snps
	e2 <- e2[!rownames(e2) %in% rownames(e1), ]
	e <- rbind(e1, e2)
	if (class(e) != "try-error") {
		genotypes <- try(get_genotypes(uniqueID, e, gtool = gtool, destinationDir = destinationDir, namingLabel = "cached.all_gwas"))
	}

	# get the ethnicity SNPs if possible
	e <- try(load(paste0(impute_me, "/ethnicity/2017-04-03_ethnicity_snps.rdata")))
	if (class(e) != "try-error") {
		genotypes <- try(get_genotypes(uniqueID, ethnicity_snps, gtool = gtool, destinationDir = destinationDir, namingLabel = "cached.ethnicity"))
	}
}

#' run_export scripts
run_export_script <- function (
    uniqueID = NULL, modules = NULL, impute_me = NULL, destinationDir = NULL, functionFile = NULL, gtool = NULL
) {
	require(jsonlite) # for toJSON function

	# print start message
	cat(paste0(Sys.time(), "\nExport test results for: ", uniqueID, "\n"))

	if (class(uniqueID) != "character") stop("UniqueID must be of class character")
	print(destinationDir)
	if (!all(file.exists(destinationDir))) stop("Given output folder was not found")

	if (is.null(modules)) {
		modules <- list.files(impute_me)
	} else {
		if (class(modules) != "character") stop("modules must be of class character")
		if (!all(file.exists(paste(impute_me, "/", modules, sep = "")))) stop ("Not all modules given were found")
	}

	print(paste("Running export script for", uniqueID))
	outputList <- list()
	# import standard pData stuff
	pDataFile <- paste(destinationDir, "/pData.txt", sep = "")
	pData <- try(read.table(pDataFile, header = TRUE, stringsAsFactors = FALSE, sep = "\t"), silent = TRUE)
	if (class(pData) == "try-error") {
		stop(paste("uniqueID", uniqueID, "was skipped due to inavailability of pData file!"))
	}
	if (nrow(pData)!=1) stop ("pData file must have 1 row")

	# check existence of cached file
	cachedFile <- paste(destinationDir, "/", uniqueID, ".cached.gz", sep = "")
	cachedData <- try(read.table(cachedFile, header = TRUE, stringsAsFactors = FALSE), silent = TRUE)
	if (class(cachedData) == "try-error") {
		print(paste("uniqueID", uniqueID, "was skipped due to inavailability of cachedData file"))
		next
	}

	# get basic stuff
	for (imp in c("uniqueID", "first_timeStamp")) {
		if (!imp %in%colnames(pData)) stop(paste("pData lacked this column:", imp))
		outputList[[imp]] <- pData[1, imp]
	}
	outputList[["current_timeStamp"]] <- as.character(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))

	# check if ethnicity is in pData, and if not save it there (because it is needed elsewhere)
	if (!"ethnicity" %in% colnames(pData)) {
		source(paste(paste0(impute_me, "/ethnicity" , "/export_script.R")))
		ethnicity <- try(export_function(uniqueID))
		if (class(ethnicity) == "try-error") {
			ethnicity <- NA
		} else {
			ethnicity <- ethnicity[["guessed_super_pop"]]
		}
		pDataFile <- paste0(destinationDir, "/pData.txt")
		pData <- try(read.table(pDataFile, header = TRUE, stringsAsFactors = FALSE, sep = "\t"), silent = TRUE)
		if (class(pData) != "try-error") {
			pData[1, "ethnicity"] <- ethnicity
			write.table(pData, file = pDataFile, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
			print(paste("Determined and saved ethnicity as:", ethnicity))
		} else {
			print(paste("Couldn't save ethnicity in pData:", ethnicity))
		}
	}

	# get remaining non-ethnicity modules
	module_count <- 0
	for (module in modules) {
		if (!file.info(paste0(impute_me, "/", module))["isdir"]) next
		if ("export_script.R" %in% list.files(paste0(impute_me, "/", module))) {
			print(paste(Sys.time(), "Running module", module, "for", uniqueID))
			if (exists("export_function")) suppressWarnings(rm("export_function"))
			source(paste(paste0(impute_me, "/", module, "/export_script.R")))
			if (!exists("export_function")) stop(paste("In module", module, "there was an export_script.R without an export_function"))
			# exp <- try(export_function(uniqueID))
			exp <- try(export_function(uniqueID, impute_me, destinationDir, functionFile, gtool))
			
			if (class(exp) == "try-error") next
			outputList[[module]] <- exp
			module_count <- module_count + 1
		}
	}

	filename <- paste0(destinationDir, "/", paste(uniqueID, "data.json", sep = "_"))

	# check if there exists previous json file, with module data that is not re-run
	# if so, include this
	if (file.exists(filename)) {
		outputList_previous <- fromJSON(filename)
		previous_unique <- outputList_previous[!names(outputList_previous) %in% names(outputList)]
		if (length(previous_unique) > 0) {
			print(paste("Inserting", length(previous_unique), "modules from existing json:", paste(names(previous_unique), collapse = ", ")))
			outputList <- c(outputList, previous_unique)
		}
	}

	# save new JSON
	JSON <- toJSON(outputList)
	f <- file(filename, "w")
	writeLines(JSON, f)
	close(f)
	return(outputList)
}

#' other functions
genes_for_good_cleaner <- function (rawdata_file) {
    print("The genes_for_good_cleaner was activated")
    if(!file.exists(rawdata_file)) stop(paste("Error in special-error-check: didn't find file at", rawdata_file))
    # Common problem 1 -  # signs in the rsids. Should remove those lines.
    cmd_special_8<-paste0("sed -i.bak6 '/#/d' ", rawdata_file)
    system(cmd_special_8)
}

get_GRS <- function(genotypes, betas){
    # this function is deprecated --- use get_GRS_2 instead
    if (class(genotypes) != "data.frame") stop(paste("genotypes must be data.frame, not", class(genotypes)))
    if (!"genotype" %in% colnames(genotypes)) stop(paste("genotypes must have a column genotype"))
    if (!all(unique(sub("[0-9].+$" ,"" ,rownames(genotypes))) %in% c("i", "rs"))) {
        stop(paste("genotypes must have rownames starting with rs. You had these:", paste(unique(sub("[0-9].+$", "", rownames(genotypes))), collapse = ", ")))
    }
    
    if (class(betas) != "data.frame") stop(paste("genotypes must be data.frame, not", class(betas)))
    necessary_columns <- c("effect_allele", "non_effect_allele", "Beta")
    if (!all(necessary_columns %in% colnames(betas))) stop(paste("betas must have a column", paste(necessary_columns, collapse = ", ")))
    if (!all(unique(sub("[0-9].+$" ,"" ,rownames(betas))) %in% c("i", "rs")))stop("betas must have rownames starting with rs")
    if (!all(rownames(betas) %in% rownames(genotypes))) stop("all SNPs in betas must be present in genotypes")
    
    geneticRiskScore <- 0
    for (snp in rownames(betas)) {
        genotype <- strsplit(genotypes[snp,], "/")[[1]]
        effect_allele <- betas[snp,"effect_allele"]
        non_effect_allele <- betas[snp,"non_effect_allele"]
        beta <- betas[snp, "Beta"]	
        geneticRiskScore <- geneticRiskScore + sum(genotype %in% effect_allele) * beta
    }
    return(geneticRiskScore)
}

get_GRS_2 <- function(snp_data, mean_scale = TRUE, unit_variance = TRUE, verbose = FALSE) {
    # snp_data  a data frame with genotype, effect sizes and information on effect/non-effect allele. Optionally also information about minor allele frequency and minor/major allele (for use with mean scaling etc)
    # mean_scale  logical. If TRUE the GRS output is scaled so that the average person, by MAF-information, will have a score of 0
    # unit_variance logical. If TRUE the GRS output is scaled so that 68% of everyone, by MAF/HWE-information, are within 1 unit of 0 (=1 SD)
    
    if (class(snp_data) != "data.frame") stop(paste("snp_data must be data.frame, not", class(snp_data)))
    if ("Beta" %in% colnames(snp_data) & !"effect_size" %in% colnames(snp_data)) {
        warning("No 'effect_size' column was found, as is necessary per 2017-03-14 - but a 'Beta' column was renamed to 'effect_size'. Do fix in the future")
        colnames(snp_data)[colnames(snp_data) %in% "Beta"] <- "effect_size"
    }
    necessary_columns <- c("genotype", "effect_allele", "non_effect_allele", "effect_size")
    if (!all(necessary_columns %in% colnames(snp_data))) {
        stop(paste("snp_data was lacking necessary columns", paste(necessary_columns[!necessary_columns %in% colnames(snp_data)], collapse=", ")))
    }
    if (!all(unique(sub("[0-9].+$", "", rownames(snp_data))) %in% c("i", "rs"))) stop("snp_data must have rownames starting with rs")
    if (class(snp_data[, "effect_size"])!="numeric") stop("Class of the effect_size column in the snp_data object must be numeric")
    
    
    if (class(mean_scale)!="logical") stop(paste("mean_scale must be logical, not", class(mean_scale)))
    if (length(mean_scale)!=1) stop(paste("mean_scale must be length 1"))
    if (mean_scale) {
        necessary_columns_2 <- c("minor_allele", "major_allele", "minor_allele_freq")
        if (!all(necessary_columns_2 %in% colnames(snp_data))) stop(paste("in mean-scaling, snp_data must have columns", paste(necessary_columns_2, collapse=", ")))
    }
    
    if (class(unit_variance)!="logical") stop(paste("unit_variance must be logical, not", class(unit_variance)))
    if (length(unit_variance)!=1) stop(paste("unit_variance must be length 1"))
    if (!mean_scale & unit_variance) stop("Cannot use unit_variance if not also using mean_scale")
    
    if (class(verbose)!="logical") stop(paste("verbose must be logical, not", class(verbose)))
    if (length(verbose)!=1) stop(paste("verbose must be length 1"))
    
    snp_data[, "personal_score"] <- NA
    snp_data[, "population_score_average"] <- NA
    snp_data[, "population_score_sd"] <- NA
    snp_data[, "score_diff"] <- NA
    missing_snps <- vector()
    missing_major_minor_snps <- vector()
    missing_effect_info_snps <- vector()
    
    for (snp in rownames(snp_data)) {
        # check for missing genotype
        if (is.na(snp_data[snp, "genotype"])) {
            missing_snps <- c(snp, missing_snps) 
            next
        }
        
        # get effect/non-effect-alleles and genotypes
        genotype <- strsplit(snp_data[snp, "genotype"], "/")[[1]]
        effect_allele <- snp_data[snp, "effect_allele"]
        non_effect_allele <- snp_data[snp, "non_effect_allele"]
        
        # check if the effect allele info is missing
        if (any(is.na(c(effect_allele, non_effect_allele))) | any(c(effect_allele, non_effect_allele) %in% "?")) {
            missing_effect_info_snps <- c(missing_effect_info_snps, snp)
            next
        }
        
        # check if the genotype is part of these
        if (!all(genotype %in% c(effect_allele, non_effect_allele))) {
            missing_effect_info_snps <- c(missing_effect_info_snps, snp)
            next
        }
        
        # get effect_size  
        effect_size <- snp_data[snp, "effect_size"]	
        if (is.na(effect_size)) {
            missing_effect_info_snps <- c(missing_effect_info_snps, snp)
            next
        }
        personal_effect_allele_count <- sum(genotype %in% effect_allele)
        snp_data[snp, "personal_score"] <- personal_effect_allele_count * effect_size
        
        # if mean scale, then also calculate the average score in this population (based on MAF) 
        if (mean_scale) {
            # get major/minor/maf info
            major_allele <- snp_data[snp, "major_allele"]
            minor_allele <- snp_data[snp, "minor_allele"]
            minor_allele_freq <- snp_data[snp, "minor_allele_freq"]
            
            # check if they are missing
            if (is.na(major_allele) | is.na(minor_allele) | is.na(minor_allele_freq) | minor_allele=="?" | major_allele=="?") {
                missing_major_minor_snps <- c(snp, missing_major_minor_snps) 
                next
            }
            
            # check if major-minor and effect-non-effect are consistent, and get effect_allele_freq
            if (minor_allele == effect_allele & major_allele == non_effect_allele) {
                effect_allele_freq <- minor_allele_freq
            } else if (minor_allele == non_effect_allele & major_allele == effect_allele) {
                effect_allele_freq <- 1 - minor_allele_freq
            } else {
                stop(paste("discrepancy between effect/non-effect allele and major/minor allele for SNP", snp))
            }
            
            # Calculate what the average score is for this population and also the difference with personal score
            average_effect_allele_count <- effect_allele_freq * 2
            snp_data[snp, "population_score_average"] <- average_effect_allele_count * effect_size
            snp_data[snp, "score_diff"] <- snp_data[snp, "personal_score"] - snp_data[snp, "population_score_average"]
            
            # calculate the extent of possible variance of the score
            # in other words -- Z-scores. The population mean will always be zero... but now we can ensure that 68% (1 SD) is within "1" and 95% (2 SD) is within "2"...
            if (unit_variance) {
                frac_0 <- (1-effect_allele_freq)^2
                frac_1 <- (1-effect_allele_freq)*(effect_allele_freq)*2
                frac_2 <- (effect_allele_freq)^2
                mean <- (frac_1 * 1 * effect_size + frac_2 * 2 * effect_size)
                sigma <- (0*effect_size - mean)^2 * frac_0 + (1*effect_size - mean)^2 * frac_1 + (2*effect_size - mean)^2 * frac_2 
                population_sd <- ( sigma)^0.5
                snp_data[snp, "population_score_sd"] <- population_sd
            }
        }  
    }  
    
    # round values
    for (col in c("personal_score", "population_score_average", "population_score_sd", "score_diff")) {
        snp_data[, col] <- signif (snp_data[, col], 2)
    }
    
    # follow up on the warning message 
    if (length(missing_snps) > 0 & verbose) {
        warning(paste("Note, for", length(missing_snps), "SNPs, we found missing genotypes. This can cause errors particularly if the data is not mean centered. These were skipped:", paste(missing_snps, collapse=", ")))
    }
    if (length(missing_major_minor_snps) > 0 & verbose) {
        warning(paste("Note, for", length(missing_major_minor_snps), "SNPs, we found missing major/minor/freq-allele information. These SNPs were skipped:", paste(missing_major_minor_snps, collapse=", ")))  
    }
    if (length(missing_effect_info_snps) > 0 & verbose) {
        warning(paste("Note, for", length(missing_effect_info_snps), "SNPs, we found wrong or missing information on what was effect-allele and what was non-effect-allele. They were skipped:", paste(missing_effect_info_snps, collapse=", ")))
        
    }
    return(snp_data)
}