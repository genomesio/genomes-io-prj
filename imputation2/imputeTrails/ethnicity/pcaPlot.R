library(plotly)
library(jsonlite)

# get required data
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    msg <- paste(
        "You need to give all these args:",
        "1: output json file",
        "2: impute-me/ethnicity/2017-04-03_ethnicity_pca.rdata file",
        "3: impute-me/ethnicity/2017-04-03_ethnicity_descriptions.txt file"
    )
    stop(paste0(msg, "\n"), call. = FALSE)
}

json_file <- args[1]
ethnicity_pca_file <- args[1]
ethnicity_desc_file <- args[1]

# load ethnicity from json file
d1 <- fromJSON(json_file)
d2 <- d1[["ethnicity"]]
you <- data.frame(d2)
you <- you[, 1:4]
colnames(you) <- c("super_pop", "pos_PC1", "pos_PC2", "pos_PC3")
you$pop <- "YOU"
you$gender <- NA
you$pos_PC4 <- 0
you$pos_PC5 <- 0
you <- you[, c(5, 1, 6, 2, 3, 4, 7, 8)]

hint_message <- ""
if ("guessed_super_pop" %in% names(d2)) {
	d3 <- d2[["guessed_super_pop"]]
	col <- c('light-blue', 'green', 'red', 'purple', 'orange')
	names(col) <- c('AFR', 'AMR', 'EAS', 'EUR', 'SAS')
	if(!is.na(d3)){
		if(d3 %in% names(col)){
			hint_message <- paste0("Your genome (with ", d3, " as predicted super population) locates in the ", col[d3], " cluster.")
		}
	}
}

# add user ethnicity to reference PCA data
load(ethnicity_pca_file)
pca <- rbind(pca_data, you)

# pick some colours for each super population (first dilute their alpha a little)
ethnicity_desc <- read.table(ethnicity_desc_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
colours <- ethnicity_desc[, "Col"]
names(colours) <- ethnicity_desc[, "PopulationDescription"]
# also get the description of each populations
pca[, "pop_long"] <- ethnicity_desc[pca[, "pop"], "PopulationDescription"]

# extract relevant data
pca[, "sizes"] <- c(rep(0.5, nrow(pca)-1), 2)
pca[, "x"] <- pca[, "pos_PC1"]
pca[, "y"] <- pca[, "pos_PC2"]
pca[, "z"] <- pca[, "pos_PC3"]

# generate plot
plot_ly(
	pca, x = ~x, y = ~y, z = ~z, type = "scatter3d", mode = "markers",
	color= ~pop_long, colors = colours,
	showlegend = FALSE, size = ~sizes, marker = list(symbol = 'circle', sizemode = 'diameter'),
	sizes = c(4, 10), hoverinfo = 'text',  text = pca[, "pop_long"]
) %>%
	layout(
	    title = '',
	    scene = list(
	        xaxis = list(
	            title = "PC1",
	            gridcolor = 'rgb(255, 255, 255)',
	            gridwidth = 2
	        ),
	        yaxis = list(
	            title = "PC2",
	            gridcolor = 'rgb(255, 255, 255)',
	            gridwith = 2
	        ),
	        zaxis = list(
	            title = "PC3",
	            gridcolor = 'rgb(255, 255, 255)',
	            gridwith = 2
	        )
	    ),
	    paper_bgcolor = 'rgb(243, 243, 243)',
	    plot_bgcolor = 'rgb(243, 243, 243)'
	)

message(paste(
    "Your genome is indicated as a slightly larger black dot in this plot, you may have to zoom in to see it.",
    hint_message,
    "Super population colors: AFRICAN light-blue, AMERICAN green, EAST ASIAN red, EUROPEAN purple, SOUTH ASIAN orange",
    sep = "\n"
))
