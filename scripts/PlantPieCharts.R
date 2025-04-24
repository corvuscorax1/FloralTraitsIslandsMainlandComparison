################################################################################
##           Floral traits across American biogeographical regions            ##
##               Plant pie charts - distribution of families                  ##
################################################################################

# 1. Libraries ####

# Ensure necessary libraries are installed and loaded
libraries <- c("RColorBrewer")

installedLibs <- libraries %in% rownames(installed.packages())
if (any(!installedLibs)) {
  install.packages(libraries[!installedLibs])
}

invisible(lapply(libraries, library, character.only = TRUE))
rm(installedLibs, libraries)

#####

# 2. Set Working Directory ####

setwd("../data")

#####

# 3. Data ####

# Load and clean plant data
plants <- read.csv("2_raw/plantTraitsIntroduced.csv", header=TRUE, sep=",", row.names=1)

# Load network metadata
netInfo <- read.csv("2_raw/networkMetadata.csv", header=TRUE, sep=",")

#####

# 4. Frequency of families ####

# Extract top 11 most frequent plant families
plantClades <- sort(table(plants$family), decreasing=TRUE)[1:11]
allPlants <- c("Others", names(plantClades))

#####

# 5. Pie charts ####

# Define color scheme for plant families
plantColor <- data.frame(clade=allPlants)
plantColor$color <- c("ghostwhite", hcl.colors(11, palette="PuBuGn"))

# Define biogeographical regions
regions <- c("Lowland South America", "Andes", "Central & North America", "Caribbean")

# Initialize data frame to store frequency of families per region
plantFreqReg <- data.frame(clade=NA, freq=NA, region=NA)

# Plotting parameters
quartz(4, 15); par(mfrow=c(4, 1), mar=c(1, 1, 1, 1) + 1.1)

# Generate pie charts per region
for (t in 1:4) {
  # Subset the network IDs for the current region
  reg <- subset(netInfo$network, netInfo$region == regions[t])
  
  # Initialize an empty vector to store plant families for the current region
  clade <- vector()
  
  # Loop through each plant to check if its Network.ID belongs to the current region
  for (i in 1:nrow(plants)) {
    # Split the Network.ID into individual components
    networkIds <- as.numeric(unlist(strsplit(as.character(plants$network[i]), ",")))
    
    # Check if any of the Network.ID components match the region IDs
    if (any(networkIds %in% reg)) {
      clade <- c(clade, as.character(plants$family[i]))
    }
  }
  
  # Replace plant families that are not in the 'allPlants' vector with "Others"
  clade[!clade %in% allPlants[-1]] <- "Others"
  
  # Create a frequency table of the plant families
  freqTable <- table(clade)
  
  # Update frequencies in the plantColor data frame
  plantColor$freq <- sapply(plantColor$clade, function(cladeName) {
    if (cladeName %in% names(freqTable)) {
      freqTable[cladeName]
    } else {
      0
    }
  })
  
  # Add the region to the plantColor data frame
  plantColor$region <- regions[t]
  
  # Create the pie chart for the current region
  pie(plantColor$freq, col=plantColor$color, 
      main=paste0(regions[t], "\nn = ", length(clade)), labels=NA)
  
  # Append the frequency information for the current region to the plantFreqReg data frame
  plantFreqReg <- rbind(plantFreqReg, plantColor[c("clade", "freq", "region")])
}

# Save pie charts as PDF

# Specify the path to the folder
folderPath <- "../figures/2_pieChart"

# Check if the folder exists; if not, create it
if (!dir.exists(folderPath)) {
  dir.create(folderPath)
  message("Folder created at: ", folderPath)
} else {
  message("Folder already exists at: ", folderPath)
}

invisible(dev.copy2pdf(file="../figures/2_pieChart/PieChartPlants.pdf", 
                       useDingbats=FALSE, family="sans"))
invisible(dev.off())

# Create legend for the pie charts
quartz(4, 6)

pdf(file="../figures/2_pieChart/plantLegend.pdf", 
    width=3, height=4.5, useDingbats=FALSE, family="sans")
# Initialize a blank plot area
plot.new()
legend("topright", legend=plantColor$clade, fill=plantColor$color)
invisible(dev.off())

# clean up again
rm(plantColor, allPlants, clade, i, networkIds, reg, t)

#####

# 6. Test for regional differences ####

# 1. Remove NAs from plantFreqReg
plantFreqReg <- plantFreqReg[complete.cases(plantFreqReg), ]

# 2. Create a frequency matrix
freqMatrix <- matrix(
  plantFreqReg$freq,  
  nrow = 12,         
  ncol = 4,          
  dimnames = list(Family=plantFreqReg$clade[1:12], Region=regions)
)

# 3. Run Fisher's exact test for all regions
fisherTestAll <- fisher.test(freqMatrix, simulate.p.value=TRUE, B=10000)

# 4. Define region combinations for pairwise tests
regionComb <- combn(regions, 2, simplify = FALSE)

# 5. Initialize a list to store results
fisheResults <- list(fisherTestAll=fisherTestAll)

# 6. Perform pairwise Fisher's exact tests
summaryFishers <- capture.output(print(fisherTestAll))
for (regionsPair in regionComb) {
  regionIndices <- which(regions %in% regionsPair)
  testResult <- fisher.test(freqMatrix[, regionIndices], simulate.p.value = TRUE, B = 10000)
  
  regionHeader <- paste0("$`Pairwise test between ", regionsPair[1], " and ", regionsPair[2], "`")
  summaryFishers <- c(summaryFishers, regionHeader, capture.output(print(testResult)))
}

# 7. Write results to a text file

# Specify the path to the folder
folderPath <- "../results/textFiles/pieChart"

# Check if the folder exists; if not, create it
if (!dir.exists(folderPath)) {
  dir.create(folderPath, recursive=TRUE)
  message("Folder created at: ", folderPath)
} else {
  message("Folder already exists at: ", folderPath)
}

writeLines(summaryFishers, con="../results/textFiles/pieChart/FisherTestsPlants.txt")

# Clean up
rm(list = ls())

#####
