################################################################################
##           Floral traits across American biogeographical regions            ##
##                  Plant traits and hummingbird diversity                    ##
################################################################################

# 1. Load Required Libraries ####

# Libraries
libraries <- c(
  "dplyr", "phytools", "RColorBrewer", "tidyverse", "MCMCglmm", "emmeans",
  "greekLetters", "car", "multcompView", "ape", "multcomp", "FD", "reshape2",
  "caper", "cowplot", "FSA"
)

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

# 3. Get hummingbird data, plant data and network Information ####

# Network information
networkInfo <- read.csv("2_raw/networkMetadata.csv", header = TRUE, sep = ",", 
                        row.names=1)

# Plant traits
plants <- read.csv("2_raw/plantTraitsIntroduced.csv", header=TRUE, sep=",",
                   row.names=1)

# hummingbird data
hummer <- read.csv2("2_raw/hummingbirdTraits.csv", header=TRUE, sep=";", dec=".",
                    row.name=1)

# Vector with the biogeographical regions in the order that I want them in
bioGeo <- c("Lowland South America", "Andes", "Central & North America",
            "Caribbean")

#####

# 4. Calculate hummingbird functional trait space ####

# First create a presence matrix for the four biogeographic regions based on networks
# Initialize empty vectors to store collected information
lapply(c("regVec", "specVec", "famVec", "morphVec", "morphVec2",
         "morphVec3", "netVec"), 
       function(x) assign(x, vector(), envir = .GlobalEnv))

# We create a loop, where we repeat each step for each biogeographical region
for(h in 1:length(bioGeo)){
  
  # A vector with the identities of the networks for the biogeographical region
  # we are currently looking at
  ma <- networkInfo[networkInfo$region == bioGeo[h], ]$network
  
  # Now we will go through each row of the data frame "corolla"
  for(i in 1:nrow(hummer)){
    
    # first we check for each row, which networks are listed
    nets <- as.numeric(strsplit(as.character(hummer$network[i]),
                                split = ",")[[1]])
    
    # now we check, if these network(s) are in the specific biogeographical region
    nn <- length(match(ma, nets)[complete.cases(match(networkInfo[networkInfo$region == 
                                                                    bioGeo[h], ]$network, 
                                                      nets))])
    
    # if the network(s) are in the biogeographical region, then we collect the 
    # morpholgical data and some other data for that specific row of "corolla"
    # Plant species
    if(nn > 0){specVec <- append(specVec,
                                 hummer$species[i])} 
    
    # Biogeographical region
    if(nn > 0){regVec <- append(regVec, bioGeo[h])} 
    
    # Bill length
    if(nn > 0){morphVec <- append(morphVec, hummer$billLength[i])} 
    
    # Bill curvature
    if(nn > 0){morphVec2 <- append(morphVec2, hummer$billCurvature[i])} 
    
    # Body weight
    if(nn > 0){morphVec3 <- append(morphVec3, hummer$bodyWeight[i])} 
    
    # Plant family
    if(nn > 0){famVec <- append(famVec, hummer$family[i])} 
    
    # The networks
    if(nn > 0){netVec <- append(netVec, hummer$network[i])} 
    
  } # end i
  
  # Now we join all the information we just collected for one region together
  # in a data frame and then we continue with the next region
  hummerRegion <- data.frame(family=famVec,
                             region=regVec,
                             species=specVec,
                             length=morphVec,
                             curvature=morphVec2,
                             weight=morphVec3,
                             network=netVec)
  
} # end h

# clean
rm(famVec, h, i, ma, morphVec, morphVec2, nets, nn, regVec, specVec,
   morphVec3, netVec)

# Transform into matrix
hummerMatrix <- dcast(hummerRegion, species ~ region, length)
row.names(hummerMatrix) <- hummerMatrix[,1]
hummerMatrix <- hummerMatrix[,-1]

# Separate Matrix, based on network observations
# Separate the rows for unique combinations of region and networks
hummerRegionExpanded <- as.data.frame(hummerRegion %>%
                        # Split the "networks" column into separate rows
                        separate_rows(network, sep = ",") %>%
                        # Remove leading/trailing spaces from network IDs
                        mutate(network = as.integer(trimws(network))) %>%
                        # Check if the network is valid for the corresponding region in networkInfo
                        semi_join(networkInfo, by = c("region" = "region", "network" = "network")) %>%
                        # Combine region and network into a unique identifier
                        mutate(regionNetwork = paste(region, network, sep = "_")) %>%
                        # Optionally arrange for clarity
                        arrange(region, network)
                        )

# Now let's add some abundance data
# Loop through each bird and the network it's currently in
hummerRegionExpanded$abundance <- NA

for (i in 1:nrow(hummerRegionExpanded)) {
  
  # Step 1: Get the bird name
  bird <- hummerRegionExpanded$species[i]
  
  # Step 2: Get the current network and the region
  network <- hummerRegionExpanded$network[i] 
  region <- hummerRegionExpanded$region[i]
  
  # Step 3: Load the network data
  
    # Load the interaction network, handle errors if the file is missing
    interactionNet <- tryCatch({
      suppressWarnings(read.csv(paste0("1_networks/", network, "_", 
                                        networkInfo[networkInfo$network == network, ]$country, ".csv"),
                                sep=";", header=TRUE, row.names=1))
    }, error = function(e) {
      return(NULL)
    })
    
    # Skip if network file not found
    if (is.null(interactionNet)) next
    
    # transform data
    colnames(interactionNet) <- gsub("\\.", " ", colnames(interactionNet))  # Adjust bird names format
    
    # Step 4: get the interaction data of the bird in this net to get an idea of abundance
    hummerRegionExpanded$abundance[i] <- sum(interactionNet[, colnames(interactionNet) == bird])
    
}

# clean
rm(i, network, region, bird)
    
# now change to Matrix
hummerMatrixExpanded <- dcast(hummerRegionExpanded, species ~ regionNetwork, 
                              value.var="abundance", fun.aggregate=sum)
row.names(hummerMatrixExpanded) <- hummerMatrixExpanded[,1]
hummerMatrixExpanded <- hummerMatrixExpanded[,-1]

# get rid of duplicates (i.e. species in multiple regions)
hummerRegion <- hummerRegion[!duplicated(hummerRegion$species), ]

# make a trait data frame, with species as row names
hummerTraits <- hummerRegion[, c(4:6)]
rownames(hummerTraits) <- hummerRegion$species

# let's see if we have NA values. If there are any, replace them by the mean
# for the respective genus.
# Step 1: Loop over each column
for (col in names(hummerTraits)) {
  
  # Step 2: Find rows with NA in the current column
  naRows <- which(is.na(hummerTraits[[col]]))
  
  # Step 3: Extract the first word from the rownames of the NA rows
  genusName <- sapply(strsplit(rownames(hummerTraits)[naRows], " "), `[`, 1)
  
  # Step 4: Loop over each unique first word to calculate the mean and replace NA values
  for (genus in unique(genusName)) {
    
    # Find all rows where the first word matches
    matchingRows <- which(sapply(strsplit(rownames(hummerTraits), " "), `[`, 1) == genus)
    
    # Calculate the mean of the non-NA values in the current column for the matching rows
    meanTrait <- mean(hummerTraits[[col]][matchingRows], na.rm=TRUE)
    
    # Replace NA values with the calculated mean
    hummerTraits[[col]][matchingRows] <- ifelse(is.na(hummerTraits[[col]][matchingRows]), 
                                                 meanTrait, 
                                                 hummerTraits[[col]][matchingRows])
  }

}

# clean
rm(col, genus, genusName, matchingRows, meanTrait, naRows)

# Curvature of Ensifera ensifera is negative. I.e. the bill is curved upwards for that species
# inverse it
hummerTraits[rownames(hummerTraits) == "Ensifera ensifera", ]$curvature <- 2.99

# arrange species in trait and occurence data in the same order
hummerMatrix <- hummerMatrix[c(match(rownames(hummerTraits), rownames(hummerMatrix))), ]
hummerMatrix <- t(hummerMatrix)

# For network based Version as well
hummerMatrixExpanded <- hummerMatrixExpanded[c(match(rownames(hummerTraits), rownames(hummerMatrixExpanded))), ]
hummerMatrixExpanded <- t(hummerMatrixExpanded)

# calculate functional trait space
hummerResults <- dbFD(x=log(hummerTraits), a=sqrt(hummerMatrix), w.abun = TRUE,
                      stand.x = TRUE, corr='none', m="max", stand.FRic=TRUE, 
                      calc.CWM=TRUE, print.pco=TRUE)

# calculate functional trait space based on Networks
hummerResultsNetwork <- dbFD(x=log(hummerTraits), a=sqrt(hummerMatrixExpanded), w.abun=TRUE,
                      stand.x=TRUE, corr='none', m="max", stand.FRic=TRUE, 
                      calc.CWM=TRUE, print.pco=TRUE)

# Store Functional Diversity of Hummingbirds per network
FDhummersNetwork <- as.data.frame(hummerResultsNetwork$FDis)
# Convert row names to a column, then use separate
FDhummersNetwork$region <- str_split_fixed(rownames(FDhummersNetwork), "_", 2)[, 1]
colnames(FDhummersNetwork) <- c("FDis", "region")
# Per region
FDhummersRegion <- as.data.frame(hummerResults$FDis)

# Test for differences between hummingbird FDis of the regions
varianceHFD <- leveneTest(FDis ~ region, data=FDhummersNetwork)
normalityHFD <- shapiro.test(FDhummersNetwork$FDis)

# Perform One way anova
anovaHFD <- aov(log(FDis) ~ region, data=FDhummersNetwork)

# post hoc
tukeyHFD <- TukeyHSD(anovaHFD)
# bonferroni adjusted
pairwise.t.test(FDhummersNetwork$FDis, FDhummersNetwork$region, p.adjust.method="bonferroni")
# Generate letters for groupings based on the p-value threshold
lettersTukeyHFD <- data.frame("Letters"=multcompLetters(extract_p(tukeyHFD)$"region")$"Letters")

# Create a data frame with the regions and their corresponding significance letters
lettersTukeyHFD <- data.frame(
  region = row.names(lettersTukeyHFD),
  label = lettersTukeyHFD
)

# Arrange as will appear on X-Axis
# Convert 'region' to a factor with specified levels and arrange
lettersTukeyHFD <- lettersTukeyHFD %>%
  mutate(region = factor(region, levels = bioGeo)) %>%
  arrange(region)

# Text for X axis Labels
bioGeoModified <- c("Lowland\nSouth America", "Andes", "Central & \nNorth America", "Caribbean")

# Summary of stats
results <- as.data.frame(summary(anovaHFD)[[1]])
resultsText <- bquote(F[.(results$Df[1]) * "," ~ .(results$Df[2])] == .(round(results$'F value'[1], 2)) * "," ~ 
                        .(ifelse(round(results$`Pr(>F)`[1], 2) < 0.05, "p < 0.05", "p > 0.05")))

# Specify the path to the folder
folderPath <- "../figures/4_results/Figure 4"

# Check if the folder exists; if not, create it
if (!dir.exists(folderPath)) {
  dir.create(folderPath, recursive=TRUE)
  message("Folder created at: ", folderPath)
} else {
  message("Folder already exists at: ", folderPath)
}

# Add pch
FDhummersNetwork$pch <- ifelse(FDhummersNetwork$region == "Lowland South America", 21,
                               ifelse(FDhummersNetwork$region == "Andes", 22,
                                      ifelse(FDhummersNetwork$region == "Caribbean", 23, 24)))
# save the figure
png("../figures/4_results/Figure 4/Fig. 4_Supplement_FDis_Hummer_Ylab.png", 
    height=21, width=29.7, units="cm", res=900)

ggplot(FDhummersNetwork, aes(x=region, y=log(FDis), fill=region)) +
  
  geom_violin(trim=FALSE, alpha=0.8) + 
  geom_jitter(width=0.25, size=2, pch=FDhummersNetwork$pch, 
              bg="black") +
  geom_boxplot(width=0.25, fill=adjustcolor("ghostwhite", alpha.f=0.8), 
               outlier.shape=NA) + 
  
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.border=element_rect(colour="black", fill=NA, linewidth=1.5)) +
  
  scale_x_discrete(limits=bioGeo, labels=bioGeoModified) +  
  scale_fill_manual(values=c("beige", "bisque4", "#f1a35f", "cornflowerblue"), 
                    limits=c(bioGeo),  
                    name="") +
  
  scale_y_continuous(
    breaks = seq(-1.5, 1, by=1), 
    labels = scales::label_number(accuracy=0.01) # Multiply by 100, no "%"
  ) +
  
  # Add the significance letters
  geom_text(data=lettersTukeyHFD, aes(x=region, y=max(log(FDhummersNetwork$FDis)+0.5),
                                  label=Letters), size=8, vjust=1) +
  
  labs(x="", y="Log functional diversity of hummingbirds (FDis)",
       title=resultsText) +
  theme(legend.position="none", 
        legend.margin=margin(0, 0, 0, 0), 
        axis.title=element_text(size=25), 
        axis.text.x=element_text(size=25), 
        axis.text.y=element_text(size=20, angle=90, hjust=0.5),
        plot.title=element_text(size=20)) 

dev.off()

# combine the Trait axes and information about where the respective species are
hummer <- hummer[c(match(rownames(hummerTraits), hummer$species)), ]
hummerAxes <- cbind(hummerResults$x.axes[,1:3], t(hummerMatrix), hummerTraits[,1:3],
                    hummer[, c(1, 3)])

# Variance explained by the respective axes
axesVariance <- round(hummerResults$x.values/sum(hummerResults$x.values) * 100, 2)

# correlation of traits with axes
(fit1 <- envfit(hummerResults$x.axes, hummerTraits))

# Extract vector coordinates and labels from the envfit object
vectors <- as.data.frame(fit1$vectors$arrows)
vectors$trait <- rownames(vectors)

# Store the basic hummingbird functional trait space per region
# list with plots
plots <- list()

# Loop through the regions to create plots
for (region in bioGeo) {
  
  # Trait space of avian seed dispersers
  plots[[region]] <- ggplot(data=hummerResults$x.axes[, 1:2], aes(x=A1, y=A2)) +
    
  # Polygon, showing the space occupied by the birds
  geom_polygon(data=as.data.frame(hummerResults$x.axes[, 1:2][chull(hummerResults$x.axes[, 1:2]), ]), 
               aes(x=A1, y=A2), fill="ghostwhite", color="grey50", alpha=0.4, 
               linetype="dashed") +
    
    # Regional trait space
    geom_polygon(data=as.data.frame(hummerAxes[hummerAxes[region] > 0, 1:2][c(chull(hummerAxes[hummerAxes[region] > 0, 1:2])), ]), 
                 aes(x=A1, y=A2), 
                 fill=ifelse(region == "Lowland South America", "beige",
                             ifelse(region == "Andes", "bisque4",
                                    ifelse(region == "Caribbean", "cornflowerblue", "#f1a35f"))), 
                 color="black", lty=1, alpha=0.4) +
    
   # Add the species of the current island
   geom_point(data=hummerAxes[hummerAxes[region] > 0, 1:2],  
              aes(x=A1, y=A2), shape=ifelse(region == "Lowland South America", 21,
                                             ifelse(region == "Andes", 22,
                                                    ifelse(region == "Caribbean", 23, 24))),
              size=2, bg=ifelse(region == "Lowland South America", "beige",
                                ifelse(region == "Andes", "bisque4",
                                       ifelse(region == "Caribbean", "cornflowerblue", "#f1a35f"))),
              col="black") +
    
    # Add vectors
    #geom_segment(data=vectors, 
    #             aes(x=0, y=0, xend=A1, yend=A2), 
    #             arrow=arrow(length=unit(0.2, "cm")), 
    #             color="black", linewidth=0.5) + # Add labels for vectors
    #geom_text(data=vectors, 
    #          aes(x=A1, y=A2, label=trait), 
    #          nudge_x=c(0.2, 0.3, 0.2), nudge_y=0.1, 
    #          color="black", size=4) +
    
    # Community centroid (all)
    #geom_point(data=as.data.frame(t(colMeans(hummerAxes[, 1:3]))),
    #           aes(x=A1, y=A2), bg="grey90", size=4, shape=25) +
    
    # Graphic paramenters
    xlab(paste0("PCoA Axis 1 (", axesVariance[1], "%)")) +
    ylab(paste0("PCoA Axis 2 (", axesVariance[2], "%)")) +
    theme_bw() +
    theme(
      legend.title = element_blank(),
      legend.position = "none",
      axis.title = element_text(size = 15), 
      axis.text.x = element_text(size = 0), 
      axis.text.y = element_text(size = 0),
      plot.title = element_text(size = 12.5),  # Adjust as desired for the plot title size
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour="black", fill=NA, linewidth=1.5)
    ) +

    # Title
    ggtitle(bquote(.(paste0(region))))

}

# Specify the path to the folder
# Specify the path to the folder
folderPath <- "../figures/5_functionalSpace"

# Check if the folder exists; if not, create it
if (!dir.exists(folderPath)) {
  dir.create(folderPath, recursive=TRUE)
  message("Folder created at: ", folderPath)
} else {
  message("Folder already exists at: ", folderPath)
}

# Plot and store some figures
png("../figures/5_functionalSpace/HummingbirdFunctionalSpace_NewColors.png",
    height=21, width=29.7, units="cm", res=900)
# Arrange the plots in a grid
plot_grid(plots[[1]], plots[[2]], plots[[3]], plots[[4]], 
          nrow=2, ncol=2)
dev.off()

# Now find the hummingbird that is in the center (the "average" hummingbird)
# Calculate the Euclidean distance of each row from (0,0)
distances <- apply(hummerAxes[, c(1:2)], 1, function(row) sqrt(sum(row^2)))

# Find the index of the minimum distance
closestToZero <- which.min(distances)

# Get the species closest to (0,0)
centerSpecies <- rownames(hummerAxes)[closestToZero]

# Print the result
centerSpecies

# Clean up
rm(closestToZero, distances, fit1, hummerRegion, hummerTraits,
   hummerMatrix, hummer, centerSpecies, region, plots, vectors,
   folderPath, resultsText, varianceHFD, tukeyHFD, normalityHFD, lettersTukeyHFD,
   anovaHFD)

#####

# 5. Data collection originality and trait matching per plant ####

# Plant df per region

# First add the biogeographical regions to the plants data frame
# Initialize an empty data frame to store results
plantsRegion <- data.frame()

# Loop through each row of the original plants data frame
for (i in 1:nrow(plants)) {
  # Get the network values for each row
  nets <- as.numeric(unlist(strsplit(as.character(plants$network[i]), split = ",")))
  
  # Find the unique regions corresponding to these networks
  regions <- unique(networkInfo$region[networkInfo$network %in% nets])
  
  # If there's only one unique region, we add the row once
  if (length(regions) == 1) {
    rowAdd <- plants[i, ]
    rowAdd$region <- regions
    plantsRegion <- rbind(plantsRegion, rowAdd)
  } else {
    # Otherwise, add the row multiple times for each region
    for (region in regions) {
      rowAdd <- plants[i, ]
      rowAdd$region <- region
      plantsRegion <- rbind(plantsRegion, rowAdd)
    }
  }
}

# Clean up
rm(region, regions, nets, i, rowAdd)

# Go through species per region and calculate metrics
# plant functional "niche breadth"
plantsRegion$nicheBreadthRegion <- NA
plantsRegion$partnersRegion <- NA
plantsRegion$originalityRegion <- NA
plantsRegion$meanCurvature <- NA
plantsRegion$meanLength <- NA
plantsRegion$meanWeight <- NA

# Through regions
for(region in bioGeo){
  
  # Subset data
  regionData <- plantsRegion[plantsRegion$region == region, ]
  
  # new data frame, where plant coordinates will be stored
  # i.e. the coordinates of each plant in the hummers space
  plantCentroidRegion <- data.frame(A1 = numeric(length(regionData$species)), 
                              A2 = numeric(length(regionData$species)), 
                              A3 = numeric(length(regionData$species)), 
                              row.names = regionData$species)

  # With the original plant data, we can get the name of the plant species
  for(i in 1:nrow(regionData)){
    
    # First, we check for each row, which networks are listed
    nets <- as.numeric(strsplit(as.character(regionData$network[i]),
                                split = ",")[[1]])
    
    # region nets
    regionNets <- networkInfo[networkInfo$region == region, ]$network
    
    # keep only the nets that are in the region
    nets <- nets[nets %in% regionNets]
    
    # What are the species that are in this plant family?
    plantSpecies <- as.character(regionData$species[i])
    
    # Initialize an empty vector to store unique hummingbirds
    hummer <- c()
    
    # Now get each network loaded
    for(t in 1:length(nets)){
      
      # Use tryCatch to handle file reading errors
      interactionNet <- tryCatch({
        suppressWarnings(read.csv(paste0("1_networks/", nets[t], "_", 
                                         networkInfo[nets[t], ]$country, ".csv"),
                                  sep=";", header=TRUE, row.names=1))
      },
      error = function(e) {
        # If an error occurs, return NULL or handle it
        return(NULL)
      })
      
      # Check if interactionNet is NULL (indicating the file was not found)
      if (is.null(interactionNet)) {
        next # Skip to the next network in this case
      }
      
      # Bird names of networks need to fit format of hummerTraits
      colnames(interactionNet) <- gsub("\\.", " ",  colnames(interactionNet))
      
      # Now match the species names between the network file and the plant morphology file
      if(plantSpecies %in% row.names(interactionNet)){
        
        # Find the hummingbirds the plants interact with
        hummerTemp <- colnames(interactionNet)[which(interactionNet[which(row.names(interactionNet) == plantSpecies), ] > 0)]
        
        # Append the found hummingbirds to the hummer vector
        hummer <- unique(c(hummer, hummerTemp))  # Ensure uniqueness
        
      } # end if plantSpecies in interactionNet
      
    } # end t loop
    
    # Extract birds from PCoA system
    hummerCoordinates <- hummerAxes[rownames(hummerAxes) %in% hummer, c(1:3)]
    
    # Store the centroid for each plant species
    if (length(hummerCoordinates)>1) hummerCoordinates <- colMeans(hummerCoordinates)
    plantCentroidRegion[i, ] <- hummerCoordinates
    
    # After processing all networks, check how many unique hummingbirds are present
    if (length(hummer) <= 1) {
      # If there's only one or no species, mark it as "only one interaction"
      regionData[i, ]$nicheBreadthRegion <- "one partner" # or use a custom value/message
      
      # Count interaction partners
      regionData[i, ]$partnersRegion <- length(hummer)
      
      # Mean bird traits
      regionData[i, ]$meanCurvature <- mean(hummerAxes[rownames(hummerAxes) %in% hummer, ]$curvature)
      regionData[i, ]$meanLength <- mean(hummerAxes[rownames(hummerAxes) %in% hummer, ]$length)
      regionData[i, ]$meanWeight <- mean(hummerAxes[rownames(hummerAxes) %in% hummer, ]$weight)
      
    } else {
      
      # Create distance matrix
      distanceMatrix <- vegdist(hummerCoordinates, "euclidean")
      
      # Build spanning tree 
      spanningTree <- spantree(distanceMatrix)
      
      # Get a value for nicheBreadth
      regionData[i, ]$nicheBreadthRegion <- sum(spanningTree$dist)
      
      # Count interaction partners
      regionData[i, ]$partnersRegion <- length(hummer)
      
      # Mean bird traits
      regionData[i, ]$meanCurvature <- mean(hummerAxes[rownames(hummerAxes) %in% hummer, ]$curvature)
      regionData[i, ]$meanLength <- mean(hummerAxes[rownames(hummerAxes) %in% hummer, ]$length)
      regionData[i, ]$meanWeight <- mean(hummerAxes[rownames(hummerAxes) %in% hummer, ]$weight)
      
    }
    
  } # end i loop
  
  # Now calculate minimum niche breadth for species with only one intercation partner
  regionData$nicheBreadthRegion <- ifelse(regionData$nicheBreadthRegion == "one partner",
                                as.numeric(min(regionData$nicheBreadthRegion))/2,
                                as.numeric(regionData$nicheBreadthRegion))
  
  # Calculate plant originality
  # First, calculate the community centroid
  plantCentroidRegion <- data.frame(rbind(plantCentroidRegion, 
                                          communityCentroid=as.vector(colMeans(hummerAxes[hummerAxes[region] > 0, 1:3]))))
  
  # Then calculate how far each plant is away
  distanceMatrixPlantCentroid <- as.matrix(dist(plantCentroidRegion, diag=TRUE, upper=TRUE))   
  for (i in 1:dim(distanceMatrixPlantCentroid)[1]){distanceMatrixPlantCentroid[i, i] <- NA}   #put NA in each cell where a species meets itself          
  regionData$originalityRegion <- as.data.frame(distanceMatrixPlantCentroid[, dim(distanceMatrixPlantCentroid)[1]])[, 1][c(1:nrow(regionData))]
  
  # Finally add to the original data frame
  for(plant in regionData$species){ 
    plantsRegion[plantsRegion$region == region & plantsRegion$species == plant, ]$nicheBreadthRegion <- regionData[regionData$species == plant, ]$nicheBreadthRegion
    plantsRegion[plantsRegion$region == region & plantsRegion$species == plant, ]$partnersRegion <- regionData[regionData$species == plant, ]$partnersRegion
    plantsRegion[plantsRegion$region == region & plantsRegion$species == plant, ]$originalityRegion <- regionData[regionData$species == plant, ]$originalityRegion

  }

} 

# Clean again
rm(i, t, region, nets, regionNets, plantSpecies, plant, hummerTemp, hummerCoordinates, hummer,
   distanceMatrix, spanningTree, plantCentroidRegion, distanceMatrixPlantCentroid,
   regionData, interactionNet, meanCurvature, meanLength, meanWeight)

# if there are no interaction partners, it means that there is no data. Lose such species
plantsRegion <- plantsRegion[!plantsRegion$partnersRegion == 0, ]

# log transform metrics
plantsRegion$log10niche <- log10(plantsRegion$nicheBreadthRegion)
plantsRegion$log10originality <- log10(plantsRegion$originalityRegion)

# Scale within each region to account for different number of available partners
# Loop over each region
for (region in bioGeo) {
  # Subset rows for the current region
  subsetDf <- plantsRegion %>% filter(region == !!region)
  
  # Scale the two columns (e.g., var1 and var2) and add as new columns
  subsetDf <- subsetDf %>%
    mutate(
      scaleNiche = scale(log10niche)[, 1],
      scaleOriginality = scale(log10originality)[, 1]
    )
  
  # Reassemble by adding the modified subset back to the original dataframe
  # Using row-binding with rows not in the region, plus the modified subset
  plantsRegion <- plantsRegion %>%
    filter(region != !!region) %>%
    bind_rows(subsetDf)
}

# clean
rm(region, subsetDf)

######

# 6. Functional diversification of plants towards hummingbirds ####

# lists to store data
regionResults <- list(); regionNets <- list()

# I will calculate this only for already published networks
networkInfoPubl <- networkInfo[!networkInfo$network %in% c(95:104), ]

# loop with lesser antillean centroid
# Loop through each region to identify networks
for (region in bioGeo) {
  
  # Step 1: Get all unique networks for the region
  regionNets[[region]] <- sort(unique(as.numeric(unlist(strsplit(as.character(networkInfoPubl[networkInfoPubl$region == region, ]$network), split = ",")))))
  
  # Initialize a vector to store originality means for the current region
  originalityRegion <- c() # completely unweighted
  originalityRegion_1 <- c() # only weight when calculating mean (interaction numbers)
  originalityRegion_2 <- c() # only weight when calculating mean (interaction strength)
  originalityRegion_3 <- c() # weight coordinates and mean (interaction numbers)
  originalityRegion_4 <- c() # weight coordinates and mean (interaction strength)
  
  # Initialize a vector to store originality means for the current region (but with network centroids)
  originalityRegionNet <- c() # completely unweighted
  originalityRegion_1Net <- c() # only weight when calculating mean (interaction numbers)
  originalityRegion_2Net <- c() # only weight when calculating mean (interaction strength)
  originalityRegion_3Net <- c() # weight coordinates and mean (interaction numbers)
  originalityRegion_4Net <- c() # weight coordinates and mean (interaction strength)
  
  # not all nets are possible
  successfulNets <- c()  # Track successful `t` indices
  
  # Step 2: Load each network for the region
  for (t in 1:length(regionNets[[region]])) {
    
    # Load the interaction network, handle errors if the file is missing
    interactionNet <- tryCatch({
      suppressWarnings(read.csv(paste0("1_networks/", regionNets[[region]][[t]], "_", 
                                       networkInfo[networkInfo$network == regionNets[[region]][[t]], ]$country, ".csv"),
                                sep=";", header=TRUE, row.names=1))
    }, error = function(e) {
      return(NULL)
    })
    
    # Skip if network file not found
    if (is.null(interactionNet)) next
    
    # transform data
    #interactionNet <- sqrt(interactionNet)
    
    # Step 3: Create a list to store interactions for each plant
    interactionsNetwork <- list(); 
    interactionsNetwork_3 <- list(); 
    interactionsNetwork_4 <- list(); 
    
    interactionNumber <- list();
    
    colnames(interactionNet) <- gsub("\\.", " ", colnames(interactionNet))  # Adjust bird names format
    
    # Step 4: Process each plant in the network
    for (plant in row.names(interactionNet)) {
      
      # Identify hummingbirds interacting with the plant
      hummers <- colnames(interactionNet)[which(interactionNet[which(row.names(interactionNet) == plant), ] > 0)]
      
      # Retrieve coordinates of hummingbirds
      hummerCoordinates <- hummerAxes[rownames(hummerAxes) %in% hummers, 1:3]
      # Interaction strength
      hummerStrength_3 <- interactionNet[which(row.names(interactionNet) == plant), , drop = FALSE]
      hummerStrength_4 <- hummerStrength_3/sum(interactionNet)
      
      # Weight coordinates
      # Initialize a matrix to store the weighted results, starting as a copy of hummerAxes
      weightedHummerAxes_3 <- hummerCoordinates
      weightedHummerAxes_4 <- hummerCoordinates
      
      # Loop through each hummingbird (column) in the strength vector
      for (hummer in hummers) {
        
        # Check if the hummingbird exists as a row in hummerAxes
        if (hummer %in% rownames(hummerCoordinates)) {
          
          # Extract the row as a numeric vector to multiply it by the strength
          weightedHummerAxes_3[hummer, ] <- hummerCoordinates[hummer, ] * hummerStrength_3[plant, hummer]
          weightedHummerAxes_4[hummer, ] <- hummerCoordinates[hummer, ] * hummerStrength_4[plant, hummer]
          
        }
      }
      
      # Calculate centroid if there are multiple hummingbirds
      if (nrow(hummerCoordinates) > 1) hummerCoordinates <- colMeans(hummerCoordinates)
      if (nrow(weightedHummerAxes_3) > 1) weightedHummerAxes_3 <- colMeans(weightedHummerAxes_3)
      if (nrow(weightedHummerAxes_4) > 1) weightedHummerAxes_4 <- colMeans(weightedHummerAxes_4)
      
      # Store the centroid in the interactionsNetwork list
      interactionsNetwork[[plant]] <- hummerCoordinates
      interactionsNetwork_3[[plant]] <- weightedHummerAxes_3
      interactionsNetwork_4[[plant]] <- weightedHummerAxes_4
      
      interactionNumber[[plant]] <- sum(interactionNet[which(row.names(interactionNet) == plant), which(interactionNet[which(row.names(interactionNet) == plant), ] > 0)])
      
    }
    
    # Convert interactionsNetwork list to a data frame
    interactionsNetwork <- as.data.frame(do.call("rbind", lapply(interactionsNetwork, unlist)))
    interactionsNetwork_3 <- as.data.frame(do.call("rbind", lapply(interactionsNetwork_3, unlist)))
    interactionsNetwork_4 <- as.data.frame(do.call("rbind", lapply(interactionsNetwork_4, unlist)))
    
    # Number of interactions per plant species
    interactionNumber <- as.data.frame(do.call("rbind", lapply(interactionNumber, unlist)))
    
    # Now calcualte metrics
    interactionsNetwork_1 <- colSums(interactionsNetwork * interactionNumber$V1) # weighted by numbers
    interactionsNetwork_2 <- colSums(interactionsNetwork * (interactionNumber$V1/sum(interactionNet))) # weighted by strength
    
    interactionsNetwork_3 <- colSums(interactionsNetwork_3 * interactionNumber$V1) # all weighted by numbers
    interactionsNetwork_4 <- colSums(interactionsNetwork_4 * (interactionNumber$V1/sum(interactionNet))) # all weighted by strength
    
    # Normal - i.e. not weighted at all
    interactionsNetwork <- colMeans(interactionsNetwork)
    
    # Step 5: Add community centroid to interactionsNetwork (Region all)
    communityCentroid <- colMeans(hummerAxes[hummerAxes[region] > 0, 1:3])
    # Step 5: Add community centroid to interactionsNetwork (Region only network)
    communityCentroidNet <- colMeans(hummerAxes[rownames(hummerAxes) %in% colnames(interactionNet), 1:3])
    
    # Add centroid to interactions
    interactionsNetwork <- rbind(interactionsNetwork, communityCentroid)
    interactionsNetworkNet <- rbind(interactionsNetwork, communityCentroidNet)
    interactionsNetwork_1 <- rbind(interactionsNetwork_1, communityCentroid)
    interactionsNetwork_1Net <- rbind(interactionsNetwork_1, communityCentroidNet)
    interactionsNetwork_2 <- rbind(interactionsNetwork_2, communityCentroid)
    interactionsNetwork_2Net <- rbind(interactionsNetwork_2, communityCentroidNet)
    interactionsNetwork_3 <- rbind(interactionsNetwork_3, communityCentroid)
    interactionsNetwork_3Net <- rbind(interactionsNetwork_3, communityCentroidNet)
    interactionsNetwork_4 <- rbind(interactionsNetwork_4, communityCentroid)
    interactionsNetwork_4Net <- rbind(interactionsNetwork_4, communityCentroidNet)
    
    # Step 6: Calculate distances to centroid for originality
    distanceMatrixPlantCentroid <- as.matrix(dist(interactionsNetwork, diag=TRUE, upper=TRUE))
    diag(distanceMatrixPlantCentroid) <- NA  # Set diagonal to NA
    distanceMatrixPlantCentroidNet <- as.matrix(dist(interactionsNetworkNet, diag=TRUE, upper=TRUE))
    diag(distanceMatrixPlantCentroidNet) <- NA  # Set diagonal to NA
    
    distanceMatrixPlantCentroid_1 <- as.matrix(dist(interactionsNetwork_1, diag=TRUE, upper=TRUE))
    diag(distanceMatrixPlantCentroid_1) <- NA  # Set diagonal to NA
    distanceMatrixPlantCentroid_1Net <- as.matrix(dist(interactionsNetwork_1Net, diag=TRUE, upper=TRUE))
    diag(distanceMatrixPlantCentroid_1Net) <- NA  # Set diagonal to NA
    
    distanceMatrixPlantCentroid_2 <- as.matrix(dist(interactionsNetwork_2, diag=TRUE, upper=TRUE))
    diag(distanceMatrixPlantCentroid_2) <- NA  # Set diagonal to NA
    distanceMatrixPlantCentroid_2Net <- as.matrix(dist(interactionsNetwork_2Net, diag=TRUE, upper=TRUE))
    diag(distanceMatrixPlantCentroid_2Net) <- NA  # Set diagonal to NA
    
    distanceMatrixPlantCentroid_3 <- as.matrix(dist(interactionsNetwork_3, diag=TRUE, upper=TRUE))
    diag(distanceMatrixPlantCentroid_3) <- NA  # Set diagonal to NA
    distanceMatrixPlantCentroid_3Net <- as.matrix(dist(interactionsNetwork_3Net, diag=TRUE, upper=TRUE))
    diag(distanceMatrixPlantCentroid_3Net) <- NA  # Set diagonal to NA
    
    distanceMatrixPlantCentroid_4 <- as.matrix(dist(interactionsNetwork_4, diag=TRUE, upper=TRUE))
    diag(distanceMatrixPlantCentroid_4) <- NA  # Set diagonal to NA
    
    distanceMatrixPlantCentroid_4Net <- as.matrix(dist(interactionsNetwork_4Net, diag=TRUE, upper=TRUE))
    diag(distanceMatrixPlantCentroid_4Net) <- NA  # Set diagonal to NA
    
    # Calculate Originality = distance for each plant from the community centroid
    originalityNetwork <- as.data.frame(distanceMatrixPlantCentroid[, dim(distanceMatrixPlantCentroid)[1]])[, 1][c(1:nrow(interactionsNetwork))]
    originalityNetwork_1 <- as.data.frame(distanceMatrixPlantCentroid_1[, dim(distanceMatrixPlantCentroid_1)[1]])[, 1][c(1:nrow(interactionsNetwork_1))]
    originalityNetwork_2 <- as.data.frame(distanceMatrixPlantCentroid_2[, dim(distanceMatrixPlantCentroid_2)[1]])[, 1][c(1:nrow(interactionsNetwork_2))]
    originalityNetwork_3 <- as.data.frame(distanceMatrixPlantCentroid_3[, dim(distanceMatrixPlantCentroid_3)[1]])[, 1][c(1:nrow(interactionsNetwork_3))]
    originalityNetwork_4 <- as.data.frame(distanceMatrixPlantCentroid_4[, dim(distanceMatrixPlantCentroid_4)[1]])[, 1][c(1:nrow(interactionsNetwork_4))]
  
    # Calculate Originality = distance for each plant from the community centroid (only at Network level)
    originalityNetworkNet <- as.data.frame(distanceMatrixPlantCentroidNet[, dim(distanceMatrixPlantCentroidNet)[1]])[, 1][c(1:nrow(interactionsNetworkNet))]
    originalityNetwork_1Net <- as.data.frame(distanceMatrixPlantCentroid_1Net[, dim(distanceMatrixPlantCentroid_1Net)[1]])[, 1][c(1:nrow(interactionsNetwork_1Net))]
    originalityNetwork_2Net <- as.data.frame(distanceMatrixPlantCentroid_2Net[, dim(distanceMatrixPlantCentroid_2Net)[1]])[, 1][c(1:nrow(interactionsNetwork_2Net))]
    originalityNetwork_3Net <- as.data.frame(distanceMatrixPlantCentroid_3Net[, dim(distanceMatrixPlantCentroid_3Net)[1]])[, 1][c(1:nrow(interactionsNetwork_3Net))]
    originalityNetwork_4Net <- as.data.frame(distanceMatrixPlantCentroid_4Net[, dim(distanceMatrixPlantCentroid_4Net)[1]])[, 1][c(1:nrow(interactionsNetwork_4Net))]
    
    # Step 7: Calculate the mean originality per plant in the network
    meanNetwork <- mean(originalityNetwork, na.rm=TRUE) / nrow(interactionNet)
    meanNetwork_1 <- mean(originalityNetwork_1, na.rm=TRUE) / nrow(interactionNet)
    meanNetwork_2 <- mean(originalityNetwork_2, na.rm=TRUE) / nrow(interactionNet)
    meanNetwork_3 <- mean(originalityNetwork_3, na.rm=TRUE) / nrow(interactionNet)
    meanNetwork_4 <- mean(originalityNetwork_4, na.rm=TRUE) / nrow(interactionNet)
       
    # Step 7: Calculate the mean originality per plant in the network (with Net)
    meanNetworkNet <- mean(originalityNetworkNet, na.rm=TRUE) / nrow(interactionNet)
    meanNetwork_1Net <- mean(originalityNetwork_1Net, na.rm=TRUE) / nrow(interactionNet)
    meanNetwork_2Net <- mean(originalityNetwork_2Net, na.rm=TRUE) / nrow(interactionNet)
    meanNetwork_3Net <- mean(originalityNetwork_3Net, na.rm=TRUE) / nrow(interactionNet)
    meanNetwork_4Net <- mean(originalityNetwork_4Net, na.rm=TRUE) / nrow(interactionNet)
    
    # Only add to originalityRegion and mark network as successful if calculation was successful
    originalityRegion <- c(originalityRegion, meanNetwork)
    originalityRegion_1 <- c(originalityRegion_1, meanNetwork_1)
    originalityRegion_2 <- c(originalityRegion_2, meanNetwork_2)
    originalityRegion_3 <- c(originalityRegion_3, meanNetwork_3)
    originalityRegion_4 <- c(originalityRegion_4, meanNetwork_4)
    
    # Same with network only
    originalityRegionNet <- c(originalityRegionNet, meanNetworkNet)
    originalityRegion_1Net <- c(originalityRegion_1Net, meanNetwork_1Net)
    originalityRegion_2Net <- c(originalityRegion_2Net, meanNetwork_2Net)
    originalityRegion_3Net <- c(originalityRegion_3Net, meanNetwork_3Net)
    originalityRegion_4Net <- c(originalityRegion_4Net, meanNetwork_4Net)
    
    successfulNets <- c(successfulNets, regionNets[[region]][[t]])
    
  }
  
  # Step 8: Store results for this region in the list with region and network info
  regionResults[[region]] <- data.frame(
    meanOriginality = originalityRegion,
    meanOriginality_1 = originalityRegion_1,
    meanOriginality_2 = originalityRegion_2,
    meanOriginality_3 = originalityRegion_3,
    meanOriginality_4 = originalityRegion_4,
    meanOriginalityNet = originalityRegionNet,
    meanOriginality_1Net = originalityRegion_1Net,
    meanOriginality_2Net = originalityRegion_2Net,
    meanOriginality_3Net = originalityRegion_3Net,
    meanOriginality_4Net = originalityRegion_4Net,
    network = successfulNets,  # Only the successful networks
    region = region
  )
  
}

# Step 9: Combine all regions into a single data frame
finalResults <- do.call("rbind", regionResults)

# log transform
finalResults[1:10] <- lapply(finalResults[1:10], log)
# region = factor
finalResults$region <- as.factor(finalResults$region)

# Add some plotting parameters: pch
finalResults$pch <- ifelse(finalResults$region == "Lowland South America", 21,
                           ifelse(finalResults$region == "Andes", 22,
                                  ifelse(finalResults$region == "Caribbean", 23, 24)))

#####

# 6.1 Compare regions ####

## First with regional centroid
# We want to compare FDis of plant communities across regions
variance <- leveneTest(meanOriginality_4 ~ region, data=finalResults)
normality <- shapiro.test(finalResults$meanOriginality_4)

# Welch Anova, Kruskal Wallis and normal Anova
oW <- oneway.test(meanOriginality_4 ~ region, data=finalResults, 
                  var.equal=FALSE)
kW <- kruskal.test(meanOriginality_4 ~ region, data=finalResults)
aN <- aov(meanOriginality_4 ~ region, data = finalResults) # weighted by interaction strength, region cetnroid

# Perform Dunn's test (as variance is not homogenous among groups)
dunnTest <- dunnTest(meanOriginality_4 ~ region, data=finalResults, method="bonferroni")

# Format post hoc results with pairwise comparisons and p-values
# Create a matrix of p-values for multcompLetters
groupNames <- unique(unlist(strsplit(dunnTest$res$Comparison, " - ")))
pMatrix <- matrix(1, nrow=length(groupNames), ncol=length(groupNames),
                   dimnames=list(groupNames, groupNames))

# Fill in the matrix with adjusted p-values
for (i in seq_len(nrow(dunnTest$res))) {
  pair <- strsplit(dunnTest$res$Comparison[i], " - ")[[1]]
  pMatrix[pair[1], pair[2]] <- dunnTest$res$P.adj[i]
  pMatrix[pair[2], pair[1]] <- dunnTest$res$P.adj[i]  # Mirror the value for the symmetric matrix
}

# Generate letters for groupings based on the p-value threshold
lettersDunn <- multcompLetters(pMatrix, threshold = 0.05)$Letters

## Then with a network centroid
# We want to compare FDis of plant communities across regions
varianceNet <- leveneTest(meanOriginality_4Net ~ region, data=finalResults)
normalityNet <- shapiro.test(finalResults$meanOriginality_4Net)

# One-way anova
aNNet <- aov(meanOriginality_4Net ~ region, data=finalResults) # weighted by interaction strength, network cetnroid

# Tukey post hoc
postHocNet <- TukeyHSD(aov(meanOriginality_4Net ~ region, finalResults))
lettersTukey <- data.frame("Letters"=multcompLetters(extract_p(TukeyHSD(aov(meanOriginality_4Net ~ region, finalResults))$"region"))$"Letters")

# To be more consistent with the other plots, I want Caribbean to be "b"
lettersTukey$Letters <- c("b", "bc", "a", "ac")

#####

# 6.2 Functional diversification regions violinplot ####

# Region centroid
# Create a data frame with the regions and their corresponding significance letters
lettersDunn <- data.frame(
  region = names(lettersDunn),
  label = lettersDunn
)

# Arrange as will appear on X-Axis
# Convert 'region' to a factor with specified levels and arrange
lettersDunn <- lettersDunn %>%
  mutate(region = factor(region, levels = bioGeo)) %>%
  arrange(region)

# Text for X axis Labels
bioGeoModified <- c("Lowland\nSouth America", "Andes", "Central & \nNorth America", "Caribbean")

# Specify the path to the folder
folderPath <- "../figures/4_results/Figure 4"

# Check if the folder exists; if not, create it
if (!dir.exists(folderPath)) {
  dir.create(folderPath)
  message("Folder created at: ", folderPath)
  } else {
  message("Folder already exists at: ", folderPath)
}

# save the figure
png("../figures/4_results/Figure 4/Fig. 4_NewColors_Ylab.png", 
    height=21, width=29.7, units="cm", res=900)

ggplot(finalResults, aes(x=region, y=meanOriginality_4, fill=region)) +
  
  geom_violin(trim=FALSE, alpha=0.8) + 
  geom_jitter(width=0.1, size=2, pch=finalResults$pch, 
              bg="black") +
  geom_boxplot(width=0.25, fill=adjustcolor("ghostwhite", alpha.f = 0.8), 
               outlier.shape=NA) + 
  
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.border=element_rect(colour="black", fill=NA, linewidth=1.5)) +
  
  scale_x_discrete(limits=bioGeo, labels=bioGeoModified) +  
  scale_fill_manual(values=c("beige", "bisque4", "#f1a35f", "cornflowerblue"), 
                    limits=c(bioGeo),  
                    name="") +
  
  scale_y_continuous(
    breaks = seq(0, min(finalResults$meanOriginality_4 - 2), by = -2), 
    labels = scales::label_number(accuracy=0.01) # Multiply by 100, no "%"
  ) +
  
  # Add the significance letters
  geom_text(data=lettersDunn, aes(x=region, y=max(finalResults$meanOriginality_4+1.5),
                              label=label), size=8, vjust=1) +
  
  labs(x="", y="Log mean distance to community centroid") +
  theme(legend.position="none", 
        legend.margin=margin(0, 0, 0, 0), 
        axis.title=element_text(size=25), 
        axis.text.x=element_text(size=25), 
        axis.text.y=element_text(size=20, angle=90, hjust=0.5),
        plot.title=element_text(size=30)) 

dev.off()

#####