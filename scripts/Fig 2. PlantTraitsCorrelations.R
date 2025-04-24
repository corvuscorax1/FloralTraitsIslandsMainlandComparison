################################################################################
##           Floral traits across American biogeographical regions            ##
##                       Plant traits across regions                          ##
################################################################################

# 1. Load Required Libraries ####

# Libraries
libraries <- c(
  "dplyr", "phytools", "RColorBrewer", "tidyverse", "emmeans",
  "car", "multcompView", "ape", "multcomp", "lme4", "MuMIn",
  "glmmTMB", "caper", "betareg", "cowplot"
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

# 3. Get plant data and store network Information ####

# Network information
networkInfo <- read.csv("2_raw/networkMetadata.csv", header = TRUE, sep = ",", 
                        row.names=1)

# Plant traits
plants <- read.csv("2_raw/plantTraitsIntroduced.csv", header=TRUE, sep=",",
                   row.names=1)

# Vector with the biogeographical regions in the order that I want them in
bioGeo <- c("Lowland South America", "Andes", "Central & North America",
            "Caribbean")

#####

# 4. Corolla and Nectar
# 4.1 Data collection corolla length and nectar concentration ####

# Dataframe, with only the species with corolla length AND nectar concentration
corollaNectar <- plants[!is.na(plants$corolla) & !is.na(plants$nectar), ]

# Additionally lose entries with nectar < 5% as they are probably not correct
corollaNectar <- corollaNectar[!corollaNectar$nectar < 5, ]

# A bunch of vectors, where I will store information about traits, region etc.
lapply(c("regVec", "specVec", "morphVec", "morphVec2", "famVec"), 
       function(x) assign(x, vector(), envir = .GlobalEnv))

# We create a loop, where we repeat each step for each biogeographical region
for(h in 1:length(bioGeo)){
  
  # A vector with the identities of the networks for the biogeographical region
  # we are currently looking at
  ma <- networkInfo[networkInfo$region == bioGeo[h], ]$network
  
  # Now go through each row of the data frame "corolla"
  for(i in 1:nrow(corollaNectar)){
    
    # first  check for each row, which networks are listed
    nets <- as.numeric(strsplit(as.character(corollaNectar$network[i]), 
                                split=",")[[1]])
    
    # Check if networks 1 or 2 are present in 'nets'
    #if(1 %in% nets || 2 %in% nets){
    #  next # Skip this iteration if networks 1 or 2 are found
    #}
    
    # now check, if these network(s) are in the specific biogeographical region
    nn <- length(match(ma, nets)
                 [complete.cases(match(networkInfo[networkInfo$region ==
                                                     bioGeo[h], ]$network, nets))])
  
    # if the network(s) are in the biogeographical region, then collect the data
    if(nn > 0){
      specVec <- append(specVec, corollaNectar$species[i]) # plant name
      regVec <- append(regVec, bioGeo[h]) # the region
      morphVec <- append(morphVec, corollaNectar$corolla[i]) # corolla
      morphVec2 <- append(morphVec2, corollaNectar$nectar[i]) # nectar
      famVec <- append(famVec, corollaNectar$family[i]) # plant family
    }
    
  } # end i
  
  # Now we join all the information we just collected for one region together
  # in a data frame and then we continue with the next region
  corollaNectarRegion <- data.frame(family=famVec,
                              region=regVec,
                              corolla=morphVec,
                              nectar=morphVec2,
                              species=specVec)
  
} # end h

# clean up
rm(famVec, h, i, ma, morphVec, morphVec2, nets, nn, regVec, specVec)

# Adapt the plant names to the format of the phylogenetic tree
corollaNectarRegion$species <- gsub(" ", "_", corollaNectarRegion$species)

# log10 transform the corolla length
corollaNectarRegion$log10corolla <- log10(corollaNectarRegion$corolla + 1)
# transform nectar to range between 0 and 1
corollaNectarRegion$nectar10 <- corollaNectarRegion$nectar / 100

# throw out any duplicate entries (some species occur in multiple networks in
# the same biogeographical region)
corollaNectarRegion <- distinct(corollaNectarRegion)

#####

# 4.2 Match data corolla and phylogenetic tree ####

# Read the plant tree
plantTree <- read.tree("4_phylogeny/PlantTree")

# Cut the tree down to the species that we have in our dataset
plantTree <- drop.tip(plantTree, setdiff(plantTree$tip.label,
                                         corollaNectarRegion$species))
# How many species were matched
setdiff(unique(corollaNectarRegion$species), plantTree$tip.label)

# if a species occurs in multiple regions, including Caribbean, keep only
# Caribbean
multipleSpecies <- names(table(corollaNectarRegion$species)[table(corollaNectarRegion$species) > 1])
# Filter multipleSpecies to keep only those species that occur in the 
# "Caribbean" region
caribbeanSpecies <- multipleSpecies[sapply(multipleSpecies, function(species) {
  "Caribbean" %in% corollaNectarRegion[corollaNectarRegion$species == species, ]$region })]
# Now keep only those rows, where the biogeographical region is Caribbean
corollaNectarRegion <- corollaNectarRegion[!(corollaNectarRegion$species 
                                             %in% caribbeanSpecies & 
                                               corollaNectarRegion$region != "Caribbean"), ]

# then we need to sort the trait data so that the species follow the order of the tree
corollaNectarRegion <- corollaNectarRegion[match(plantTree$tip.label,
                                                 corollaNectarRegion$species), ]

# Clean working environment
rm(caribbeanSpecies, multipleSpecies)

#####

# 4.3 Corolla Nectar Model ####

# PGLS - First test for phylogenetic signal in the data
# Prepare data
# Remove node labels by setting them to NULL
plantTree$node.label <- NULL
# Phylogenetic object
phyloData <- comparative.data(phy=plantTree, 
                               data=corollaNectarRegion, 
                               names.col="species") 
# PGLS Model 1
pglsModel <- pgls(nectar10 ~ log10corolla, 
                  data=phyloData, 
                  lambda="ML")

# Store results
# Specify the path to the folder
folderPath <- "../results/Figure 2"

# Check if the folder exists; if not, create it
if (!dir.exists(folderPath)) {
  dir.create(folderPath)
  message("Folder created at: ", folderPath)
} else {
  message("Folder already exists at: ", folderPath)
}

# Extract residuals corrected for phylogeny
residuals <- pglsModel$residuals

# Get the variance-covariance matrix of the tree to control for phylogeny
vcvMatrix <- vcv(plantTree)

# Apply Cholesky transformation to control for phylogeny
residualsPhylo <- chol(solve(vcvMatrix)) %*% residuals

# Save some diagnostic plots
# Create individual plots using ggplot2
# 1. Residuals vs Fitted Plot
plot1 <- ggplot(data=data.frame(Fitted=pglsModel$fitted, Residuals=residualsPhylo), 
                aes(x=Fitted, y=Residuals)) +
  geom_point(size=2) +
  geom_hline(yintercept=0, color="red") +
  ggtitle("Residuals vs Fitted") +
  theme_minimal()

# 2. Q-Q Plot
plot2 <- ggplot(data=data.frame(Residuals=residualsPhylo), aes(sample=Residuals)) +
  stat_qq(size=2, color="black") +
  stat_qq_line(color="indianred3", size=1) +
  ggtitle("Q-Q Plot of Residuals") +
  theme_minimal()

# 3. Histogram of Residuals
plot3 <- ggplot(data=data.frame(Residuals=residualsPhylo), aes(x=Residuals)) +
  geom_histogram(bins=30, fill="gray", color="black") +
  ggtitle("Histogram of Residuals") +
  theme_minimal()

png("../results/Figure 2/diagnosticsNectarCorolla.png", height=14, width=14, 
    units="cm", res=900)
plot_grid(plot1, plot2, plot3, ncol = 2)
dev.off()

# Since data is not normal, I check also with a beta-regression
# beta regression
bR <- betareg(nectar10 ~ log10corolla, corollaNectarRegion)
# Summary
summary(bR)

# Store model output
# Create the summary with structured content
summaryCorollaNectar <- paste(
  # Title
  "Corolla Length and Nectar Concentration Analysis",
  "\n===============================================\n",
  
  # Effect size (Beta coefficient) and p-value from pglsModel
  "Effect Size and Significance",
  "\n----------------------------\n",
  sprintf("Beta (ß) = %.2f", summary(pglsModel)$coefficients[2, "Estimate"]),
  sprintf("p-value = %.3f", summary(pglsModel)$coefficients[2, "Pr(>|t|)"]),
  
  # R-squared
  "\n\nR-squared",
  "\n---------\n",
  sprintf("R-squared = %.2f", summary(pglsModel)$r.squared),
  
  # Full Model Output: pglsModel
  "\n\nModel Output: PGLS",
  "\n--------------------\n",
  paste(capture.output(summary(pglsModel)), collapse="\n"),
  
  # Full Model Output: Beta regression (bR)
  "\n\nModel Output: Beta Regression",
  "\n------------------------------\n",
  paste(capture.output(summary(bR)), collapse="\n"),
  
  # End of summary
  "\n"
)

# Write the summary to a text file

# create folder
folderPath <- "../results/textFiles/Figure 2"

# Check if the folder exists; if not, create it
if (!dir.exists(folderPath)) {
  dir.create(folderPath)
  message("Folder created at: ", folderPath)
} else {
  message("Folder already exists at: ", folderPath)
}

writeLines(summaryCorollaNectar, "../results/textFiles/Figure 2/summaryNectarCorolla.txt")

# Clean up
rm(summaryCorollaNectar, phyloData, vcvMatrix, residualsPhylo, residuals,
   bR, plantTree, plot1, plot2, plot3)

#####

# 4.4 Corolla Nectar scatterplot ####

# First I create some graphics parameters, pch and color
# Add some plotting parameters: pch
corollaNectarRegion$pch <- ifelse(corollaNectarRegion$region == "Lowland South America", 21,
                                  ifelse(corollaNectarRegion$region == "Andes", 22,
                                         ifelse(corollaNectarRegion$region == "Caribbean", 23, 
                                                ifelse(corollaNectarRegion$region == "All", NA, 24))))

# Add some plotting parameters: color
corollaNectarRegion$col <- ifelse(corollaNectarRegion$region == "Lowland South America", "beige",
                                  ifelse(corollaNectarRegion$region == "Andes", "bisque4",
                                         ifelse(corollaNectarRegion$region == "Caribbean", "cornflowerblue", "#f1a35f")))

# Create a new data frame for predictions
newData <- data.frame(log10corolla = seq(min(corollaNectarRegion$log10corolla),
                                         max(corollaNectarRegion$log10corolla), 
                                         length.out=nrow(corollaNectarRegion)))

# Add predictions from the quadratic model
newData$predict <- predict(pglsModel, newData)

# Store figure
figureManuscript <- list()

figureManuscript[[1]] <- ggplot(corollaNectarRegion, aes(x=log10corolla, y=nectar10, size=4)) +
  
  # Graphic parameters of plotting area
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border=element_rect(colour="black", fill=NA, size=1.5)) +
 
  # Add Caribbean on top to make them more visible
  geom_point(data=corollaNectarRegion, bg=corollaNectarRegion$col, col = "black", pch=corollaNectarRegion$pch, 
             alpha=0.7) +
  
  geom_line(data=newData, aes(x=log10corolla, y=predict), 
            color="black", linewidth=1, 
            lty=ifelse(summary(pglsModel)$coefficients[2, ][4] < 0.05, 1, 2)) +
  
  scale_x_continuous(
    breaks = seq(0, max(corollaNectarRegion$log10corolla), by = 0.5), # Adjust as needed
    labels = scales::label_number(accuracy = 0.01) # Ensure 2 decimal places
  ) +
  
  scale_y_continuous(
    breaks = seq(0, max(corollaNectarRegion$nectar10 + 0.2), by = 0.1), # Adjust as needed
    labels = scales::label_number(accuracy = 0.01, scale=100) # Ensure 2 decimal places
  ) +
  
  # Axes
  xlab("Log corolla length (mm)") +
  ylab("Nectar concentration (%)") +
  theme(legend.title = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.title = element_text(size=30), 
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=20, angle=90, hjust=0.5),
        plot.title = element_text(size=0)) 
  
# Specify the path to the folder
folderPath <- "../figures/4_results/Figure 2"

# Check if the folder exists; if not, create it
if (!dir.exists(folderPath)) {
  dir.create(folderPath, recursive=TRUE)
  message("Folder created at: ", folderPath)
} else {
  message("Folder already exists at: ", folderPath)
}

# Plot and store
png("../figures/4_results/Figure 2/Fig. 2a) CorollaNectarScatterplot_NewColors_Ylab.png", height=21, 
    width=29.7, units="cm", res=900)
print(figureManuscript[[1]])
dev.off()

# clean a bit
rm(newData, folderPath, corollaNectar)

#####

# 4.4.1 Corolla Nectar scatterplot - Region ####

# First I create some graphics parameters, pch and color

# Read the plant tree
plantTree <- read.tree("4_phylogeny/PlantTree")

# Create a list to store models and plots
models <- list()
plotsCorolla <- list()
residuals <- list()

# Go through each region separately
for (i in seq_along(rev(bioGeo))) {
  region <- rev(bioGeo)[i]
  
  # Subset data for the current region
  regionData <- corollaNectarRegion %>% filter(region == !!region)

  # Region tree
  # Cut the tree down to the species that we have in our dataset
  regionTree <- drop.tip(plantTree, setdiff(plantTree$tip.label,
                                            regionData$species))
  # How many species were matched
  setdiff(unique(regionData$species), regionTree$tip.label)
  
  # then we need to sort the trait data so that the species follow the order of the tree
  regionData <- regionData[match(regionTree$tip.label,
                                 regionData$species), ]
  
  # Remove node labels by setting them to NULL
  regionTree$node.label <- NULL
  
  # Phylogenetic object
  phyloDataRegion <- comparative.data(phy=regionTree, 
                                      data=regionData, vcv=TRUE,
                                      names.col="species") 
  # PGLS Model 1
  models[[region]] <- pgls(nectar10 ~ log10corolla, 
                           data=phyloDataRegion, 
                           lambda="ML")
  
  # Check residuals
  # Extract residuals corrected for phylogeny
  residualsModel <- models[[region]]$residuals
  
  # Get the variance-covariance matrix of the tree to control for phylogeny
  vcvMatrix <- vcv(regionTree)
  
  # Apply Cholesky transformation to control for phylogeny
  residualsPhylo <- chol(solve(vcvMatrix)) %*% residualsModel
  
  # QQ plot
  residuals[[region]] <- ggplot(data.frame(residuals=residualsPhylo), 
                                aes(sample=residuals)) +
    geom_qq() +
    geom_qq_line(color="indianred") +
    ggtitle(paste("Q-Q Plot ", region)) +
    theme_bw()
  
  # Create a new data frame for predictions
  newData <- data.frame(log10corolla = seq(min(regionData$log10corolla),
                                           max(regionData$log10corolla), 
                                           length.out=nrow(regionData)))
  
  # Add predictions from the quadratic model
  newData$predict <- predict(models[[region]], newData)
  
  # Create a ggplot for the current region
  plotsCorolla[[region]] <- ggplot(regionData, aes(x=log10corolla, y=nectar10, shape=region, 
                                                   colour=region, fill=region, size=4)) +
    
    # No background, thick lines around plot
    theme_bw() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          panel.border = element_rect(colour="black", fill=NA, linewidth=1.5)) +
    
    # set standard plot covering the potential range of all data
    geom_point(data=corollaNectarRegion[!corollaNectarRegion$region == region, ], size=3, pch=20,
               col=rgb(40, 60, 80, alpha=25, maxColorValu=250)) +
    
    # Now actual region
    geom_point(bg=regionData$col, col="black", pch=regionData$pch, 
               alpha = 0.7) +
    
    
    geom_line(data=newData, aes(x=log10corolla, y=predict), 
              color="black", linewidth=1, 
              lty=ifelse(summary(models[[region]])$coefficients[2, ][4] < 0.05, 1, 2)) +
    
    # Axis labels
    xlab(ifelse(i == length(bioGeo), "Corolla length", "")) +
    ylab("Nectar concentration") +
    theme(legend.title = element_blank(), legend.position = "none",
          axis.title = element_text(size=18), 
          axis.text.x = element_text(size=0), 
          axis.text.y = element_text(size=0),
          plot.title = element_text(size=18)) +
    
    # Title
    ggtitle(bquote(.(region) * ": " * beta * "=" * .(round(summary(models[[region]])$coefficients[2, 1], 2)) * "," ~
                     R^2 * "=" * .(round(summary(models[[region]])$r.squared, 2)) * ", p=" * .(round(summary(models[[region]])$coefficients[2, 4], 3))))
  
}

# Specify the path to the folder
folderPath <- "../figures/4_results/Figure 2/Region"

# Check if the folder exists; if not, create it
if (!dir.exists(folderPath)) {
  dir.create(folderPath, recursive=TRUE)
  message("Folder created at: ", folderPath)
} else {
  message("Folder already exists at: ", folderPath)
}

# Display plots (example to show the first plot)
# Plot and store
png("../figures/4_results/Figure 2/Region/SupplFig. 2a) CorollaNectarScatterplotRegion_NewColors.png",
    height=29, width=18, units="cm", res=900)
# Arrange the plots in a grid
plot_grid(plotsCorolla[[1]], plotsCorolla[[2]], 
          plotsCorolla[[3]], plotsCorolla[[4]], 
          nrow=4, ncol=1)
dev.off()

# Clean working environment
rm(i, models, region, vcvMatrix, residualsPhylo, residualsModel,
   residuals, regionTree, regionData, phyloDataRegion, plantTree,
   plotsCorolla, pglsModel, newData)

#####

# 4.4.2 Corolla Nectar scatterplot - Family ####

# First I create some graphics parameters, pch and color

# Create a list to store models and plots
models <- list()
plots <- list()
residuals <- list()

# Unique families
# Count species per family and order by species number
families <- as.data.frame(corollaNectarRegion %>%
                            group_by(family) %>%
                            summarise(species_count = n()) %>%
                            arrange(desc(species_count)))$family

# Go through each region separately
for (family in families) {
  
  # Subset data for the current family (and change region to All for plotting)
  familyData <- corollaNectarRegion %>% filter(family == !!family)
  familyData$region <- "All"
  
  # Only test this for families with at least 20 species
  if(nrow(familyData) < 7) {
    
    # Create a ggplot for the current region
    plots[[family]] <- ggplot(familyData, aes(x=log10corolla, y=nectar10, shape=region, 
                                              colour=region, fill=region, size=4)) +
      
      theme_bw() + 
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            panel.border = element_rect(colour="black", fill=NA, linewidth=1.5)) +
      
      # set standard plot covering the potential range of all data
      geom_point(data=corollaNectarRegion[!corollaNectarRegion$family == family, ], size=3, pch=20,
                 col=rgb(40, 60, 80, alpha=25, maxColorValue=250)) +
      
      # Add family specific points
      geom_point(bg=familyData$col, col="black", pch=familyData$pch, 
                 alpha=0.7) +
      
      xlab("Corolla length") + 
      ylab("Nectar concentration") +
      theme(legend.title = element_blank(), legend.position="none",
            axis.title = element_text(size=10), 
            axis.text.x = element_text(size=0), 
            axis.text.y = element_text(size=0),
            plot.title = element_text(size=10)) +
      ggtitle(paste0(family, ": Not enough (n<7)"))
    
  } else {
    
    # Family tree
    
    # Model 1
    models[[family]] <- lm(nectar10 ~ log10corolla, 
                            data=familyData)
    
    # Residuals
    residuals[[family]] <- ggplot(data.frame(residuals=models[[family]]$residuals), 
                                  aes(sample=residuals)) +
      geom_qq() +
      geom_qq_line(color="indianred") +
      ggtitle(paste("Q-Q Plot ", family)) +
      theme_bw()
    
    # Create a ggplot for the current region
    plots[[family]] <- ggplot(familyData, aes(x=log10corolla, y=nectar10, shape=region, 
                                              colour=region, fill=region, size=4)) +
      
      theme_bw() + 
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            panel.border = element_rect(colour="black", fill=NA, linewidth=1.5)) +
      
      # set standard plot covering the potential range of all data
      geom_point(data=corollaNectarRegion[!corollaNectarRegion$family == family, ], size=3, pch=20,
                 col=rgb(40, 60, 80, alpha=25, maxColorValu=250)) +
      
      # Now actual region
      geom_point(bg=familyData$col, col="black", pch=familyData$pch, 
                 alpha = 0.7) +
      
      # Model line
        geom_smooth(
          method="glm",
          se = summary(models[[family]])$coefficients[2, ][4] < 0.05,
          show.legend=FALSE,
          lty=ifelse(summary(models[[family]])$coefficients[2, ][4] < 0.05, 1, 2),
          color="black",
          linewidth=1
        )
      
    
    plots[[family]] <- plots[[family]] + scale_fill_manual(values=rep(rgb(40, 60, 80, alpha=125, 
                                                                          maxColorValue=250), 5)) +
      
      xlab("Corolla length") + 
      ylab("Nectar concentration") +
      theme(legend.title = element_blank(), legend.position="none",
            axis.title = element_text(size=10), 
            axis.text.x = element_text(size=0), 
            axis.text.y = element_text(size=0),
            plot.title = element_text(size=10)) +
      
      # Title
      ggtitle(bquote(.(family) * ": " * beta * "=" * .(round(summary(models[[family]])$coefficients[2, 1], 2)) * "," ~
                       R^2 * "=" * .(round(summary(models[[family]])$r.squared, 2)) * ", p=" * .(round(summary(models[[family]])$coefficients[2, 4], 3))))
    
  }
  
}

# Specify the path to the folder
folderPath <- "../figures/4_results/Figure 2/Family"

# Check if the folder exists; if not, create it
if (!dir.exists(folderPath)) {
  dir.create(folderPath, recursive=TRUE)
  message("Folder created at: ", folderPath)
} else {
  message("Folder already exists at: ", folderPath)
}

# Display plots (example to show the first plot)
# Plot and store
# Define the ranges of family indices for separate plots
familyRanges <- list(1:15, 16:30, 31:45, 46:56)

# Set a counter
counter <- 1

# Now let's go through all the webs ans plot them
# But we want to go from south to North (i.e. from Trinidad to Saba)
for(i in seq_along(familyRanges)) {
  
  # Define the current range of families
  current <- familyRanges[[i]]
  
  # Create the file name with the counter prepended
  fileName <- paste0("../figures/4_results/Figure 2/Family/SupplFig. 2a) CorollaNectarScatterplotFamily_NewColors",
                     counter, ".png")
  
  # plot it 
  png(fileName, height=29.7, width=21, units="cm", res=900)
  
  # Arrange the plots for the current family range
  print(plot_grid(plotlist=plots[current], 
                  nrow=5, ncol=3))
  
  # close and save
  dev.off()
  
  # Increment the counter for the next iteration
  counter <- counter + 1
  
}

# Now plot significant models only for supplements
# Initialize an empty vector to store names of significant models
significantModels <- c()

# Loop through each model in the list
for (modelName in names(models)) {
  model <- models[[modelName]]
  # Extract p-values from the model summary
  p <- summary(model)$coefficients[2, ][4]
  
  # Check if any predictor has p < 0.05
  if (any(p < 0.05, na.rm = TRUE)) {
    # If significant, add model name to the vector
    significantModels <- c(significantModels, modelName)
  }
}

# Calculate rows based on the number of significant models
numPlots <- length(significantModels)
ncol <- 3
nrow <- ceiling(numPlots / ncol)

# Create a list of plots for significant models
significantPlots <- lapply(significantModels, function(modelName) plots[[modelName]])

# Define the file name for saving the plot
fileName <- "../figures/4_results/Figure 2/Family/SupplFig. 2a) CorollaNectarScatterplotFamily_Significant.png"

# Plot and save
png(fileName, height=21, width=29, units="cm", res=900)
print(plot_grid(plotlist=significantPlots, nrow=nrow, ncol=ncol))
dev.off()

# Clean working environment
rm(numPlots, nrow, ncol, significantPlots, fileName, i, family, families, current,
   counter, residuals, plots, familyRanges, familyData, corollaNectarRegion,
   model, significantModels, p, modelName, folderPath,
   models)

#####

# 5. Corolla and Color
# 5.1 Data collection corolla length and nectar concentration ####

# First create a sub dataframe, with only the species where we have info on corolla length
# AND color. Cut the data down to the relevant plants
corollaColor <- plants[!is.na(plants$corolla) & !is.na(plants$colorCode), ]

# A bunch of vectors, where I will store information about traits, region etc.
lapply(c("regVec", "specVec", "morphVec", "morphVec2", "famVec"), 
       function(x) assign(x, vector(), envir = .GlobalEnv))

# We create a loop, where we repeat each step for each biogeographical region
for(h in 1:length(bioGeo)){
  
  # A vector with the identities of the networks for the biogeographical region
  # we are currently looking at
  ma <- networkInfo[networkInfo$region == bioGeo[h], ]$network
  
  # Now we will go through each row of the data frame "corolla"
  for(i in 1:nrow(corollaColor)){
    
    # first we check for each row, which networks are listed
    nets <- as.numeric(strsplit(as.character(corollaColor$network[i]),
                                split = ",")[[1]])
    
    # Check if networks 1 or 2 are present in 'nets'
    #if(1 %in% nets || 2 %in% nets){
    #  next # Skip this iteration if networks 1 or 2 are found
    #}
    
    # now we check, if these network(s) are in the specific biogeographical region
    nn <- length(match(ma, nets)[complete.cases(match(networkInfo[networkInfo$region == 
                                                                    bioGeo[h], ]$network, 
                                                      nets))])
    
    # if the network(s) are in the biogeographical region, then collect the data
    if(nn > 0){
      specVec <- append(specVec, corollaColor$species[i]) # plant name
      regVec <- append(regVec, bioGeo[h]) # the region
      morphVec <- append(morphVec, corollaColor$corolla[i]) # corolla
      morphVec2 <- append(morphVec2, corollaColor$colorCode[i]) # color
      famVec <- append(famVec, corollaColor$family[i]) # plant family
    }
    
  } # end i
  
  # Now we join all the information we just collected for one region together
  # in a data frame and then we continue with the next region
  corollaColorRegion <- data.frame(family = famVec,
                                    region = regVec,
                                    corolla = morphVec,
                                    color = morphVec2,
                                    species = specVec)
  
} # end h

# clean a bit
rm(famVec, h, i, ma, morphVec, morphVec2, nets, nn, regVec, specVec)

# Now I adapt the plant names to the format of the phylogenetic tree
corollaColorRegion$species <- gsub(" ", "_", corollaColorRegion$species)

# And finally I log10 transform the corolla length
corollaColorRegion$log10corolla <- log10(corollaColorRegion$corolla + 1)

# Some extra data formatting is needed for the color code
# Data frame where the "clean" color will be stored
corollaColorRegionClean <- data.frame()

# Columns to split or replicate (all columns except "color" could be replicated)
columnsToReplicate <- setdiff(names(corollaColorRegion), "colorCode")

# Now we go through the raw data color data frame to get information for each plant and 
# their colors row by row
for (j in 1:nrow(corollaColorRegion)) {
  
  # Split the color column and find its length for replication
  colorValues <- as.numeric(strsplit(as.character(corollaColorRegion$color[j]), split = ",")[[1]])
  numColors <- length(colorValues)
  
  # Create a list to store replicated data for each column dynamically
  tempData <- list()
  
  # Populate temp_data for each column
  for (col in columnsToReplicate) {
    tempData[[col]] <- rep(corollaColorRegion[[col]][j], numColors)
  }
  
  # Add the "color" column as a factor
  tempData[["color"]] <- as.factor(colorValues)
  
  # Convert temp_data to a data frame for this row
  temp <- as.data.frame(tempData)
  
  # Append the information for each species to the big data frame
  corollaColorRegionClean <- rbind(corollaColorRegionClean, temp)
  
} # end j loop

# Loose the NAs 
corollaColorRegion <- corollaColorRegionClean[complete.cases(corollaColorRegionClean), ]

# clean up a bit
rm(temp, j, col, corollaColorRegionClean, tempData, numColors, colorValues,
   columnsToReplicate)

# throw out any duplicate entries (some species occur in multiple networks in
# the same biogeographical region)
corollaColorRegion <- distinct(corollaColorRegion)

#####

# 5.1.1 Data cleaning - handle duplicate species #####

# Because of the nature of the color data, we now have species with multiple
# entries. In phylogenetic analyses, it is not possible to keep these multiple
# data.
# Go through them and filter
# Prioritize the region Caribbean, as fewest observations there.
# Then prioritize non-ornithophilous and ornithophilous colors (i.e. exclude Two)
speciesFilter <- corollaColorRegion[corollaColorRegion$species %in% 
                                      corollaColorRegion$species[duplicated(corollaColorRegion$species)], ]

# Split by species
checkSpecies <- list()

# Aplly some filters now 
for (species in unique(speciesFilter$species)) {
  # Extract the data frame for the current species
  speciesData <- speciesFilter[speciesFilter$species == species, ]
  
  # Check if there is an entry with region == "Caribbean"
  if ("Caribbean" %in% speciesData$region) {
    # Keep only the "Caribbean" entry
    speciesData <- speciesData[speciesData$region == "Caribbean", ]
  }
  
  # Check if there is an entry with color == 1
  if (1 %in% speciesData$color) {
    # Keep only the entries with color == 1
    speciesData <- speciesData[speciesData$color == 1, ]
  } else if (3 %in% speciesData$color) {
    # If no entry with color == 1, keep those with color == 3
    speciesData <- speciesData[speciesData$color == 3, ]
  }
  
  # Update the list with filtered data for the current species
  checkSpecies[[species]] <- speciesData
  
}

# Convert the list back to a data frame
speciesFilter <- do.call(rbind, checkSpecies)

# Next cut duplicates from original data and then repaste the filtered version
corollaColorRegion <- corollaColorRegion[!corollaColorRegion$species %in% speciesFilter$species, ]
corollaColorRegion <- rbind(corollaColorRegion, speciesFilter)

# rename rows
row.names(corollaColorRegion) <- 1:nrow(corollaColorRegion)

# Clean up
rm(checkSpecies, species, speciesFilter)

# There are still some duplicates, but only entries where a species occurs e.g.
# simultaneously in Andes and Lowland South America
# These entries will be randomly removed when matching to the plantTree

#####

# 5.2 Match data corolla and phylogenetic tree ####

# Read the plant tree
plantTree <- read.tree("4_phylogeny/PlantTree")

# Cut the tree down to the species that we have in our dataset
plantTree <- drop.tip(plantTree, setdiff(plantTree$tip.label,
                                         corollaColorRegion$species))
# How many species were matched
setdiff(unique(corollaColorRegion$species), plantTree$tip.label)

# only keep Caribbean species if duplicated
multipleSpecies <- names(table(corollaColorRegion$species)[table(corollaColorRegion$species) > 1])
# Filter multipleSpecies to keep only those species that occur in the "Caribbean" region
caribbeanSpecies <- multipleSpecies[sapply(multipleSpecies, function(species) {
  "Caribbean" %in% corollaColorRegion[corollaColorRegion$species == species, ]$region
})]
# Now keep only those rows, where the biogeographical region is Caribbean
corollaColorRegion <- corollaColorRegion[!(corollaColorRegion$species %in% caribbeanSpecies & 
                                             corollaColorRegion$region != "Caribbean"), ]

# then we need to sort the trait data so that the species follow the order of the tree
corollaColorRegion <- corollaColorRegion[match(plantTree$tip.label,
                                               corollaColorRegion$species), ]

# Clean working environment
rm(caribbeanSpecies, multipleSpecies)

#####

# 5.3 Corolla Color phylogenetic model ####

# Prepare data
# Remove node labels by setting them to NULL
plantTree$node.label <- NULL
# Phylogenetic object
phyloData <- comparative.data(phy=plantTree, 
                              data=corollaColorRegion, 
                              names.col="species") 
# Phyl Model 1
phylModelCorolla <- phylANOVA(plantTree, corollaColorRegion$color, 
                       corollaColorRegion$log10corolla, 
                       p.adj="bonferroni",
                       nsim=100000)

# Generate ANOVA results section
resultsCorolla <- sprintf(
  "ANOVA table: Phylogenetic ANOVA\n\nResponse: Corolla length\n%-15s %-10s %-10s\n%-15s %-10.5f %-10.5f\n%-15s %-10.5f %-10.5f\n\nF value: %.2f\nP value: %.3f\n\n",
  "Source", "Sum Sq", "Mean Sq",
  "Region", phylModelCorolla$SumSq[1], phylModelCorolla$MeanSq[1],
  "Residual", phylModelCorolla$SumSq[2], phylModelCorolla$MeanSq[2],
  phylModelCorolla$F, phylModelCorolla$Pf
)

# Format post hoc results with pairwise comparisons and p-values
posthocCorolla <- "Pairwise posthoc test using method = 'bonferroni'\n\nComparison                P-value\n"
comparisons <- combn(colnames(phylModelCorolla$Pt), 2)
pValues <- mapply(function(i, j) phylModelCorolla$Pt[i, j], comparisons[1,], comparisons[2,])
comparisonNames <- apply(comparisons, 2, function(pair) paste(pair, collapse = "-"))
names(pValues) <- comparisonNames
posthocCorolla <- paste(posthocCorolla, paste(sprintf("%-24s %-10.3f", comparisonNames, pValues), collapse = "\n"))

# Generate letters for significant groupings
lettersCorolla <- multcompLetters(pValues, threshold = 0.05)$Letters
posthocCorolla <- paste(posthocCorolla, "\n\nSignificance Letters:\n\n", paste(sprintf("%-24s %s", names(lettersCorolla), lettersCorolla), collapse = "\n"))

# Combine results
outputCorolla <- paste(resultsCorolla, posthocCorolla)

# Write output to file
writeLines(outputCorolla, "../results/textFiles/Figure 2/phylANOVA(CorollaColor).txt")

# Test for homogeneity of variance
homogeneityCorolla <- leveneTest(log10corolla ~ color, data=corollaColorRegion)$'Pr(>F)'[1]

# Summary generation
summaryCorolla <- sprintf(
  "Corolla Length\n\nPhylogenetic ANOVA Results:\nF = %.2f\np = %.4f\n\nHomogeneity of variance: %s",
  phylModelCorolla$F,
  phylModelCorolla$Pf,
  ifelse(homogeneityCorolla < 0.05, "is not given", "is given")
)

# Now let's run a few different tests to confirm.
# log10corolla is not normal, and Homogeneity of variance is not given
# so let's just run a Kruskal Wallis
kT <- kruskal.test(log10corolla ~ color, data=corollaColorRegion)
# then Welch's Anova
wA <- oneway.test(log10corolla ~ color, data=corollaColorRegion, var.equal=FALSE)

# Append Kruskal-Wallis test results
summaryCorolla <- paste(
  summaryCorolla,
  sprintf(
    "\n\nKruskal-Wallis Test Results:\nChi-squared = %.2f\nDegrees of Freedom = %d\np-value = %.4f",
    kT$statistic, kT$parameter, kT$p.value
  ),
  sprintf(
    "\n\nWelch Anova Results:\nF = %.2f\nDegrees of Freedom = %d\np-value = %.4f",
    wA$statistic, wA$parameter[1], wA$p.value
  )
)

# Output summary (e.g., print, save, etc.)
writeLines(summaryCorolla, "../results/textFiles/Figure 2/summaryColorCorolla.txt") 

# Clean up
rm(summaryCorolla, resultsCorolla, pValues, posthocCorolla, outputCorolla,
   homogeneityCorolla, comparisonNames, comparisons, kT, phyloData, wA,
   phylModelCorolla, speciesData)

#####

# 5.4 Corolla Color Violinplot ####

# First I create some graphics parameters, pch and color
# Add some plotting parameters: pch
corollaColorRegion$pch <- ifelse(corollaColorRegion$region == "Lowland South America", 21,
                                  ifelse(corollaColorRegion$region == "Andes", 22,
                                         ifelse(corollaColorRegion$region == "Caribbean", 23, 24)))

# Add some plotting parameters: color
corollaColorRegion$col <- ifelse(corollaColorRegion$region == "Lowland South America", "beige",
                                  ifelse(corollaColorRegion$region == "Andes", "bisque4",
                                         ifelse(corollaColorRegion$region == "Caribbean", "cornflowerblue", "#f1a35f")))

# Add a new label for the axis of the plot
corollaColorRegion$axis <- ifelse(corollaColorRegion$color == "1", "One",
                                 ifelse(corollaColorRegion$color == "2", "Two", "Three"))

# Letters as data frame
# Create a data frame with the regions and their corresponding significance letters
lettersCorolla <- data.frame(
  color = c(names(lettersCorolla)),
  label = c(lettersCorolla)
)

# Arrange as will appear on X-Axis
# Convert 'region' column to a factor with specified levels
lettersCorolla$color <- factor(lettersCorolla$color)
lettersCorolla <- lettersCorolla %>% arrange(color)
lettersCorolla$axis <- c("One", "Two", "Three")

# Plot and store
figureManuscript[[2]] <- ggplot(corollaColorRegion, aes(x=axis, y=log10corolla, fill=axis)) +
  
  geom_violin(trim=FALSE, alpha=0.4) + 
  geom_jitter(width=0.1, size=4, pch=corollaColorRegion$pch, 
              bg=corollaColorRegion$col) +  
  geom_boxplot(width=0.2, fill=adjustcolor("ghostwhite", alpha.f = 0.7), 
               outlier.shape=NA) + 
  
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border=element_rect(colour="black", fill=NA, size=1.5)) +
   
  scale_x_discrete(limits=c("One", "Two", "Three")) +
  scale_fill_manual(values=c("ghostwhite", "#ffd8b1", "indianred3"), 
                    limits=c("One", "Two", "Three"),  
                    name="") +
  
  scale_y_continuous(
    breaks = seq(0, max(corollaColorRegion$log10corolla), by = 0.5), # Adjust as needed
    labels = scales::label_number(accuracy = 0.01) # Ensure 2 decimal places
  ) +
  
  labs(x="Floral color category", y="Log corolla length (mm)") +
  
  theme(legend.position="none", 
        legend.margin=margin(0, 0, 0, 0), 
        axis.title=element_text(size=30), 
        axis.text.x=element_text(size=20), 
        axis.text.y=element_text(size=20, angle=90),
        plot.title=element_text(size=30)) +
  
  # Add the significance letters
  geom_text(data=lettersCorolla, aes(x=axis, y=max(corollaColorRegion$log10corolla+0.5),
                                     label=label), 
            size=8, vjust=1) 

# print it
png("../figures/4_results/Figure 2/Fig. 2b) CorollaColorViolinplot_NewColors_Ylab.png", 
    height=21, width=29.7, units="cm", res=900)
print(figureManuscript[[2]])
dev.off()

#####

# 5.4.1 Corolla Color Violinplot - Region ####

# Create a list to store models and plots
models <- list()
plotsCorollaColor <- list()

# Go through each region separately
# Go through each region separately
for (i in seq_along(rev(bioGeo))) {
  region <- rev(bioGeo)[i]
  
  # Subset data for the current region
  regionData <- corollaColorRegion %>% filter(region == !!region)
  
  # Region tree
  # Cut the tree down to the species that we have in our dataset
  regionTree <- drop.tip(plantTree, setdiff(plantTree$tip.label,
                                            regionData$species))
  # How many species were matched
  setdiff(unique(regionData$species), regionTree$tip.label)
  
  # then we need to sort the trait data so that the species follow the order of the tree
  regionData <- regionData[match(regionTree$tip.label,
                                 regionData$species), ]
  
  # PGLS Model 1
  models[[region]] <- phylANOVA(regionTree, regionData$color, 
                                regionData$log10corolla, 
                                p.adj="bonferroni",
                                nsim=100000)
  
  # Extract the post hoc results
  posthocColorP <- models[[region]]$Pt
  
  # Create a vector to store p-values for pairwise comparisons
  pValues <- c()
  comparisons <- c()
  
  # Loop through pairwise comparisons to add to the string
  for (k in 1:nrow(posthocColorP)) {
    for (j in 1:ncol(posthocColorP)) {
      if (k < j) {  # To avoid duplicate comparisons and self-comparisons
        comparison <- paste(rownames(posthocColorP)[k], "vs", colnames(posthocColorP)[j])
        pValue <- posthocColorP[k, j]
        
        # Append p-values and comparison names
        pValues <- c(pValues, pValue)
        comparisons <- c(comparisons, comparison)
      }
    }
  }
  
  # Reformat comparison names to contain only a single hyphen
  comparisons <- gsub(" vs ", "-", comparisons)
  
  # Assign the reformatted names to the p-values vector
  names(pValues) <- comparisons
  
  # Generate the letters to represent significant groupings
  lettersColor <- multcompLetters(pValues, threshold=0.05)$Letters
  
  # Create a data frame with the regions and their corresponding significance letters
  lettersColor <- data.frame(
    color = c(names(lettersColor)),
    label = c(lettersColor)
  )
  
  # Arrange as will appear on X-Axis
  # Convert 'region' column to a factor with specified levels
  lettersColor$color <- factor(lettersColor$color)
  lettersColor <- lettersColor %>% arrange(color)
  lettersColor$axis <- c("One", "Two", "Three")
  
  # Create a ggplot for the current region
  plotsCorollaColor[[region]] <- ggplot(corollaColorRegion, aes(x=axis, y=log10corolla, fill=axis)) +
    
    # Transparent base layer to set plot dimensions without displaying data
    geom_blank() + 
    scale_y_continuous(limits=c(min(corollaColorRegion$log10corolla)-0.6,
                                max(corollaColorRegion$log10corolla)+0.6)) +
    
    # Add violin plot from `regionData`
    geom_violin(data=regionData, aes(x=axis, y=log10corolla, fill=axis), 
                trim=FALSE, alpha=0.4) +
    
    # Add jitter layer from `colorRegion` in grey
    geom_jitter(data=corollaColorRegion[!corollaColorRegion$region == region, ],
                aes(x=axis, y=log10corolla), 
                bg="grey", alpha=0.5, width=0.1, size=4, 
                pch=corollaColorRegion[!corollaColorRegion$region == region, ]$pch) +
    
    # Add jitter layer from `regionData` with original color and shape settings
    geom_jitter(data=regionData, aes(x=axis, y=log10corolla), 
                width=0.1, size=4, pch=regionData$pch, 
                bg=regionData$col) +
    
    # Boxplot from `regionData` to show central tendency and spread
    geom_boxplot(data=regionData, aes(x=axis, y=log10corolla), 
                 width=0.2, fill=adjustcolor("ghostwhite", alpha.f=0.7), outlier.shape=NA) +  
    
    coord_flip() +
    
    theme_bw() +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          panel.border=element_rect(colour="black", fill=NA, 
                                    size=1.5)) +
    
    scale_x_discrete(limits=c("One", "Two", "Three")) +
    scale_fill_manual(values=c("ghostwhite", "#ffd8b1", "indianred3"), 
                      limits=c("One", "Two", "Three"),  
                      name="") + 
    
    xlab("Color category") +
    ylab(ifelse(i == length(bioGeo), "Corolla length", "")) +
    theme(legend.position="none", 
          legend.margin=margin(0, 0, 0, 0), 
          axis.title=element_text(size=18), 
          axis.text.x=element_text(size=0), 
          axis.text.y=element_text(size=0),
          plot.title=element_text(size=18)) +
    
    # Add the significance letters
    geom_text(data=lettersColor, aes(x=axis, y=max(corollaColorRegion$log10corolla+0.6),
                                     label=label), size=8, vjust=1) +
    
    ggtitle(bquote(.(region) * ": " * "F=" * .(round(models[[region]]$F, 2)) * "," * 
                     " p=" * .(round(models[[region]]$Pf, 3))))
  
}

# Plot and store
png("../figures/4_results/Figure 2/Region/SupplFig. 2b) CorollaColorViolinplotRegion_NewColors.png", 
    height=29, width=18, units="cm", res=900)
# Arrange the plots in a grid
plot_grid(plotsCorollaColor[[1]], plotsCorollaColor[[2]], 
          plotsCorollaColor[[3]], plotsCorollaColor[[4]], 
          nrow=4, ncol=1)
dev.off()

# Clean up
rm(region, pValues, pValue, k, j, i, comparison, comparisons,
   regionTree, regionData, posthocColorP, models, plotsCorollaColor,
   lettersColor)

#####

# 5.4.2 Corolla Color Violinplot - Family ####

# Create a list to store models and plots
models <- list()
plots <- list()

# Count species per family and order by species number
families <- as.data.frame(corollaColorRegion %>%
                            group_by(family) %>%
                            summarise(species_count = n()) %>%
                            arrange(desc(species_count)))$family

# Go through each family separately
for (family in families) {
  
  # Subset data for the current family (and change region to All for plotting)
  familyData <- corollaColorRegion %>% filter(family == !!family)
  familyData$region <- "All"
  
  # Only test this for families with at least 20 species
  if(nrow(familyData) < 7) {
    
    # Create a ggplot for the current region
    plots[[family]] <- ggplot(corollaColorRegion, aes(x=axis, y=log10corolla)) +
      
      # Transparent base layer to set plot dimensions without displaying data
      geom_blank() + 
      scale_y_continuous(limits=c(min(corollaColorRegion$log10corolla)-0.6,
                                  max(corollaColorRegion$log10corolla)+0.6)) +
      
      # Add violin plot from `regionData`
      geom_violin(data=familyData, aes(x=axis, y=log10corolla, fill=axis), 
                  trim=FALSE, alpha=0.4) +
      
      # Add jitter layer from `colorRegion` in grey
      geom_jitter(data=corollaColorRegion[!corollaColorRegion$family == family, ], 
                  aes(x=axis, y=log10corolla), 
                  bg="grey", alpha=0.5, width=0.1, size=4, 
                  pch=corollaColorRegion[!corollaColorRegion$family == family, ]$pch) +
      
      # Add jitter layer from `regionData` with original color and shape settings
      geom_jitter(data=familyData, aes(x=axis, y=log10corolla), 
                  width=0.1, size=4, pch=familyData$pch, 
                  bg=familyData$col) +
      
      # Boxplot from `regionData` to show central tendency and spread
      geom_boxplot(data=familyData, aes(x=axis, y=log10corolla), 
                   width=0.2, fill=adjustcolor("ghostwhite", alpha.f=0.7), outlier.shape=NA) +  
      
      coord_flip() +
      
      theme_bw() +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            panel.border=element_rect(colour="black", fill=NA, size=1.5)) +
      
      scale_x_discrete(limits=c("One", "Two", "Three")) +
      scale_fill_manual(values=c("ghostwhite", "#ffd8b1", "indianred3"), 
                        limits=c("One", "Two", "Three"),  
                        name="") + 
      
      xlab("Color category") + 
      ylab("Corolla length") +
      theme(legend.title = element_blank(), legend.position="none",
            axis.title = element_text(size=10), 
            axis.text.x = element_text(size=0), 
            axis.text.y = element_text(size=0),
            plot.title = element_text(size=10)) +
      ggtitle(paste0(family, ": Not enough (n<7)"))
    
  } else if (length(unique(familyData$axis)) < 2) {
    
    # Create a ggplot for the current region
    plots[[family]] <- ggplot(corollaColorRegion, aes(x=axis, y=log10corolla)) +
      
      # Transparent base layer to set plot dimensions without displaying data
      geom_blank() + 
      scale_y_continuous(limits=c(min(corollaColorRegion$log10corolla)-0.6,
                                  max(corollaColorRegion$log10corolla)+0.6)) +
      
      # Add violin plot from `regionData`
      geom_violin(data=familyData, aes(x=axis, y=log10corolla, fill=axis), 
                  trim=FALSE, alpha=0.4) +
      
      # Add jitter layer from `colorRegion` in grey
      geom_jitter(data=corollaColorRegion[!corollaColorRegion$family == family, ], 
                  aes(x=axis, y=log10corolla), 
                  bg="grey", alpha=0.5, width=0.1, size=4, 
                  pch=corollaColorRegion[!corollaColorRegion$family == family, ]$pch) +
      
      # Add jitter layer from `regionData` with original color and shape settings
      geom_jitter(data=familyData, aes(x=axis, y=log10corolla), 
                  width=0.1, size=4, pch=familyData$pch, 
                  bg=familyData$col) +
      
      # Boxplot from `regionData` to show central tendency and spread
      geom_boxplot(data=familyData, aes(x=axis, y=log10corolla), 
                   width=0.2, fill=adjustcolor("ghostwhite", alpha.f=0.7), outlier.shape=NA) +  
      
      coord_flip() +
      
      theme_bw() +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            panel.border=element_rect(colour="black", fill=NA, size=1.5)) +
      
      scale_x_discrete(limits=c("One", "Two", "Three")) +
      scale_fill_manual(values=c("ghostwhite", "#ffd8b1", "indianred3"), 
                        limits=c("One", "Two", "Three"),  
                        name="") + 
      
      xlab("Color category") + 
      ylab("Corolla length") +
      theme(legend.title = element_blank(), legend.position="none",
            axis.title = element_text(size=10), 
            axis.text.x = element_text(size=0), 
            axis.text.y = element_text(size=0),
            plot.title = element_text(size=10)) +
      ggtitle(paste0(family, ": Only one color"))
    
  } else {
    
    # Model 1
    models[[family]] <- glm(log10corolla ~ color, 
                            data=familyData)
    
    # Letters
    ptt2 <- aov(log10corolla ~ color, familyData)
    lettersColor <- data.frame("Letters" = multcompLetters(extract_p(TukeyHSD(ptt2)$"color"))$"Letters")
    
    # Create a data frame with the regions and their corresponding significance letters
    lettersColor <- data.frame(
      color = c(row.names(lettersColor)),
      label = c(lettersColor)
    )
    
    # Define the complete set of colors you expect
    expectedColors <- 1:3
    
    # Check for missing colors
    missingColors <- setdiff(expectedColors, lettersColor$color)
    
    # If there are missing colors, add them with "na" in the Letters column
    if (length(missingColors) > 0) {
      # Create a dataframe with the missing colors and "na" in the Letters column
      missingDf <- data.frame(color=missingColors, Letters = rep("na", length(missingColors)))
      
      # Bind the missing colors to the original dataframe
      lettersColor <- rbind(lettersColor, missingDf)
    }
    
    # Arrange as will appear on X-Axis
    # Convert 'region' column to a factor with specified levels
    lettersColor$color <- factor(lettersColor$color)
    lettersColor <- lettersColor %>% arrange(color)
    lettersColor$axis <- c("One", "Two", "Three")
    
    # Create a ggplot for the current region
    plots[[family]] <- ggplot(corollaColorRegion, aes(x=axis, y=log10corolla)) +
      
      # Transparent base layer to set plot dimensions without displaying data
      geom_blank() + 
      scale_y_continuous(limits=c(min(corollaColorRegion$log10corolla)-0.6,
                                  max(corollaColorRegion$log10corolla)+0.6)) +
      
      # Add violin plot from `regionData`
      geom_violin(data=familyData, aes(x=axis, y=log10corolla, fill=axis), 
                  trim=FALSE, alpha=0.4) +
      
      # Add jitter layer from `colorRegion` in grey
      geom_jitter(data=corollaColorRegion[!corollaColorRegion$family == family, ], aes(x=axis, y=log10corolla), 
                  bg="grey", alpha=0.5, width=0.1, size=4, 
                  pch=corollaColorRegion[!corollaColorRegion$family == family, ]$pch) +
      
      # Add jitter layer from `regionData` with original color and shape settings
      geom_jitter(data=familyData, aes(x=axis, y=log10corolla), 
                  width=0.1, size=4, pch=familyData$pch, 
                  bg=familyData$col) +
      
      # Boxplot from `regionData` to show central tendency and spread
      geom_boxplot(data=familyData, aes(x=axis, y=log10corolla), 
                   width=0.2, fill=adjustcolor("ghostwhite", alpha.f=0.7), 
                   outlier.shape=NA) +  
      
      coord_flip() +
      
      theme_bw() +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            panel.border=element_rect(colour="black", fill=NA, size=1.5)) +
      
      scale_x_discrete(limits=c("One", "Two", "Three")) +
      scale_fill_manual(values=c("ghostwhite", "#ffd8b1", "indianred3"), 
                        limits=c("One", "Two", "Three"),  
                        name="") + 
      
      xlab("Color category") +
      ylab("Corolla length") +
      theme(legend.position="none", 
            legend.margin=margin(0, 0, 0, 0), 
            axis.title=element_text(size=10), 
            axis.text.x=element_text(size=0), 
            axis.text.y=element_text(size=0),
            plot.title=element_text(size=10)) +
      
      # Add the significance letters
      geom_text(data=lettersColor, aes(x=axis, y=max(corollaColorRegion$log10corolla+0.6),
                                       label=Letters), size=8, vjust=1) +
      
      ggtitle(bquote(.(family) * ": " * chi^2 * "=" * .(round(Anova(models[[family]], type = c("II","III", 2, 3), 
                                                                    vcov.= vcov(models[[family]], complete=TRUE))[[1]], 2)) * 
                       "," ~ "p=" * .(round(Anova(models[[family]], type = c("II","III", 2, 3), 
                                                  vcov.= vcov(models[[family]], complete=TRUE))[[3]], 3))))
    
  }
  
}

# Display plots (example to show the first plot)
# Plot and store
# Define the ranges of family indices for separate plots
familyRanges <- list(1:15, 16:30, 31:45, 46:60, 61:65)

# Set a counter
counter <- 1

# Now let's go through all the webs ans plot them
# But we want to go from south to North (i.e. from Trinidad to Saba)
for(i in seq_along(familyRanges)) {
  
  # Define the current range of families
  current <- familyRanges[[i]]
  
  # Create the file name with the counter prepended
  fileName <- paste0("../figures/4_results/Figure 2/Family/SupplFig. 2b) CorollaColorViolinplotFamily_NewColors",
                     counter, ".png")
  
  # plot it 
  png(fileName, height=29.7, width=21, units="cm", res=900)
  
  # Arrange the plots for the current family range
  print(plot_grid(plotlist=plots[current], 
                  nrow=5, ncol=3))
  
  # close and save
  dev.off()
  
  # Increment the counter for the next iteration
  counter <- counter + 1
  
}

# Now plot significant models only for supplements
# Initialize an empty vector to store names of significant models
significantModels <- c()

# Loop through each model in the list
for (modelName in names(models)) {
  model <- models[[modelName]]
  # Extract p-values from the model summary
  p <- round(Anova(model, type = c("II","III", 2, 3), 
                   vcov.= vcov(model, complete=TRUE))[[3]], 3)
  
  # Check if any predictor has p < 0.05
  if (any(p < 0.05, na.rm = TRUE)) {
    # If significant, add model name to the vector
    significantModels <- c(significantModels, modelName)
  }
}

# Calculate rows based on the number of significant models
numPlots <- length(significantModels)
ncol <- 3
nrow <- ceiling(numPlots / ncol)

# Create a list of plots for significant models
significantPlots <- lapply(significantModels, function(modelName) plots[[modelName]])

# Define the file name for saving the plot
fileName <- "../figures/4_results/Figure 2/Family/SupplFig. 2b) CorollaColorViolinplotFamily_Significant_NewColors.png"

# Plot and save
png(fileName, height=21, width=29, units="cm", res=900)
print(plot_grid(plotlist=significantPlots, nrow=nrow, ncol=ncol))
dev.off()

# Clean working environment
rm(significantModels, pValue, pValues, p, numPlots, nrow, ncol, modelName, missingColors,
   k, j, i, fileName, family, families, expectedColors, current, counter, comparisons,
   comparison, significantPlots, ptt2, model, missingDf, lettersColor, familyRanges,
   familyData, models, plots, plantTree, lettersCorolla, corollaColor, corollaColorRegion)

#####

# 6. Nectar and Color
# 6.1 Data collection nectar concentration ####

# Dataframe, with nectar concentration AND color 
nectarColor <- plants[!is.na(plants$nectar) & !is.na(plants$colorCode), ]
nectarColor <- nectarColor[!nectarColor$nectar < 5, ]

# A bunch of vectors, where I will store information about traits, region etc.
lapply(c("regVec", "specVec", "morphVec", "morphVec2", "famVec"), 
       function(x) assign(x, vector(), envir = .GlobalEnv))

# Function to collect data for a specific region
collectRegionData <- function(region, regionNetworks, data) {
  for (i in seq_len(nrow(data))) {
    nets <- as.numeric(strsplit(as.character(data$network[i]), split = ",")[[1]])
    
    # Check if any network in the row matches the networks for the region
    if (any(nets %in% regionNetworks)) {
      specVec <<- append(specVec, data$species[i])
      regVec <<- append(regVec, region)
      morphVec <<- append(morphVec, data$nectar[i])
      morphVec2 <<- append(morphVec2, data$colorCode[i])
      famVec <<- append(famVec, data$family[i])
    }
  }
}

# Loop through each biogeographical region
for (region in bioGeo) {
  regionNetworks <- networkInfo$network[networkInfo$region == region]
  collectRegionData(region, regionNetworks, nectarColor)
}

# Combine collected data into a single data frame
nectarColorRegion <- data.frame(
  family = famVec,
  region = regVec,
  nectar = morphVec,
  color = morphVec2,
  species = specVec
)

# Adjust species names and scale nectar concentration
nectarColorRegion$species <- gsub(" ", "_", nectarColorRegion$species)
nectarColorRegion$nectar10 <- nectarColorRegion$nectar / 100

# clean
nectarColorRegion <- as.data.frame(nectarColorRegion %>%
  separate_rows(color, sep = ",") %>%
  mutate(color = as.factor(as.numeric(color))) %>%
  distinct(family, region, color, species, .keep_all = TRUE))

# clean up a bit
rm(specVec, regVec, regionNetworks, region, morphVec2, morphVec,
   famVec)

#####

# 6.1.1 Data cleaning - handle duplicate species #####

# Because of the nature of the color data, we now have species with multiple
# entries. In phylogenetic analyses, it is not possible to keep these multiple
# data.
# Go through them and filter
# Prioritize the region Caribbean, as fewest observations there.
# Then prioritize non-ornithophilous and ornithophilous colors (i.e. exclude Two)
speciesFilter <- nectarColorRegion[nectarColorRegion$species %in% 
                                     nectarColorRegion$species[duplicated(nectarColorRegion$species)], ]

# Split by species
checkSpecies <- list()

# Aplly some filters now 
for (species in unique(speciesFilter$species)) {
  # Extract the data frame for the current species
  speciesData <- speciesFilter[speciesFilter$species == species, ]
  
  # Check if there is an entry with region == "Caribbean"
  if ("Caribbean" %in% speciesData$region) {
    # Keep only the "Caribbean" entry
    speciesData <- speciesData[speciesData$region == "Caribbean", ]
  }
  
  # Check if there is an entry with color == 1
  if (1 %in% speciesData$color) {
    # Keep only the entries with color == 1
    speciesData <- speciesData[speciesData$color == 1, ]
  } else if (3 %in% speciesData$color) {
    # If no entry with color == 1, keep those with color == 3
    speciesData <- speciesData[speciesData$color == 3, ]
  }
  
  # Update the list with filtered data for the current species
  checkSpecies[[species]] <- speciesData
  
}

# Convert the list back to a data frame
speciesFilter <- do.call(rbind, checkSpecies)

# Next cut duplicates from original data and then repaste the filtered version
nectarColorRegion <- nectarColorRegion[!nectarColorRegion$species %in% speciesFilter$species, ]
nectarColorRegion <- rbind(nectarColorRegion, speciesFilter)

# rename rows
row.names(nectarColorRegion) <- 1:nrow(nectarColorRegion)

# Clean up
rm(checkSpecies, species, speciesData, speciesFilter)

# There are still some duplicates, but only entries where a species occurs e.g.
# simultaneously in Andes and Lowland South America
# These entries will be randomly removed when matching to the plantTree

#####

# 6.2 Match data nectar and phylogenetic tree ####

# Read the plant tree
plantTree <- read.tree("4_phylogeny/PlantTree")

# Cut the tree down to the species that we have in our dataset
plantTree <- drop.tip(plantTree, setdiff(plantTree$tip.label,
                                         nectarColorRegion$species))
# How many species were matched
setdiff(unique(nectarColorRegion$species), plantTree$tip.label)

# Prioritize Caribbean 
multipleSpecies <- names(table(nectarColorRegion$species)[table(nectarColorRegion$species) > 1])
# Filter multipleSpecies to keep only those species that occur in the "Caribbean" region
caribbeanSpecies <- multipleSpecies[sapply(multipleSpecies, function(species) {
  "Caribbean" %in% nectarColorRegion[nectarColorRegion$species == species, ]$region
})]
# Now keep only those rows, where the biogeographical region is Caribbean
nectarColorRegion <- nectarColorRegion[!(nectarColorRegion$species %in% caribbeanSpecies & 
                                           nectarColorRegion$region != "Caribbean"), ]

# then we need to sort the trait data so that the species follow the order of the tree
nectarColorRegion <- nectarColorRegion[match(plantTree$tip.label,
                                        nectarColorRegion$species), ]

# Clean working environment
rm(caribbeanSpecies, multipleSpecies)

#####

# 6.3 Nectar color phylogenetic model ####

# Prepare data
# Remove node labels by setting them to NULL
plantTree$node.label <- NULL
# Phylogenetic object
phyloData <- comparative.data(phy=plantTree, 
                              data=nectarColorRegion, 
                              names.col="species") 
# Phyl Model 1
phylModelNectar <- phylANOVA(plantTree, nectarColorRegion$color, 
                       nectarColorRegion$nectar10, 
                       p.adj="bonferroni",
                       nsim=100000)

# Generate ANOVA results section for Nectar
resultsNectar <- sprintf(
  "ANOVA table: Phylogenetic ANOVA\n\nResponse: Nectar concentration\n%-15s %-10s %-10s\n%-15s %-10.5f %-10.5f\n%-15s %-10.5f %-10.5f\n\nF value: %.2f\nP value: %.3f\n\n",
  "Source", "Sum Sq", "Mean Sq",
  "Region", phylModelNectar$`Sum Sq`[1], phylModelNectar$`Mean Sq`[1],
  "Residual", phylModelNectar$`Sum Sq`[2], phylModelNectar$`Mean Sq`[2],
  phylModelNectar$F, phylModelNectar$Pf
)

# Format post hoc results with pairwise comparisons and p-values for Nectar
posthocNectar <- "Pairwise posthoc test using method = 'bonferroni'\n\nComparison                P-value\n"
comparisons <- combn(colnames(phylModelNectar$Pt), 2)
pValues <- mapply(function(i, j) phylModelNectar$Pt[i, j], comparisons[1,], comparisons[2,])
comparisonNames <- apply(comparisons, 2, function(pair) paste(pair, collapse = "-"))
names(pValues) <- comparisonNames
posthocNectar <- paste(posthocNectar, paste(sprintf("%-24s %-10.3f", comparisonNames, pValues), collapse = "\n"))

# Generate letters for significant groupings for Nectar
lettersNectar <- multcompLetters(pValues, threshold = 0.05)$Letters
posthocNectar <- paste(posthocNectar, "\n\nSignificance Letters:\n\n", paste(sprintf("%-24s %s", names(lettersNectar), lettersNectar), collapse = "\n"))

# Combine results for Nectar
outputNectar <- paste(resultsNectar, posthocNectar)

# Write output to file for Nectar
writeLines(outputNectar, "../results/textFiles/Figure 2/phylANOVA(NectarColor).txt")

# Test for homogeneity of variance for Nectar
homogeneityNectar <- leveneTest(nectar10 ~ color, data=nectarColorRegion)$`Pr(>F)`[1]

# Summary generation for Nectar
summaryNectar <- sprintf(
  "Nectar Concentration\n\nPhylogenetic ANOVA Results:\nF = %.2f\np = %.4f\n\nHomogeneity of variance: %s",
  phylModelNectar$F,
  phylModelNectar$Pf,
  ifelse(homogeneityNectar < 0.05, "is not given", "is given")
)

# Now let's run a few different tests to confirm for Nectar.
# log10nectar is not normal, and Homogeneity of variance is not given
# so let's just run a Kruskal Wallis for Nectar
kT <- kruskal.test(nectar10 ~ color, data=nectarColorRegion)
# then Welch's Anova for Nectar
wA <- oneway.test(nectar10 ~ color, data=nectarColorRegion, var.equal=FALSE)

# Append Kruskal-Wallis test results for Nectar
summaryNectar <- paste(
  summaryNectar,
  sprintf(
    "\n\nKruskal-Wallis Test Results:\nChi-squared = %.2f\nDegrees of Freedom = %d\np-value = %.4f",
    kT$statistic, kT$parameter, kT$p.value
  ),
  sprintf(
    "\n\nWelch Anova Results:\nF = %.2f\nDegrees of Freedom = %d\np-value = %.4f",
    wA$statistic, wA$parameter[1], wA$p.value
  )
)

# Output summary for Nectar (e.g., print, save, etc.)
writeLines(summaryNectar, "../results/textFiles/Figure 2/summaryColorNectar.txt")

# Clean up
rm(summaryNectar, resultsNectar, pValues, posthocNectar, outputNectar,
   homogeneityNectar, comparisonNames, comparisons, kT, phyloData, wA,
   phylModelNectar, collectRegionData)

#####

# 6.4 Nectar color Violinplot ####

# First I create some graphics parameters, pch and color
# Add some plotting parameters: pch
nectarColorRegion$pch <- ifelse(nectarColorRegion$region == "Lowland South America", 21,
                                 ifelse(nectarColorRegion$region == "Andes", 22,
                                        ifelse(nectarColorRegion$region == "Caribbean", 23, 24)))

# Add some plotting parameters: color
nectarColorRegion$col <- ifelse(nectarColorRegion$region == "Lowland South America", "beige",
                                 ifelse(nectarColorRegion$region == "Andes", "bisque4",
                                        ifelse(nectarColorRegion$region == "Caribbean", "cornflowerblue", "#f1a35f")))

# Add a new label for the axis of the plot
nectarColorRegion$axis <- ifelse(nectarColorRegion$color == "1", "One",
                                  ifelse(nectarColorRegion$color == "2", "Two", "Three"))

# Letters as data frame
# Create a data frame with the regions and their corresponding significance letters
lettersNectar <- data.frame(
  color = c(names(lettersNectar)),
  label = c(lettersNectar)
)

# Arrange as will appear on X-Axis
# Convert 'region' column to a factor with specified levels
lettersNectar$color <- factor(lettersNectar$color)
lettersNectar <- lettersNectar %>% arrange(color)
lettersNectar$axis <- c("One", "Two", "Three")

# Create the figure now
figureManuscript[[3]] <- ggplot(nectarColorRegion, aes(x=axis, y=nectar10, fill=axis)) +
  
  geom_violin(trim=FALSE, alpha=0.4) + 
  
  geom_jitter(width=0.1, size=4, pch=nectarColorRegion$pch, 
              bg=nectarColorRegion$col) +  
  geom_boxplot(width=0.2, fill=adjustcolor("ghostwhite", alpha.f = 0.7), 
               outlier.shape=NA) + 
  
  coord_flip() +
  
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border=element_rect(colour="black", fill=NA, size=1.5)) +
  
  scale_x_discrete(limits=c("One", "Two", "Three")) +
  scale_fill_manual(values=c("ghostwhite", "#ffd8b1", "indianred3"), 
                    limits=c("One", "Two", "Three"),  
                    name="") +
  
  scale_y_continuous(
    breaks = seq(0, max(nectarColorRegion$nectar10 + 0.2), by = 0.1), 
    labels = scales::label_number(accuracy = 0.01, scale = 100) # Multiply by 100, no "%"
  ) +
  
  labs(x="Floral color category", y="Nectar concentration (%)") +
  theme(legend.position="none", 
        legend.margin=margin(0, 0, 0, 0), 
        axis.title=element_text(size=30), 
        axis.text.x=element_text(size=20), 
        axis.text.y=element_text(size=20, angle=90, hjust=0.5),
        plot.title=element_text(size=0)) +

    # Add the significance letters
    geom_text(data=lettersNectar, aes(x=axis, y=max(nectarColorRegion$nectar10+0.1),
                                   label=label), 
          size=8, vjust=1) 

# Plot and store
png("../figures/4_results/Figure 2/Fig. 2c) NectarColorViolinplot_NewColors_Ylab.png", height = 21, width = 29.7, 
    units = "cm", res = 900)
print(figureManuscript[[3]])
dev.off()

#####

# 6.4.1 Nectar Color Violinplot - Region ####

# Create a list to store models and plots
models <- list()
plotsNectarColor <- list()

# Go through each region separately
# Go through each region separately
for (i in seq_along(rev(bioGeo))) {
  region <- rev(bioGeo)[i]
  
  # Subset data for the current region
  regionData <- nectarColorRegion %>% filter(region == !!region)
  
  # Region tree
  # Cut the tree down to the species that we have in our dataset
  regionTree <- drop.tip(plantTree, setdiff(plantTree$tip.label,
                                            regionData$species))
  # How many species were matched
  setdiff(unique(regionData$species), regionTree$tip.label)
  
  # then we need to sort the trait data so that the species follow the order of the tree
  regionData <- regionData[match(regionTree$tip.label,
                                 regionData$species), ]
  
  # PGLS Model 1
  models[[region]] <- phylANOVA(regionTree, regionData$color, 
                                regionData$nectar10, 
                                p.adj="bonferroni",
                                nsim=100000)
  
  # Extract the post hoc results
  posthocNectarP <- models[[region]]$Pt
  
  # Create a vector to store p-values for pairwise comparisons
  pValues <- c()
  comparisons <- c()
  
  # Loop through pairwise comparisons to add to the string
  for (k in 1:nrow(posthocNectarP)) {
    for (j in 1:ncol(posthocNectarP)) {
      if (k < j) {  # To avoid duplicate comparisons and self-comparisons
        comparison <- paste(rownames(posthocNectarP)[k], "vs", colnames(posthocNectarP)[j])
        pValue <- posthocNectarP[k, j]
        
        # Append p-values and comparison names
        pValues <- c(pValues, pValue)
        comparisons <- c(comparisons, comparison)
      }
    }
  }
  
  # Reformat comparison names to contain only a single hyphen
  comparisons <- gsub(" vs ", "-", comparisons)
  
  # Assign the reformatted names to the p-values vector
  names(pValues) <- comparisons
  
  # Generate the letters to represent significant groupings
  lettersNectar <- multcompLetters(pValues, threshold=0.05)$Letters
  
  # Create a data frame with the regions and their corresponding significance letters
  lettersNectar <- data.frame(
    color = c(names(lettersNectar)),
    label = c(lettersNectar)
  )
  
  # Arrange as will appear on X-Axis
  # Convert 'region' column to a factor with specified levels
  lettersNectar$color <- factor(lettersNectar$color)
  lettersNectar <- lettersNectar %>% arrange(color)
  lettersNectar$axis <- c("One", "Two", "Three")
  
  # Create a ggplot for the current region
  plotsNectarColor[[region]] <- ggplot(nectarColorRegion, aes(x=axis, y=nectar10)) +
    
    # Transparent base layer to set plot dimensions without displaying data
    geom_blank() + 
    scale_y_continuous(limits=c(min(nectarColorRegion$nectar10)-0.1,
                                max(nectarColorRegion$nectar10)+0.1)) +
    
    # Add violin plot from `regionData`
    geom_violin(data=regionData, aes(x=axis, y=nectar10, fill=axis), 
                trim=FALSE, alpha=0.4) +
    
    # Add jitter layer from `colorRegion` in grey
    geom_jitter(data=nectarColorRegion[!nectarColorRegion$region == region, ],
                aes(x=axis, y=nectar10), 
                bg="grey", alpha=0.5, width=0.1, size=4, 
                pch=nectarColorRegion[!nectarColorRegion$region == region, ]$pch) +
    
    # Add jitter layer from `regionData` with original color and shape settings
    geom_jitter(data=regionData, aes(x=axis, y=nectar10), 
                width=0.1, size=4, pch=regionData$pch, 
                bg=regionData$col) +
    
    # Boxplot from `regionData` to show central tendency and spread
    geom_boxplot(data=regionData, aes(x=axis, y=nectar10), 
                 width=0.2, fill=adjustcolor("ghostwhite", alpha.f=0.7), outlier.shape=NA) +  
    
    coord_flip() +
    
    theme_bw() +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          panel.border=element_rect(colour="black", fill=NA, 
                                    size=1.5)) +
    
    scale_x_discrete(limits=c("One", "Two", "Three")) +
    scale_fill_manual(values=c("ghostwhite", "#ffd8b1", "indianred3"), 
                      limits=c("One", "Two", "Three"),  
                      name="") + 
    
    xlab("Color category") +
    ylab(ifelse(i == length(bioGeo), "Nectar concentration", "")) +
    theme(legend.position="none", 
          legend.margin=margin(0, 0, 0, 0), 
          axis.title=element_text(size=18), 
          axis.text.x=element_text(size=0), 
          axis.text.y=element_text(size=0),
          plot.title=element_text(size=18)) +
    
    # Add the significance letters
    geom_text(data=lettersNectar, aes(x=axis, y=max(nectarColorRegion$nectar10+0.1),
                                     label=label), size=8, vjust=1) +
    
    ggtitle(paste0(region, ": ", "F=", 
                   round(models[[region]]$F, 2), ",",
                   " p=", round(models[[region]]$Pf, 3)))
  
}

# Plot and store
png("../figures/4_results/Figure 2/Region/SupplFig. 2c) CorollaNectarViolinplot_NewColors.png", 
    height=29, width=18, units="cm", res=900)
# Arrange the plots in a grid
plot_grid(plotsNectarColor[[1]], plotsNectarColor[[2]], 
          plotsNectarColor[[3]], plotsNectarColor[[4]], 
          nrow=4, ncol=1)
dev.off()

#####

# 6.4.2 Nectar color Violinplot - Family ####

# Create a list to store models and plots
models <- list()
plots <- list()

# Count species per family and order by species number
families <- as.data.frame(nectarColorRegion %>%
                            group_by(family) %>%
                            summarise(species_count = n()) %>%
                            arrange(desc(species_count)))$family

# Go through each family separately
for (family in families) {
  
  # Subset data for the current family (and change region to All for plotting)
  familyData <- nectarColorRegion %>% filter(family == !!family)
  familyData$region <- "All"
  
  # Only test this for families with at least 20 species
  if(nrow(familyData) < 7) {
    
    # Create a ggplot for the current region
    plots[[family]] <- ggplot(nectarColorRegion, aes(x=axis, y=nectar10)) +
      
      # Transparent base layer to set plot dimensions without displaying data
      geom_blank() + 
      scale_y_continuous(limits=c(min(nectarColorRegion$nectar10)-0.1,
                                  max(nectarColorRegion$nectar10)+0.1)) +
      
      # Add violin plot from `regionData`
      geom_violin(data=familyData, aes(x=axis, y=nectar10, fill=axis), 
                  trim=FALSE, alpha=0.4) +
      
      # Add jitter layer from `colorRegion` in grey
      geom_jitter(data=nectarColorRegion[!nectarColorRegion$family == family, ], 
                  aes(x=axis, y=nectar10), 
                  bg="grey", alpha=0.5, width=0.1, size=4, 
                  pch=nectarColorRegion[!nectarColorRegion$family == family, ]$pch) +
      
      # Add jitter layer from `regionData` with original color and shape settings
      geom_jitter(data=familyData, aes(x=axis, y=nectar10), 
                  width=0.2, size=4, pch=familyData$pch, 
                  bg=familyData$col) +
      
      # Boxplot from `regionData` to show central tendency and spread
      geom_boxplot(data=familyData, aes(x=axis, y=nectar10), 
                   width=0.2, fill=adjustcolor("ghostwhite", alpha.f=0.7), 
                   outlier.shape=NA) +  
      
      coord_flip() +
      
      theme_bw() +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            panel.border=element_rect(colour="black", fill=NA, size=1.5)) +
      
      scale_x_discrete(limits=c("One", "Two", "Three")) +
      scale_fill_manual(values=c("ghostwhite", "#ffd8b1", "indianred3"), 
                        limits=c("One", "Two", "Three"),  
                        name="") + 
      
      xlab("Color category") + 
      ylab("Nectar concentration") +
      theme(legend.title = element_blank(), legend.position="none",
            axis.title = element_text(size=10), 
            axis.text.x = element_text(size=0), 
            axis.text.y = element_text(size=0),
            plot.title = element_text(size=10)) +
      ggtitle(paste0(family, ": Not enough (n<7)"))
    
  } else if (length(unique(familyData$axis)) < 2) {
    
    # Create a ggplot for the current region
    plots[[family]] <- ggplot(nectarColorRegion, aes(x=axis, y=nectar10)) +
      
      # Transparent base layer to set plot dimensions without displaying data
      geom_blank() + 
      scale_y_continuous(limits=c(min(nectarColorRegion$nectar10)-0.1,
                                  max(nectarColorRegion$nectar10)+0.1)) +
      
      # Add violin plot from `regionData`
      geom_violin(data=familyData, aes(x=axis, y=nectar10, fill=axis), 
                  trim=FALSE, alpha=0.4) +
      
      # Add jitter layer from `colorRegion` in grey
      geom_jitter(data=nectarColorRegion[!nectarColorRegion$family == family, ], 
                  aes(x=axis, y=nectar10), 
                  bg="grey", alpha=0.5, width=0.1, size=4, 
                  pch=nectarColorRegion[!nectarColorRegion$family == family, ]$pch) +
      
      # Add jitter layer from `regionData` with original color and shape settings
      geom_jitter(data=familyData, aes(x=axis, y=nectar10), 
                  width=0.1, size=4, pch=familyData$pch, 
                  bg=familyData$col) +
      
      coord_flip() +
      
      theme_bw() +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            panel.border=element_rect(colour="black", fill=NA, size=1.5)) +
      
      scale_x_discrete(limits=c("One", "Two", "Three")) +
      scale_fill_manual(values=c("ghostwhite", "#ffd8b1", "indianred3"), 
                        limits=c("One", "Two", "Three"),  
                        name="") + 
      
      xlab("Color category") + 
      ylab("Nectar concentration") +
      theme(legend.title = element_blank(), legend.position="none",
            axis.title = element_text(size=10), 
            axis.text.x = element_text(size=0), 
            axis.text.y = element_text(size=0),
            plot.title = element_text(size=10)) +
      ggtitle(paste0(family, ": Only one color"))
    
  } else {
    
    # Model 1
    models[[family]] <- glm(nectar10 ~ color, 
                            data=familyData)
    
    # Letters
    ptt2 <- aov(nectar10 ~ color, familyData)
    lettersNectar <- data.frame("Letters" = multcompLetters(extract_p(TukeyHSD(ptt2)$"color"))$"Letters")
    
    # Create a data frame with the regions and their corresponding significance letters
    lettersNectar <- data.frame(
      color = c(row.names(lettersNectar)),
      label = c(lettersNectar)
    )
    
    # Define the complete set of colors you expect
    expectedColors <- 1:3
    
    # Check for missing colors
    missingColors <- setdiff(expectedColors, lettersNectar$color)
    
    # If there are missing colors, add them with "na" in the Letters column
    if (length(missingColors) > 0) {
      # Create a dataframe with the missing colors and "na" in the Letters column
      missingDf <- data.frame(color=missingColors, Letters = rep("na", length(missingColors)))
      
      # Bind the missing colors to the original dataframe
      lettersNectar <- rbind(lettersNectar, missingDf)
    }
    
    # Arrange as will appear on X-Axis
    # Convert 'region' column to a factor with specified levels
    lettersNectar$color <- factor(lettersNectar$color)
    lettersNectar <- lettersNectar %>% arrange(color)
    lettersNectar$axis <- c("One", "Two", "Three")
    
    # Create a ggplot for the current region
    plots[[family]] <- ggplot(nectarColorRegion, aes(x=axis, y=nectar10)) +
      
      # Transparent base layer to set plot dimensions without displaying data
      geom_blank() + 
      scale_y_continuous(limits=c(min(nectarColorRegion$nectar10)-0.1,
                                  max(nectarColorRegion$nectar10)+0.1)) +
      
      # Add violin plot from `regionData`
      geom_violin(data=familyData, aes(x=axis, y=nectar10, fill=axis), 
                  trim=FALSE, alpha=0.4) +
      
      # Add jitter layer from `colorRegion` in grey
      geom_jitter(data=nectarColorRegion[!nectarColorRegion$family == family, ], 
                  aes(x=axis, y=nectar10), 
                  bg="grey", alpha=0.5, width=0.1, size=4, 
                  pch=nectarColorRegion[!nectarColorRegion$family == family, ]$pch) +
      
      # Add jitter layer from `regionData` with original color and shape settings
      geom_jitter(data=familyData, aes(x=axis, y=nectar10), 
                  width=0.1, size=4, pch=familyData$pch, 
                  bg=familyData$col) +
      
      # Boxplot from `regionData` to show central tendency and spread
      geom_boxplot(data=familyData, aes(x=axis, y=nectar10), 
                   width=0.2, fill=adjustcolor("ghostwhite", alpha.f=0.7), 
                   outlier.shape=NA) +  
      
      coord_flip() +
      
      theme_bw() +
      theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
            panel.border=element_rect(colour="black", fill=NA, size=1.5)) +
      
      scale_x_discrete(limits=c("One", "Two", "Three")) +
      scale_fill_manual(values=c("ghostwhite", "#ffd8b1", "indianred3"), 
                        limits=c("One", "Two", "Three"),  
                        name="") + 
      
      xlab("Color category") +
      ylab("Nectar concentration") +
      theme(legend.position="none", 
            legend.margin=margin(0, 0, 0, 0), 
            axis.title=element_text(size=10), 
            axis.text.x=element_text(size=0), 
            axis.text.y=element_text(size=0),
            plot.title=element_text(size=10)) +
      
      # Add the significance letters
      geom_text(data=lettersNectar, aes(x=axis, y=max(nectarColorRegion$nectar10+0.1),
                                       label=Letters), size=8, vjust=1) +
      
      ggtitle(bquote(.(family) * ": " * chi^2 * "=" * .(round(Anova(models[[family]], type = c("II","III", 2, 3), 
                                                                    vcov.= vcov(models[[family]], complete=TRUE))[[1]], 2)) * 
                       "," ~ "p=" * .(round(Anova(models[[family]], type = c("II","III", 2, 3), 
                                                  vcov.= vcov(models[[family]], complete=TRUE))[[3]], 3))))
    
  }
  
}

# Display plots (example to show the first plot)
# Plot and store
# Define the ranges of family indices for separate plots
familyRanges <- list(1:15, 16:30, 31:45, 46:55)

# Set a counter
counter <- 1

# Now let's go through all the webs ans plot them
# But we want to go from south to North (i.e. from Trinidad to Saba)
for(i in seq_along(familyRanges)) {
  
  # Define the current range of families
  current <- familyRanges[[i]]
  
  # Create the file name with the counter prepended
  fileName <- paste0("../figures/4_results/Figure 2/Family/SupplFig. 2c) NectarColorViolinplotFamily_NewColors",
                     counter, ".png")
  
  # plot it 
  png(fileName, height=29.7, width=21, units="cm", res=900)
  
  # Arrange the plots for the current family range
  print(plot_grid(plotlist=plots[current], 
                  nrow=5, ncol=3))
  
  # close and save
  dev.off()
  
  # Increment the counter for the next iteration
  counter <- counter + 1
  
}

# Now plot significant models only for supplements
# Initialize an empty vector to store names of significant models
significantModels <- c()

# Loop through each model in the list
for (modelName in names(models)) {
  model <- models[[modelName]]
  # Extract p-values from the model summary
  p <- round(Anova(model, type = c("II","III", 2, 3), 
                   vcov.= vcov(model, complete=TRUE))[[3]], 3)
  
  # Check if any predictor has p < 0.05
  if (any(p < 0.05, na.rm = TRUE)) {
    # If significant, add model name to the vector
    significantModels <- c(significantModels, modelName)
  }
}

# Calculate rows based on the number of significant models
numPlots <- length(significantModels)
ncol <- 3
nrow <- ceiling(numPlots / ncol)

# Create a list of plots for significant models
significantPlots <- lapply(significantModels, function(modelName) plots[[modelName]])

# Define the file name for saving the plot
fileName <- "../figures/4_results/Figure 2/Family/SupplFig. 2c) NectarColorViolinplotFamily_Significant_NewColors.png"

# Plot and save
png(fileName, height=21, width=29, units="cm", res=900)
print(plot_grid(plotlist=significantPlots, nrow=2, ncol=ncol))
dev.off()

# Clean working environment
rm(list = ls())

#####
