################################################################################
##           Floral traits across American biogeographical regions            ##
##                       Plant traits across regions                          ##
################################################################################

# 1. Load Required Libraries ####

# Libraries
libraries <- c(
  "dplyr", "phytools", "RColorBrewer", "tidyverse", "emmeans",
  "car", "multcompView", "ape", "multcomp", "cowplot", "greekLetters",
  "MCMCglmm"
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

# 3. Plant data and network Information ####

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

# 4. Corolla
# 4.1 Data collection corolla length ####

# Filter the main dataset to retain only species with corolla length data
corolla <- plants %>% drop_na(corolla)

# Initialize empty vectors to store collected information
lapply(c("regVec", "specVec", "morphVec", "famVec"), 
       function(x) assign(x, vector(), envir = .GlobalEnv))

# Loop through each biogeographical region in bioGeo
for(h in seq_along(bioGeo)){
  
  # Get the network identities for the current region
  ma <- networkInfo$network[networkInfo$region == bioGeo[h]]
  
  # Filter rows in "corolla" where network overlaps with current region's networks
  corollaFiltered <- corolla %>%
    rowwise() %>%
    filter(any(as.numeric(strsplit(as.character(network), ",")[[1]]) %in% ma))
  
  # Append collected data for the region
  regVec <- c(regVec, rep(bioGeo[h], nrow(corollaFiltered)))
  specVec <- c(specVec, corollaFiltered$species)
  morphVec <- c(morphVec, corollaFiltered$corolla)
  famVec <- c(famVec, corollaFiltered$family)
}

# Create the final data frame with the collected information
corollaRegion <- data.frame(family=famVec,
                            region=regVec,
                            corolla=morphVec,
                            species=specVec)

# Clean up environment
rm(famVec, h, ma, morphVec, regVec, specVec, corollaFiltered)

# Now I adapt the plant names to the format of the phylogenetic tree
corollaRegion$species <- gsub(" ", "_", corollaRegion$species)

# And finally I log10 transform the corolla length
corollaRegion$log10corolla <- log10(corollaRegion$corolla + 1)

# throw out potential duplicate entries (some species occur in multiple networks in
# the same biogeographical region)
corollaRegion <- distinct(corollaRegion)

#####

# 4.2 Match data corolla and phylogenetic tree ####

# Read the plant tree
plantTree <- read.tree("4_phylogeny/PlantTree")

# Cut the tree down to the species that we have in our dataset
# Use old species, as many more match
plantTree <- drop.tip(plantTree, setdiff(plantTree$tip.label,
                                           corollaRegion$species))
# How many species were matched
setdiff(unique(corollaRegion$species), plantTree$tip.label)

# Next step is to reorder the trait dataframe by putting the species names in the 
# same order as they appear in the phylogenetic tree. This step will also exclude 
# multiple entries of plant names (plants may occur in multiple biogeographical regions).
# If a species occurs also in the Caribbean make sure to keep it, seeing how
# ther are very few species there to begin with
multipleSpecies <- names(table(corollaRegion$species)[table(corollaRegion$species) > 1])
# Filter multipleSpecies to keep only those species that occur in the "Caribbean" region
caribbeanSpecies <- multipleSpecies[sapply(multipleSpecies, function(species) {
  "Caribbean" %in% corollaRegion[corollaRegion$species == species, ]$region
})]
# Now keep only those rows, where the biogeographical region is Caribbean
corollaRegion <- corollaRegion[!(corollaRegion$species %in% caribbeanSpecies & 
                                   corollaRegion$region != "Caribbean"), ]

# then we need to sort the trait data so that the species follow the order of the tree
corollaRegion <- corollaRegion[match(plantTree$tip.label,
                                     corollaRegion$species), ]

# Clean working environment
rm(caribbeanSpecies, multipleSpecies)

#####

# 4.3 Corolla phylANOVA ####

# First check the phylogenetic signal
phyloSignalCorolla <- phylosig(plantTree, corollaRegion$log10corolla, 
                               method="lambda", test=TRUE)
lambdaCorolla <- phyloSignalCorolla$lambda

# I use log transformed corolla length
phylModelCorolla <- phylANOVA(plantTree, corollaRegion$region, 
                              corollaRegion$log10corolla, 
                              p.adj="bonferroni",
                              nsim=100000)

# Store results
# Specify the path to the folder
folderPath <- "../results/textFiles/Figure 3"

# Check if the folder exists; if not, create it
if (!dir.exists(folderPath)) {
  dir.create(folderPath)
  message("Folder created at: ", folderPath)
} else {
  message("Folder already exists at: ", folderPath)
}

# Generate ANOVA results section
resultsCorolla <- sprintf(
  "ANOVA table: Phylogenetic ANOVA\n\nResponse: Corolla length\n%-15s %-10s %-10s\n%-15s %-10.5f %-10.5f\n%-15s %-10.5f %-10.5f\n\nF value: %.2f\nP value: %.3f\n\n",
  "Source", "Sum Sq", "Mean Sq",
  "Region", phylModelCorolla$`Sum Sq`[1], phylModelCorolla$`Mean Sq`[1],
  "Residual", phylModelCorolla$`Sum Sq`[2], phylModelCorolla$`Mean Sq`[2],
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
writeLines(outputCorolla, "../results/textFiles/Figure 3/phylANOVA(Corolla).txt")

# Test for homogeneity of variance
homogeneityCorolla <- leveneTest(log10corolla ~ region, data = corollaRegion)$`Pr(>F)`[1]

# Summary generation
summaryCorolla <- sprintf(
  "Corolla Length\nPhylogenetic signal λ = %.2f\nThe phylogenetic signal λ suggests a %s phylogenetic effect (p = %.3f).\n\nPhylogenetic ANOVA Results:\nF = %.2f\np = %.4f\n\nHomogeneity of variance: %s",
  lambdaCorolla,
  ifelse(lambdaCorolla < 0.3, "weak", ifelse(lambdaCorolla < 0.7, "moderate", "strong")),
  phyloSignalCorolla$P,
  phylModelCorolla$F,
  phylModelCorolla$Pf,
  ifelse(homogeneityCorolla < 0.05, "is not given", "is given")
)

# Now let's run a few different tests to confirm.
# log10corolla is not normal, but Homogeneity of variance is given
# so let's just run a Kruskal Wallis
kT <- kruskal.test(log10corolla ~ region, data = corollaRegion)

# Append Kruskal-Wallis test results
summaryCorolla <- paste(
  summaryCorolla,
  sprintf(
    "\n\nKruskal-Wallis Test Results:\nChi-squared = %.2f\nDegrees of Freedom = %d\np-value = %.4f",
    kT$statistic, kT$parameter, kT$p.value
  )
)

# Output summary (e.g., print, save, etc.)
writeLines(summaryCorolla, "../results/textFiles/Figure 3/summaryCorolla.txt")

# Clean up
rm(phyloSignalCorolla, homogeneityCorolla, lambdaCorolla,
   outputCorolla, posthocCorolla, resultsCorolla, summaryCorolla, 
   comparisons, pValues, kT, folderPath, comparisonNames)

#####

# 4.4 Corolla Violinplot ####

# First copy this dataframe and append it below itself but change region to "All"
corollaRegionPlot <- rbind(corollaRegion, transform(corollaRegion, region = "All"))
# Add some plotting parameters
corollaRegionPlot$pch <- ifelse(corollaRegionPlot$region == "Lowland South America", 21,
                                ifelse(corollaRegionPlot$region == "Andes", 22,
                                       ifelse(corollaRegionPlot$region == "Caribbean", 23, 
                                              ifelse(corollaRegionPlot$region == "All", NA, 24))))

# Create a data frame with the regions and their corresponding significance letters
lettersCorolla <- data.frame(
  region = c(names(lettersCorolla), "All"),
  label = c(lettersCorolla, NA)
)

# Arrange as will appear on X-Axis
# Convert 'region' to a factor with specified levels and arrange
lettersCorolla <- lettersCorolla %>%
  mutate(region = factor(region, levels = c(bioGeo, "All"))) %>%
  arrange(region)

# Text for X axis Labels
bioGeoModified <- c("Lowland\nSouth America", "Andes", "Central & \nNorth America", "Caribbean", "All")

# Store figure
figureManuscript <- list()

figureManuscript[[1]] <- ggplot(corollaRegionPlot, aes(x=region, y=log10corolla, fill=region)) +
  
  geom_violin(trim=FALSE, alpha=0.8) + 
  geom_jitter(width=0.1, size=2, pch=corollaRegionPlot$pch, 
              bg="black") +
  geom_boxplot(width=0.25, fill=adjustcolor("ghostwhite", alpha.f = 0.8), 
               outlier.shape=NA) + 
  
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.border=element_rect(colour="black", fill=NA, linewidth=1.5)) +
  
  scale_x_discrete(limits=c(bioGeo, "All"), labels=bioGeoModified) +  
  
  scale_fill_manual(values=c(
                            alpha("beige", 0.4), # #eb4450
                            alpha("bisque4", 0.4), # #f8704d
                            alpha("#f1a35f", 0.4),
                            alpha("cornflowerblue", 0.4), # #e9d58f
                            alpha("grey80", 0.4)
                            ),
                    limits=c(bioGeo, "All"),
                    name="") +
  
  scale_y_continuous(
    breaks = seq(0, max(corollaRegionPlot$log10corolla), by = 0.5), # Adjust as needed
    labels = scales::label_number(accuracy = 0.01) # Ensure 2 decimal places
  ) +
  
  labs(x="", y="Log corolla length (mm)") +
  theme(legend.position="none", 
        legend.margin=margin(0, 0, 0, 0), 
        axis.title=element_text(size=30), 
        axis.text.x=element_text(size=0), 
        axis.text.y=element_text(size=20, angle=90, hjust=0.5),
        plot.title=element_text(size=0)) +
  
  # add a vertical line between the category "all" and the others
  geom_vline(xintercept=which(bioGeo == bioGeo[4]) + 0.5, 
             linetype="dotted", 
             color="black", linewidth=1) +
  # Add the significance letters
  geom_text(data=lettersCorolla, aes(x=region, y=max(corollaRegionPlot$log10corolla+0.5),
                                   label=label), 
            size=8, vjust=1) 

# Specify the path to the folder
folderPath <- "../figures/4_results/Figure 3"

# Check if the folder exists; if not, create it
if (!dir.exists(folderPath)) {
  dir.create(folderPath)
  message("Folder created at: ", folderPath)
} else {
  message("Folder already exists at: ", folderPath)
}

# Plot and store
png("../figures/4_results/Figure 3/Fig. 3a) PhyloAnovaCorollaViolin_NewColors_Ylab.png", 
    height=21, width=29.7, units="cm", res=900)
print(figureManuscript[[1]])
dev.off()

# Clean working environment
rm(corolla, corollaRegion, folderPath, plantTree,
   phylModelCorolla)

#####

# 5. Nectar
# 5.1 Data collection nectar concentration ####

# First create a sub dataframe, with only the species where we have info on nectar concentration
# Cut the data down to the relevant plants
nectar <- plants %>%
  drop_na(nectar)

# Lose all entries with nectar < 5% as they are most likely wrong measurements
nectar <- nectar[!nectar$nectar < 5, ]

# Initialize empty vectors to store collected information
lapply(c("regVec", "specVec", "morphVec", "famVec"), 
       function(x) assign(x, vector(), envir = .GlobalEnv))

# Loop through each biogeographical region in bioGeo
for(h in seq_along(bioGeo)){
  
  # Get the network identities for the current region
  ma <- networkInfo$network[networkInfo$region == bioGeo[h]]
  
  # Filter rows in "corolla" where network overlaps with current region's networks
  nectarFiltered <- nectar %>%
    rowwise() %>%
    filter(any(as.numeric(strsplit(as.character(network), ",")[[1]]) %in% ma))
  
  # Append collected data for the region
  regVec <- c(regVec, rep(bioGeo[h], nrow(nectarFiltered)))
  specVec <- c(specVec, nectarFiltered$species)
  morphVec <- c(morphVec, nectarFiltered$nectar)
  famVec <- c(famVec, nectarFiltered$family)
}

# Create the final data frame with the collected information
nectarRegion <- data.frame(family=famVec,
                            region=regVec,
                            nectar=morphVec,
                            species=specVec)

# Clean up environment
rm(famVec, h, ma, morphVec, regVec, specVec, nectarFiltered)

# Now I adapt the plant names to the format of the phylogenetic tree
nectarRegion$species <- gsub(" ", "_", nectarRegion$species)

# throw out any duplicate entries (some species occur in multiple networks in
# the same biogeographical region)
nectarRegion <- distinct(nectarRegion)

#####

# 5.2 Match data nectar and phylogenetic tree ####

# Read the plant tree
plantTree <- read.tree("4_phylogeny/PlantTree")

# Cut the tree down to the species that we have in our dataset
plantTree <- drop.tip(plantTree, setdiff(plantTree$tip.label,
                                         nectarRegion$species))
# How many species were matched
setdiff(unique(nectarRegion$species), plantTree$tip.label)

# Next step is to reorder the trait dataframe by putting the species names in the 
# same order as they appear in the phylogenetic tree. This step will also exclude 
# multiple entries of plant names.
# As some plants occur in multiple biogeographical regions, I will try to figure out
# which ones are also in the Caribbean (if) and make sure to keep those, seeing how
# we have very few species there to begin with
multipleSpecies <- names(table(nectarRegion$species)[table(nectarRegion$species) > 1])
# Filter multipleSpecies to keep only those species that occur in the "Caribbean" region
caribbeanSpecies <- multipleSpecies[sapply(multipleSpecies, function(species) {
  "Caribbean" %in% nectarRegion[nectarRegion$species == species, ]$region
})]
# Now keep only those rows, where the biogeographical region is Caribbean
nectarRegion <- nectarRegion[!(nectarRegion$species %in% caribbeanSpecies & 
                                 nectarRegion$region != "Caribbean"), ]

# then we need to sort the trait data so that the species follow the order of the tree
nectarRegion <- nectarRegion[match(plantTree$tip.label,
                                   nectarRegion$species), ]

# Clean working environment
rm(caribbeanSpecies, multipleSpecies)

#####

# 5.3 Nectar phylANOVA ####

# First check the phylogenetic signal
phyloSignalNectar <- phylosig(plantTree, nectarRegion$nectar, method="lambda",
                              test=TRUE)
lambdaNectar <- phyloSignalNectar$lambda

# Phylogenetic Anova
phylModelNectar <- phylANOVA(plantTree, nectarRegion$region, 
                             nectarRegion$nectar,
                             p.adj="bonferroni", 
                             nsim=100000)

# Store it:
# Generate ANOVA results section
resultsNectar <- sprintf(
  "ANOVA table: Phylogenetic ANOVA\n\nResponse: Nectar\n%-15s %-10s %-10s\n%-15s %-10.5f %-10.5f\n%-15s %-10.5f %-10.5f\n\nF value: %.2f\nP value: %.3f\n\n",
  "Source", "Sum Sq", "Mean Sq",
  "Region", phylModelNectar$`Sum Sq`[1], phylModelNectar$`Mean Sq`[1],
  "Residual", phylModelNectar$`Sum Sq`[2], phylModelNectar$`Mean Sq`[2],
  phylModelNectar$F, phylModelNectar$Pf
)

# Format post hoc results with pairwise comparisons and p-values
posthocNectar <- "Pairwise posthoc test using method = 'bonferroni'\n\nComparison                P-value\n"
comparisons <- combn(colnames(phylModelNectar$Pt), 2)
pValues <- mapply(function(i, j) phylModelNectar$Pt[i, j], comparisons[1,], comparisons[2,])
comparisonNames <- apply(comparisons, 2, function(pair) paste(pair, collapse = "-"))
names(pValues) <- comparisonNames
posthocNectar <- paste(posthocNectar, paste(sprintf("%-24s %-10.3f", comparisonNames, pValues), collapse = "\n"))

# Generate letters for significant groupings
lettersNectar <- multcompLetters(pValues, threshold = 0.05)$Letters
posthocNectar <- paste(posthocNectar, "\n\nSignificance Letters:\n\n", paste(sprintf("%-24s %s", names(lettersNectar), lettersNectar), collapse = "\n"))

# Combine results
outputNectar <- paste(resultsNectar, posthocNectar)

# Write output to file
writeLines(outputNectar, "../results/textFiles/Figure 3/phylANOVA(Nectar).txt")

# Test for homogeneity of variance
homogeneityNectar <- leveneTest(nectar ~ region, data=nectarRegion)$`Pr(>F)`[1]

# Summary generation
summaryNectar <- sprintf(
  "Nectar Concentration\nPhylogenetic signal λ = %.2f\nThe phylogenetic signal λ suggests a %s phylogenetic effect (p = %.3f).\n\nPhylogenetic ANOVA Results:\nF = %.2f\np = %.4f\n\nHomogeneity of variance: %s",
  lambdaNectar,
  ifelse(lambdaNectar < 0.3, "weak", ifelse(lambdaNectar < 0.7, "moderate", "strong")),
  phyloSignalNectar$P,
  phylModelNectar$F,
  phylModelNectar$Pf,
  ifelse(homogeneityNectar < 0.05, "is not given", "is given")
)

# Now let's run a few different tests to confirm.
# log10nectar is not normal, and Homogeneity of variance is not given
# so let's start with a Kruskal Wallis
kT <- kruskal.test(nectar ~ region, data=nectarRegion)
# then Welch's Anova
wA <- oneway.test(nectar ~ region, data=nectarRegion, var.equal=FALSE)

# Append Kruskal-Wallis test results
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

# Savce
writeLines(summaryNectar, "../results/textFiles/Figure 3/summaryNectar.txt")

# Clean up
rm(summaryNectar, resultsNectar, pValues, posthocNectar, outputNectar, lambdaNectar,
   homogeneityNectar, comparisonNames, wA, plantTree, phyloSignalNectar, phylModelNectar,
   kT, comparisons)

#####

# 5.4 Nectar Violinplot ####

# First I copy this dataframe and append it below itself but change region to "All"
nectarRegionPlot <- rbind(nectarRegion, transform(nectarRegion, region = "All"))
# Add some plotting parameters
nectarRegionPlot$pch <- ifelse(nectarRegionPlot$region == "Lowland South America", 21,
                                ifelse(nectarRegionPlot$region == "Andes", 22,
                                       ifelse(nectarRegionPlot$region == "Caribbean", 23, 
                                              ifelse(nectarRegionPlot$region == "All", NA, 24))))

# Create a data frame with the regions and their corresponding significance letters
lettersNectar <- data.frame(
  region = c(names(lettersNectar), "All"),
  label = c(lettersNectar, NA)
)

# Arrange as will appear on X-Axis
# Convert 'region' to a factor with specified levels and arrange
lettersNectar <- lettersNectar %>%
  mutate(region = factor(region, levels = c(bioGeo, "All"))) %>%
  arrange(region)

# plot the figure
figureManuscript[[2]] <- ggplot(nectarRegionPlot, aes(x=region, y=nectar, fill=region)) +
  
  geom_violin(trim=FALSE, alpha=0.8) +
  geom_jitter(width=0.1, size=2, pch=nectarRegionPlot$pch, 
              bg="black") +  
  geom_boxplot(width=0.25, fill=adjustcolor("ghostwhite", alpha.f=0.8), 
               outlier.shape=NA) + 
  
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.border=element_rect(colour="black", fill=NA, size=1.5)) +
  
  scale_x_discrete(limits=c(bioGeo, "All"), labels=bioGeoModified) +
  
  scale_fill_manual(values=c(
                            alpha("beige", 0.4), # #eb4450
                            alpha("bisque4", 0.4), # #f8704d
                            alpha("#f1a35f", 0.4),
                            alpha("cornflowerblue", 0.4), # #e9d58f
                            alpha("grey80", 0.4)
                            ),
                    limits=c(bioGeo, "All"),  
                    name="") +
  
  scale_y_continuous(
    breaks = seq(0, max(nectarRegionPlot$nectar + 20), by = 20), # Adjust as needed
    labels = scales::label_number(accuracy = 0.01) # Ensure 2 decimal places
  ) +
  
  labs(x="", y="Nectar concentration (%)") +
  theme(legend.position="none", 
        legend.margin=margin(0, 0, 0, 0), 
        axis.title=element_text(size=30), 
        axis.text.x=element_text(size=0), 
        axis.text.y=element_text(size=20, angle=90, hjust=0.5),
        plot.title=element_text(size=0)) +
  
  # add a vertical line between the category "all" and the others
  geom_vline(xintercept=which(bioGeo == bioGeo[4]) + 0.5, 
             linetype="dotted", 
             color="black", linewidth=1) +
  
  # Add the significance letters
  geom_text(data=lettersNectar, aes(x=region, y=max(nectarRegionPlot$nectar+20),
                                     label=label), 
            size=8, vjust=1) 

# Plot and store
png("../figures/4_results/Figure 3/Fig. 3b) PhyloAnovaNectarViolin_NewColors_Ylab.png", height=21, width=29.7, 
    units="cm", res=900)
print(figureManuscript[[2]])
dev.off()

# Clean working environment
rm(nectar, nectarRegion)

#####

# 6. Color
# 6.1 Data collection color ####

# First create a sub dataframe, with only the species where we have info on color
# Cut the data down to the relevant plants
color <- plants %>%
  drop_na(colorCode)

# Now let's make sure that there is no "ghost" data
color <- color[-which(color$colorCode == ""), ]

# Initialize empty vectors to store collected information
lapply(c("regVec", "specVec", "morphVec", "famVec"), 
       function(x) assign(x, vector(), envir = .GlobalEnv))

# Loop through each biogeographical region in bioGeo
for(h in seq_along(bioGeo)){
  
  # Get the network identities for the current region
  ma <- networkInfo$network[networkInfo$region == bioGeo[h]]
  
  # Filter rows in "corolla" where network overlaps with current region's networks
  colorFiltered <- color %>%
    rowwise() %>%
    filter(any(as.numeric(strsplit(as.character(network), ",")[[1]]) %in% ma))
  
  # Append collected data for the region
  regVec <- c(regVec, rep(bioGeo[h], nrow(colorFiltered)))
  specVec <- c(specVec, colorFiltered$species)
  morphVec <- c(morphVec, colorFiltered$colorCode)
  famVec <- c(famVec, colorFiltered$family)
}

# Create the final data frame with the collected information
colorRegion <- data.frame(family=famVec,
                            region=regVec,
                            color=morphVec,
                            species=specVec)

# Clean up environment
rm(famVec, h, ma, morphVec, regVec, specVec, colorFiltered)

# Now I adapt the plant names to the format of the phylogenetic tree
colorRegion$species <- gsub(" ", "_", colorRegion$species)

# Separate rows based on color values and repeat other columns accordingly
colorRegionClean <- colorRegion %>%
  # Split the 'color' column by commas and separate into rows
  separate_rows(color, sep = ",") %>%
  # Convert color to numeric and factor
  mutate(color = as.factor(as.numeric(color)))

# Remove duplicates (a species can have several entries of color, as it may be reported
# in different networks)
colorRegion <- colorRegionClean %>%
  distinct(family, region, color, species, .keep_all = TRUE)

# Create a binary variable of white to test if there are more white plants in one of the regions
colorRegion$colorWhite <- ifelse(colorRegion$color == 1, 1, 0)

# clean up
rm(colorRegionClean)

#####

# 6.2 Match data color and phylogenetic tree ####

# Read the plant tree
plantTree <- read.tree("4_phylogeny/PlantTree")

# Cut the tree down to the species that we have in our dataset
plantTree <- drop.tip(plantTree, setdiff(plantTree$tip.label,
                                         colorRegion$species))
# How many species were matched
setdiff(unique(colorRegion$species), plantTree$tip.label)

# If occuring in multiple regions, keep only Caribeean
multipleSpecies <- names(table(colorRegion$species)[table(colorRegion$species) > 1])
# Filter multipleSpecies to keep only those species that occur in the "Caribbean" region
caribbeanSpecies <- multipleSpecies[sapply(multipleSpecies, function(species) {
  "Caribbean" %in% colorRegion[colorRegion$species == species, ]$region
})]
# Now keep only those rows, where the biogeographical region is Caribbean
colorRegion <- colorRegion[!(colorRegion$species %in% caribbeanSpecies & 
                               colorRegion$region != "Caribbean"), ]

# then we need to sort the trait data so that the species follow the order of the tree
colorRegion <- colorRegion[match(plantTree$tip.label,
                                 colorRegion$species), ]

# clean up objects that are not needed
rm(caribbeanSpecies, multipleSpecies)

#####

# 6.3 Color model MCMCglmm #####

# Convert categorical color to numeric
numericColors <- as.numeric(factor(colorRegion$color, levels = c("1", "2", "3")))

# Get lambda
phyloSignalColor <- phylosig(plantTree, numericColors, method="lambda", test=TRUE)
lambdaColor <- phyloSignalColor$lambda

# subset data
colorRegionModel <- as.data.frame(colorRegion[, c("species", "region", 
                                                  "colorWhite", "color")])

# Remove node labels by setting them to NULL
plantTree$node.label <- NULL

# Convert phylogenetic tree into a covariance matrix 
Ainv <- inverseA(plantTree)$Ainv 

# Run the GLMM with phylogenetic correction 
# First for white
phylModelColorWhite <- MCMCglmm(
  colorWhite ~ region,  # Fixed effects (region as the predictor)
  family = "categorical",  # For binary outcome
  ginverse = list(species = Ainv),  # Phylogenetic correction
  data = colorRegionModel,  # Your data frame
  prior = list(R = list(V = 1, nu = 0.002)),  # Prior settings
  nitt = 13000, burnin = 3000, thin = 10  # MCMC settings
)

# Specify the path to the folder
folderPath <- "../results/Figure 3"

# Check if the folder exists; if not, create it
if (!dir.exists(folderPath)) {
  dir.create(folderPath)
  message("Folder created at: ", folderPath)
} else {
  message("Folder already exists at: ", folderPath)
}

# Extract the MCMC samples
mcmcSamplesWhite <- as.mcmc(phylModelColorWhite$Sol)  # Fixed effects

# Plot and store White model diagnostics
png("../results/Figure 3/MCMCglmmDiagnosticWhite.png", height=15, width=15,
    units = "cm", res = 900)
plot(mcmcSamplesWhite)  # Trace plots for fixed effects
dev.off()

# Extract and store information
emmeansWhite <- emmeans(phylModelColorWhite, ~ region, data=colorRegionModel)

# Create a data frame from the pairs results
pairsResultsWhite <- as.data.frame(pairs(emmeansWhite))

# Add a significance level column
pairsResultsWhite$significance <- ifelse(pairsResultsWhite$lower.HPD > 0, 0.04, 
                                         ifelse(pairsResultsWhite$upper.HPD < 0, 
                                                0.04, 0.1))

# Get a named vector with just the comparisons
# Create a named vector from the pairsResults_with_p data frame
vectorResultsWhite <- setNames(pairsResultsWhite$significance, 
                               gsub("\\s*-\\s*", "-",
                                    pairsResultsWhite$contrast))

# Create Letters to highlight which groups differ
lettersColorWhite <- multcompLetters(vectorResultsWhite, threshold=0.05)$Letters

# Store model output
# Capture the summary output of the model
outputColorWhite <- capture.output(summary(phylModelColorWhite))

# Capture the comparisons of the regions
comparisonColorWhite <- capture.output(pairsResultsWhite)

# Create the summary with structured content
summaryColor <- c(
  "Floral color code",
  "",  # Add an empty line
  paste("Phylogenetic signal:", greeks("lambda"), "=", 
        round(lambdaColor, 2)),  # Combined into one line
  "\nColor category White (non-ornithophilous colors)",
  "\nModel Output:",
  outputColorWhite,
  "\nComparing regions:",
  "",
  comparisonColorWhite,
  "",
  ""
)

# Save
writeLines(summaryColor, 
           "../results/textFiles/Figure 3/summaryColor_MCMCglmm.txt")

# clean up
rm(emmeansWhite, pairsResultsWhite, 
   mcmcSamplesWhite, numericColors, vectorResultsWhite,outputColorWhite,
   comparisonColorWhite, lettersColorWhite, summaryColor, 
   phylModelColorWhite, colorRegionModel, Ainv)

#####

# 6.3.1 Color model GLM #####

# As the phylogenetic signal is moderate and implementing a binomial
# phylogenetic model is highly complex and it does not run well on this data
# a GLM is justified

# Fisher's exact test
# Frequency table of color by region
colorTable <- with(colorRegion[order(factor(colorRegion$region, 
                                          levels=unique(colorRegion$region)),
                                          colorRegion$colorWhite), ], 
                   table(colorWhite, region))

# Now we can run the Fisher's exact test for everything at once
fisherRegions <- fisher.test(colorTable, simulate.p.value=TRUE, 
                             B=10000)

# Save output
writeLines(capture.output(print(fisherRegions)), 
           con = file("../results/textFiles/Figure 3/fisherTest(Color).txt"))

# GLMM
colorRegion$region <- factor(colorRegion$region, levels=bioGeo)
colorModelWhite <- glm(colorWhite ~ region, colorRegion,
                         family=binomial)

# Save some diagnostic plots
png("../results/Figure 3/diagnosticsColorRegion.png", height=14, width=14, 
    units="cm", res=900)
performance::check_model(colorModelWhite)
dev.off()

# Store output
# Capture the summary output of the model
outputColorWhite <- capture.output(summary(colorModelWhite))

# Assess whether region has significant effect
coefficientsRegion <- Anova(colorModelWhite, type="III")

# Check overdispersion
overdispersion <- sum(residuals(colorModelWhite, type="pearson")^2) / df.residual(colorModelWhite)

# Pairwise Fisher's exact test
# Frequency table of color by region
colorTableWhite <- with(colorRegion[order(factor(colorRegion$region, 
                                            levels=unique(colorRegion$region)),
                                     colorRegion$colorWhite), ], 
                   table(colorWhite, region))

# Define the regions and initialize storage for p-values
regions <- colnames(colorTableWhite)
regionPairs <- combn(regions, 2, simplify = FALSE)
pValues <- numeric(length(regionPairs))
comparisons <- character(length(regionPairs))

# Perform Fisher's Exact Test for each pair
for (i in seq_along(regionPairs)) {
  pair <- regionPairs[[i]]
  
  # Extract counts for the two regions
  counts_1 <- colorTableWhite[, pair[1]]
  counts_2 <- colorTableWhite[, pair[2]]
  
  # Perform Fisher's test on a 2x2 matrix with counts for the pair of regions
  fisherTest <- fisher.test(rbind(counts_1, counts_2))
  
  # Store the p-value and the comparison label
  pValues[i] <- fisherTest$p.value
  comparisons[i] <- paste(pair, collapse = " vs ")
}

# Apply Bonferroni correction for multiple testing
pAdj <- p.adjust(pValues, method = "bonferroni")

# Create a data frame to display the results
regionCompared <- data.frame(
  comparison = comparisons,
  pValue = pValues,
  pAdj = pAdj
)

# Get letters for different groups
lettersColorWhite <- (emmeans(colorModelWhite, pairwise ~ region, 
          adjust="bonferroni") %>%
              cld(Letters = letters))
lettersColorWhite <- tolower(setNames(gsub("\\ ", "", lettersColorWhite$.group),
                              lettersColorWhite$region))

# Let's also see with Fisher's Exact tests
# Extract the adjusted p-values and assign names
pValuesAdj <- setNames(regionCompared$pAdj, regionCompared$comparison)

# Replace " vs " with "-" in the comparison names
names(pValuesAdj) <- gsub(" vs ", "-", names(pValuesAdj))

# Now apply multcompLetters
lettersColorFisher <- multcompLetters(pValuesAdj, threshold = 0.05)$Letters

# Both approaches tell me the same thing.
# I use the Fisher vector for the figure

# Create a slimmer and more concise summary
summaryColor <- c(
  "Floral Color Analysis\n",
  
  # Phylogenetic signal
  paste("Phylogenetic signal:", greeks("lambda"), "=", round(lambdaColor, 2)),
  paste("Phylogenetic effect:", ifelse(lambdaColor < 0.3, "Weak", ifelse(lambdaColor < 0.7, "Moderate", "Strong")), "p = ", 
        round(phyloSignalColor$P, 3)),
  
  "\n\nRegion Effect (Color White):",
  capture.output(coefficientsRegion),
  
  "\n\nModel Output (Color White):",
  outputColorWhite,
  
  "\n\nRegion Effect:",
  ifelse(coefficientsRegion$`Pr(>Chisq)` < 0.05, "Region significantly affects floral color.", "No significant effect of region on floral color."),
  
  "\n\nModel Fit (Overdispersion):",
  paste("Model is", 
        ifelse(overdispersion < 0.8, "underdispersed", 
               ifelse(overdispersion <= 1.2, "well-fitted", 
                      ifelse(overdispersion <= 1.5, "slightly overdispersed", 
                             ifelse(overdispersion <= 2, "moderately overdispersed", "strongly overdispersed"))))),
  
  "\n\nPost-Hoc Region Comparisons:",
  capture.output(regionCompared),
  
  "\n\nRegional Differences (Groups):",
  "\n\nEmmeans",
  capture.output(lettersColorWhite),
  "\n\nFisher Exact",
  capture.output(lettersColorFisher)
)

# Save
writeLines(summaryColor, "../results/textFiles/Figure 3/summaryColor.txt")

# clean up
rm(phyloSignalColor, lambdaColor, colorTable, fisherRegions, colorModelWhite, 
   outputColorWhite, coefficientsRegion, overdispersion, colorTableWhite, regionPairs, 
   pValues, comparisons, pair, counts_1, counts_2, fisherTest, pAdj, regionCompared, 
   summaryColor, pValuesAdj, plantTree, i, regions, lettersColorWhite)

#####

# 6.4 Color stacked Barplot ####

# First I copy this dataframe and append it below itself but change region to "All"
tablePct <- as.data.frame(rbind(colorRegion, 
                                 transform(colorRegion, region="All")) %>%
                             count(color, region) %>%       
                             group_by(region) %>%
                             mutate(pct=prop.table(n) * 100)) 
tablePct$colorCode <- ifelse(tablePct$color == 1, "One",
                               ifelse(tablePct$color == 2, "Two", "Three"))
# Match letters to the regions and create a new column in table_pct
tablePct$lettersWhite <- lettersColorFisher[tablePct$region]

# then create a stacked barplot
figureManuscript[[3]] <- tablePct %>%
  
  ggplot() + 
  aes(region, pct, fill = rev(color)) +
  
  geom_bar(colour="black", stat="identity", alpha=0.4, size=0.8) +
  geom_text(aes(label=paste0(sprintf(colorCode))), size=6,
            position=position_stack(vjust=0.5)) +
  
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border=element_rect(colour="black", fill=NA, size=1.5)) +
  labs(x="", y="Floral color category (%)") +
  
  # arrange by region and highlight color categories
  scale_x_discrete(limits = c(bioGeo, "All"), labels=bioGeoModified) +
  scale_fill_manual(values = c("indianred3", "#ffd8b1", "ghostwhite"), 
                    name="", labels = c("One", "Two", "Three")) +
  
  scale_y_continuous(
    breaks = seq(0, 100, by = 25), # Adjust as needed
    labels = scales::label_number(accuracy = 0.01) # Ensure 2 decimal places
  ) +
  
  # Adjust the theme to remove legend and modify axis title and text sizes
  theme(legend.position="", legend.margin=margin(0, 0, 0, 0), 
        legend.key = element_rect(fill="white", colour="black"),
        axis.title = element_text(size=30), 
        axis.text.x = element_text(size=0), 
        axis.text.y = element_text(size=20, angle=90, hjust=0.5),
        plot.title=element_text(size=0)) +
  
  # Add a vertical line between the category "all" and the others
  geom_vline(xintercept = which(bioGeo == bioGeo[4]) + 0.5,
             linetype="dotted", 
             color="black", linewidth=1) +
  
  # Add the letters to highlight differences between groups
  geom_text(label=tablePct$lettersWhite, y=103, cex=8) 

# Plot and store white
png("../figures/4_results/Figure 3/Fig. 3c) GLMColorBarWhite_Ylab.png", 
    height=21, width=29.7,
    units="cm", res=900)
print(figureManuscript[[3]])
dev.off()

#####

# 7 Figure Legend ####

png("../figures/4_results/Figure 3/Legend.png",
    height=21, width=29.7, units="cm", res=900)

ggplot(nectarRegionPlot, aes(x=region, y=nectar, fill=region)) +
  
  geom_violin(trim=FALSE, alpha=0.8) +
  
  theme_bw() +
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        panel.border=element_rect(colour="black", fill=NA, size=1.5)) +
  
  scale_x_discrete(limits=c(bioGeo, "All"), labels=bioGeoModified) +
  
  labs(x="", y="Nectar concentration") +
  theme(legend.position="none", 
        legend.margin=margin(0, 0, 0, 0), 
        axis.title=element_text(size=30), 
        axis.text.x=element_text(size=25), 
        axis.text.y=element_text(size=0)) 

dev.off()

#####

# Clean
rm(list = ls())
