################################################################################
##           Floral traits across American biogeographical regions            ##
##                    Phylogenetic tree of plant species.                     ##
################################################################################

# 1. Libraries ####

# Ensure necessary libraries are installed and loaded
libraries <- c(
  "V.PhyloMaker2", "ggtree", "dplyr", "RColorBrewer", "ggplot2"
  )

installedLibs <- libraries %in% rownames(installed.packages())
if (any(!installedLibs)) {
  install.packages(libraries[!installedLibs])
}

invisible(lapply(libraries, library, character.only=TRUE))
rm(installedLibs, libraries)

#####

# 2. Set Working Directory ####

setwd("../data")

#####

# 3. Data ####

# Plant traits
plants <- read.csv("2_raw/plantTraitsIntroduced.csv", header=TRUE, sep=",",
                   row.names=1)

# Now keep only what's needed
# The tree we have here matches well with the old species names and not so perfectly with
# the new species names
plantSpecies <- plants[ ,c(2, 1)] # keep only species (2) and family (1) columns
colnames(plantSpecies) <- c("species", "family")
# Create the "genus" column by extracting the first word of the "species" column
plantSpecies$genus <- sapply(strsplit(plantSpecies$species, " "), `[`, 1)
# Reorder the columns to place "genus" between "species" and "family"
plantSpecies <- plantSpecies[, c("species", "genus", "family")]
# Add two empty columns for "species.relative" and "genus.relative"
plantSpecies$species.relative <- ""
plantSpecies$genus.relative <- ""
# Inspect the updated data frame
head(plantSpecies)

# Quick check families
sort(unique(plantSpecies$family))

# get rid of duplicated species names
plantSpecies <- plantSpecies %>% distinct(species, .keep_all=TRUE)

#####

# 4. Create and Store the Tree ####

# Specify the path to the folder
folderPath <- "4_phylogeny"

# Check if the folder exists; if not, create it
if (!dir.exists(folderPath)) {
  dir.create(folderPath)
  message("Folder created at: ", folderPath)
} else {
  message("Folder already exists at: ", folderPath)
}

# Create plant tree
plantTree <- phylo.maker(sp.list=plantSpecies, tree=GBOTB.extended.TPL,
                         nodes=nodes.info.1.TPL, scenarios="S3")

# Export the tree
write.tree(plantTree$scenario.3, "4_phylogeny/PlantTree")

# And store the metadata (i.e., species list as it was used for creation of tree. can later be used to match traits etc)
metadata <- plantTree$species.list
metadata$species <- gsub(" ", "_", metadata$species)
write.csv(metadata, file="4_phylogeny/PlantMetaData.csv")

#####

# 5. Draw tree and save to PNG ####

# As the nodes are not necessarily clearly marked, we add the information
# Make dataframe for clade nodes
cladesDf <- data.frame(
  
  clade=unique(metadata$family),
  node=NA
  
)

#Find the most recent common ancestor for each clade
for (i in 1:length(cladesDf$clade)) {
  
  cladesDf$node[i] <- MRCA(
    plantTree$scenario.3,
    metadata$species[metadata$family == cladesDf$clade[i]]
  )
  
}

# Make some colors (see colors chosen for PlantPieCharts)
plantColor <- data.frame(clade=c(names(sort(table(plantSpecies$family), 
                                        decreasing=TRUE))[1:11], "Other"))

# Or get colors automatically from an existing palette
plantColor$color <- c(hcl.colors(n=11, palette="PuBuGn"), "transparent") 
# Add color to clade.df
cladesDf$color <- plantColor$color[match(cladesDf$clade, plantColor$clade)]
cladesDf$color[is.na(cladesDf$color)] <- "transparent"

# Displaying all names of plant families makes the tree very messy
# Throw out a few clades that aren't very abundant anyway
cladeSelect <- data.frame(clade=c(names(sort(table(plantSpecies$family), 
                                               decreasing=TRUE))), 
                          count=as.numeric(sort(table(plantSpecies$family), 
                                                  decreasing=TRUE)))
cladesOut <- cladeSelect[which(cladeSelect$count < 2), ]$clade

# Plot the figure nicely

# Specify the path to the folder
folderPath <- "../figures/3_phylogeneticTree"

# Check if the folder exists; if not, create it
if (!dir.exists(folderPath)) {
  dir.create(folderPath)
  message("Folder created at: ", folderPath)
} else {
  message("Folder already exists at: ", folderPath)
}

png("../figures/3_phylogeneticTree/PhylogeneticTreePlantsAmerica.png", height=20, width=20, 
    units="cm", res=900)

ggtree(plantTree$scenario.3, layout="circular", linetype=NA) %<+% metadata +
  
  geom_highlight(data=cladesDf, 
                 aes(node=node, fill=color),
                 alpha=1,
                 align="right",
                 extend=0.04,
                 show.legend=FALSE) +
  
  geom_cladelab(data=cladesDf[!cladesDf$clade %in% cladesOut, ],
                mapping=aes(node=node, label=clade),
                fontsize=2,
                align="TRUE",
                angle="auto",
                offset=0.04,
                offset.text=0.01) +
  
  geom_tree(linewidth=0.25) +
  geom_tippoint() +
  scale_fill_identity()  

# Save the plot
dev.off()

#####
