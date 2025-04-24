################################################################################
##           Floral traits across American biogeographical regions            ##
##                           Map of the study region                          ##
################################################################################

# 1. Load Required Libraries ####
required_libraries <- c("raster", "sf", "colorspace")

# Install missing libraries
missing_libs <- required_libraries[!(required_libraries %in% rownames(installed.packages()))]
if (length(missing_libs) > 0) {
  install.packages(missing_libs)
}

# Load libraries
invisible(lapply(required_libraries, library, character.only = TRUE))
rm(missing_libs, required_libraries)

#####

# 2. Set Working Directory ####

setwd("../data")

#####

# 3. Map ####
# Load topography raster
topo <- raster("3_mapdata/DTM2minAme.tif")

# Load shapefiles
andes <- shapefile("3_mapdata/Mountains.shp")

# Load network metadata
netInfo <- read.csv("2_raw/networkMetadata.csv", header = TRUE, sep = ",", 
                    row.names=1)

# Define spatial extent (latitude and longitude)
CP <- as(extent(-150, -30, -40, 50), "SpatialPolygons")
elev <- crop(topo, CP)  # Crop elevation data to extent

# Terrain Analysis
slope <- terrain(elev, opt = 'slope')
aspect <- terrain(elev, opt = 'aspect')
hill <- hillShade(slope, aspect, 40, 270)

# Prepare coordinates for plotting
coords <- cbind(netInfo$longitude, netInfo$latitude)
coord2 <- SpatialPointsDataFrame(coords, 
                                 data = data.frame(ID = netInfo$network),
                                 proj4string = CRS("+proj=longlat +datum=WGS84"))

# Define plotting symbols and colors based on biogeographical regions
pch <- rep(21, nrow(netInfo))
biogeogr <- as.character(netInfo$region)

# Customize symbols
pch[biogeogr == "Central & North America"] <- 24
pch[biogeogr == "Andes"] <- 22
pch[biogeogr == "Caribbean"] <- 23

# Customize colors
col <- rep("beige", nrow(netInfo)) # "#eb4450"
col[biogeogr == "Central & North America"] <- "#f1a35f"
col[biogeogr == "Andes"] <- "bisque4" # #f8704d
col[biogeogr == "Caribbean"] <- "cornflowerblue" # #e9d58f

# Plotting the map
quartz(11, 10)  # Use quartz on macOS; replace with windows() on Windows
plot(hill, col = "#0d182e", legend = FALSE, axes = FALSE, box = FALSE, alpha = 0.4)
plot(elev, col = rev(grey(1:100 / 100, alpha = 0.5)), add = TRUE, legend = FALSE)
plot(coord2, add = TRUE, pch = pch, bg = col)

# Specify the path to the folder
folderPath <- "../figures/1_map"

# Check if the folder exists; if not, create it
if (!dir.exists(folderPath)) {
  dir.create(folderPath)
  message("Folder created at: ", folderPath)
} else {
  message("Folder already exists at: ", folderPath)
}

# Save the plot
invisible(dev.copy2pdf(file="../figures/1_map/MapAmericaSimplified_NewColors.pdf", 
                       useDingbats=FALSE, family="sans"))
invisible(dev.off())

# Legend
quartz(4, 6)

pdf(file="../figures/1_map/legendBiogeographicalRegions_NewColors.pdf", 
    width=3.5, height=3, useDingbats=FALSE, family="sans")
# Initialize a blank plot area
plot.new()
# Add the legend to the top right
legend("topright", legend=c("Caribbean", "Central & North America",
                            "Andes", "Lowland South America"), 
       pch=c(23, 24, 22, 21),
       pt.bg=c("cornflowerblue", "#f1a35f", "bisque4", "beige"),
       col="black")
# Close the PDF device
dev.off()

#####