#################################################################################################
###
### Goal: Network Eigengene Analysis for the traits: Bc_24hpi and Ch_22hpi
###
### Made by: Cynthia Soto 
### Date: Feb 18, 2022
### Latest update: Apr 18, 2022
###
###

# BASE CONFIGURATION #############################################################################

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);
# Clear object lists
rm(list = ls());

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
library(here);    ## allows you to set the top level of your project folder as “here” and to specify where things live relative to that location   
here::here();     # Top level dir: /data/run/cyntsc/Project_athal_wgcna

# Load the expression data and trait data saved in the first part
lnames = load(file = "Athal_Infected_MatrixE_simulated5000.RData");                  
lnames = load(file = "Athal_Infected_Module_Identification_MatrixE_simulated5000.RData");   
lnames = load(file = "Athal_Infected_MergedMods_MatrixE_simulated5000.RData");
lnames = load(file = "Athal_Infected_DissTOM_MatrixE_simulated5000.RData");
lnames = load(file = "Athal_Infected_DynamicMods_MatrixE_simulated5000.RData");

nGenes = ncol(athalData3);
nSamples = nrow(athalData3);

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
# dissTOM[1:10,1:5]
plotTOM = dissTOM^9;
plotTOM[1:10,1:5]
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, 
        moduleColors, 
        terrainColors = TRUE,
        main = "Arabidopsis infected Ntw Heatmap, all genes")

#########################     Visualizing the network of Eigengene   ###################################
##
## Use the Eigengenes as representative profiles and quantify module similarity by Eigengene correlation
## 
## 

# Recalculate module Eigengenes
MEs = moduleEigengenes(athalData3, moduleColors)$eigengenes
# Isolate B24 or Ch22 from the clinical traits
##################  First B24  ##################
Trait = as.data.frame(datTraits2$B_24hpi);
names(Trait) = "Bc_24"
hm_title = 'Botrytis 24 hpi'

##################  First Ch22  ##################
Trait = as.data.frame(datTraits2$Ch_22hpi);
names(Trait) = "Ch_22"
hm_title = 'Colletotrichum 22 hpi'

#################################################

# Add the Trait to existing module eigengenes
MET = orderMEs(cbind(MEs, Trait))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)

# Generates a summary plot of the eigengene network. 
# It is usually informative to add a clinical trait (or multiple traits) to the eigengenes to see how the traits fit into the eigengene network.
# the function produces a dendrogram of the eigengenes and trait(s), and a heatmap of their relationship
plotEigengeneNetworks(MET, paste("Eigengene Ntw for ", hm_title), marDendro = c(0,4,1,2), marHeatmap = c(4,4,1,1), cex.lab = 0.8, xLabelsAngle = 90)

# Split the dendrogram and heatmap plots
# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene dendrogram for ", hm_title), marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, paste("Eigengene adjacency heatmap for ", hm_title), marHeatmap = c(9,5,2,1),
                      xLabels = names(MEs), plotDendrograms = FALSE,  xLabelsAngle = 90)
