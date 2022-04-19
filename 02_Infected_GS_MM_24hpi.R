#################################################################################################
###
### Goal: Analysis of Gene Significance (GS) and Module Membership (MM) for Bc24hpi Y Ch22hpi
###
### Made by: Cynthia Soto 
### Date: Feb 18, 2022
### Latest update: Apr 18, 2022
###
##################################   BASE CONFIGURATION   ########################################

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

# Load the expression data and trait data saved in the first script
lnames = load(file = "Athal_Infected_MatrixE_simulated5000.RData");                  
lnames = load(file = "Athal_Infected_Module_Identification_MatrixE_simulated5000.RData");   
lnames = load(file = "Athal_Infected_MergedMods_MatrixE_simulated5000.RData");

nGenes = ncol(athalData3);
nSamples = nrow(athalData3);

##############     Quantifying module–trait associations     ##############     

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(athalData3, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits2, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

## suitable representation in a heatmap
sizeGrWindow(9, 12)
par(mar= c(3, 10, 2, 1))    # margen inferior - y(izq) + margen superior - y(der)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), " (",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot

#library(tidyverse)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits2),
               yLabels = names(MEs),
               #yLabels =str_remove(names(MEs),'ME'),     # substring(names(MEs), 3)
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),   # greenWhiteRed
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait by fungi and hpi"))

###################################################################################################
##                     Intra-modular analysis: identifying genes with high GS and MM
##
##                 Relating modules to external information and identifying important genes   
##                            Gene Significance (GS) and Module Membership (MM)
##                      ANALYSIS FOCUSED IN THE SAMPLES for BOTRYTIS CINEREA A LAS 24 HPI
##
###################################################################################################

# Define a variable Bc_24 that contains the Bc_24 reference of datTrait2 metadata
Bc_24 = as.data.frame(datTraits2$B_24hpi);
names(Bc_24) = "Bc_24hpi";
# names (colors) of the modules
modNames = substring(names(MEs), 3);
geneModuleMembership = as.data.frame(cor(athalData3, MEs, use = "p"));                #pearson corr (continuous data, linear correlation, no outliers supported)
                                                                                      #spearman rank corr (nearly all the same assumptions as p, but it doesn’t rely on normality, and your data can be ordinal as well)
                                                                                      #spearman corr (non-parametric as spearman, use with ordinal or continuous data)
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));#Calculates Student asymptotic p-value for given correlations
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(athalData3, Bc_24, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Bc_24), sep="");
names(GSPvalue) = paste("p.GS.", names(Bc_24), sep="");

##############  Intra-modular analysis: identifying genes with high GS and MM  ##############

sizeGrWindow(7, 7);
par(mfrow = c(2, 2))
#par(mfrow = c(1, 1))

####    High correlated
module = "magenta4"
module = "turquoise"
module = "brown4"
module = "darkviolet"

column = match(module, modNames);
moduleGenes = moduleColors == module;

# panel superior (1,2)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("MM in", module, "module"),
                   ylab = "GS for body Bc_24",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
# panel inferior (3,4)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("MM in", module, "module"),
                   ylab = "GS for body Ch_22",
                   cex.main = 1, cex.lab = 1, cex.axis = 1, col = module,
                   abline.color = 1)

###################################################################################################
##                     Intra-modular analysis: identifying genes with high GS and MM
##
##                 Relating modules to external information and identifying important genes   
##                            Gene Significance (GS) and Module Membership (MM)
##                      ANALYSIS FOCUSED IN THE SAMPLES for COLLETROTRICHUM A LAS 22 HPI
##
###################################################################################################

# Define variable Ch_22 containing the Ch22 column of datTrait2
Ch_22 = as.data.frame(datTraits2$Ch_22hpi);
names(Ch_22) = "Ch_22hpi"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(athalData3, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(athalData3, Ch_22, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Ch_22), sep="");
names(GSPvalue) = paste("p.GS.", names(Ch_22), sep="");

##############  Intramodular analysis: identifying genes with high GS and MM  ##############

# If required a 2x2 frame can be design with the main modules 
sizeGrWindow(7, 7);
par(mfrow = c(2, 2))

####    High correlated
module = "darkred"
module = "coral2"
module = "mediumpurple4"
module = "grey"

column = match(module, modNames);
moduleGenes = moduleColors == module;

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("MM in", module, "module"),
                   ylab = "GS for body Ch_22",
                   main = paste("Module Membership & Gene Significance\n"),
                   cex.main = 1, cex.lab = 1, cex.axis = 1, col = module,
                   abline.color = 1)

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("MM in", module, "module"),
                   ylab = "GS for body Ch_22",
                   cex.main = 1, cex.lab = 1, cex.axis = 1, col = module,
                   abline.color = 1)

########################  Summary output of network analysis results  #########################  
##
## merge this statistical information with gene annotation and write out a
## file that summarizes the most important results and can be inspected in standard spreadsheet
##
###############################################################################################  

# Genes to annotate included in the module
names(athalData3)                                  # these are all genes included in the analysis
names(athalData3)[moduleColors=="magenta4"]        # these are all genes included in the magenta4 module

# Arabidopsis annotation file avaulable in Bioconductor 
# http://bioconductor.org/packages/2.4/Arabidopsis_thaliana.html

# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("BSgenome.Athaliana.TAIR.04232008")

#lnames = load("Athal_Infected_AnnotFile.RData")
data_annot <- read.csv(here("data","geneGO_metadata.csv"), header = TRUE, sep = '\t');
names(data_annot)
data_annot[1:10,1:5]
annot <- data_annot[,2:5]      #subset keeping certain columns
names(annot)
annot[1:10,1:4]

#Genes del modulo que deseamos anotar 
genes = names(athalData3)
length(genes)
genes2annot = match(genes, annot$Genes)
genes2annot[1:100]
#The following is the number or genes without annotation:
sum(is.na(genes2annot))
# Should return 0.
#save(genes, annot, genes2annot, file="Athal_Infected_AnnotFile.RData")

#The modules will be ordered by their significance for <Bc_24> or <Ch_22>, with the most significant ones to the left.

# Create the starting data frame (table)
geneInfo0 = data.frame(GenesID = genes,
                       GO_term = annot$GO_term[genes2annot],
                       GO_ID = annot$GO_id[genes2annot],
                       TAIR_ID = annot$TAIR_id[genes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

dim(geneInfo0)
geneInfo0[1:10,2:7]


# Order modules by their significance (GS) for Ch_22
modOrder = order(-abs(cor(MEs, Bc_24, use = "p")));
# Order modules by their significance (GS) for Ch_22
# modOrder = order(-abs(cor(MEs, Ch_22, use = "p")));

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

# Bc_24: Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Bc_24hpi));
# Ch_22: Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
# geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Ch_22hpi));

geneInfo = geneInfo0[geneOrder, ]
dim(geneInfo)
geneInfo[1:10,2:12]

write.csv(geneInfo, file = "results/geneInfo_GSignificance_Bc24hpi_simulation5000.csv")
#view(geneInfo)

