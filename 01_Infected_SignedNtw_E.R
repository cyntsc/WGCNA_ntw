######################################################################################
###
### WGCNA: Signed-Ntw with dynamic cut off (pearson) 
###
### Step 1/4: TOM calculation and adjacency matrix ( SIGNED NTW / Pearson corr)
###
### Made by: Cynthia Soto 
### Date: March 23, 2021
### L.Update: April 18, 2022
###
### DATA ASSUMPIONS:
### 1) Data are RNA-Seq expression profiles (raw-counts) from 14 plants (Arabidopsis) infected by some pathogen fungus.
### 2) Raw counts were estimated with HT-Seq/Solomon.
### 3) Data were scaled to log2(TPM) standardized values.

#######################    CLEAR OBJECT LIST AND IMAGE CANVAS   #######################
rm(list = ls());
if(!is.null(dev.list())) dev.off()

#####################  LIBRARIES REQUIRED FOR RUNNING THE CODE   ######################
#is.element("preprocessCore", packages)
packages <- c( "WGCNA", "matrixStats", "splines", "dynamicTreeCut", "flashClust", "foreach", "doParallel", "Hmisc", "survival",
               "lattice", "here", "devtools","utils", "tidyverse", "BiocManager")
packagecheck <- match( packages, utils::installed.packages()[,1] )
packagestoinstall <- packages[ is.na( packagecheck ) ]
if( length( packagestoinstall ) > 0L ) {
  utils::install.packages( packagestoinstall,repos = "http://cran.csiro.au"
  )
} else {
  print( "All requested packages already installed" )
}
# suppressPackageStartupMessages() use a special #+ comment to silence messages and/or warnings for the chunk that holds a chatty.
# We don’t want to suppress messages and warnings, in general, because they are an important part of the reprex. 
for( package in packages ) {
  suppressPackageStartupMessages(
    library( package, character.only = TRUE, quietly = TRUE )
  )
}
#BiocManager::install(c("GO.db", "preprocessCore", "impute")); 
#source("http://bioconductor.org/biocLite.R") 
#biocLite("impute")

##########################################################################################################
#
#                                   IMPORT LIBRARIES
#
##########################################################################################################
library(WGCNA);             # The package for weighted correlation network analysis. 
library(tidyverse);         # To processing data.
library(dynamicTreeCut);    # For detection of clusters in hierarchical clustering dendrograms.
library(flashClust);        # Fast implementation of hierarchical clustering.
library(lattice);           # High-level data visualization tool, with an emphasis on multivariate data. 
library(survival);          # In applied statistics, survival analysis studies the random processes related to the death of living organisms and the failure of mechanical systems (Kaplan-Meier and Aalen-Johansen curves, Cox models and accelerated failure time parametric models).
library(Hmisc);             # Functions for high-level graphics, functions for computing sample size and power, simulation, importing and annotating datasets, imputing missing values, etc.)
library(here);              # Set the top level of your project folder as “here” and to specify where things live relative to that location.   

#getwd();        
here::here();     # Top level dir: /data/run/cyntsc/Project_athal_wgcna

#############  Allow multi-treads depending the number of processors available in your PC ############### 
# If you are a linux user: $nproc o $lscpu -e
allowWGCNAThreads();
ALLOW_WGCNA_THREADS = 8;                # MSU Server=30; if not given, the number of processors online (as reported by system configuration) will be used
# Initial variables 
options(stringsAsFactors = FALSE);    # indicates whether strings in a data frame should be treated as factor variables or as just plain strings
enableWGCNAThreads();

##########################################################################################################
#
#                   LOAD DATA AND META-DATA; PREPROCESS AND FORMAT DATA
#
##########################################################################################################
# LOAD DATA 
athalData3 <- read.csv(here("data", "matrix_E_infected_simulated5000.csv"), header=TRUE,  sep=',');
athalData3[1:5,1:6];
dim(athalData3);  
# rearrange columns to match meta-data
col_order <- c("Genes", "Bc12",   "Bc12.1", "Bc18",   "Bc18.1", "Bc24",
               "Bc24.1", "Ch22",   "Ch22.1", "Ch22.2", "Ch22.3",
               "Ch40",   "Ch40.1", "Ch40.2", "Ch40.3");
athalData3 <- athalData3[, col_order];
## You see that genes are listed in a column named "Genes" and samples are in columns
athalData3[1:5,1:6];

#################    Create variables and format data to match the format WGCNA needs #################### 
row.names(athalData3) = athalData3$Genes;   # Pull the gene names 
SubGeneNames=rownames(athalData3);
typeof(SubGeneNames);
SubGeneNames[0:10];

athalData3$Genes = NULL                      # Set genes column as index row
athalData3 = as.data.frame(t(athalData3))    # transpose dataset, now samples are rows and genes are columns
athalData3[,1:10]
rownames(athalData3)                         # now the row names are the samples

################# Run this to check if there are gene outliers   ###################
gsg = goodSamplesGenes(athalData3, verbose = 3)
gsg$allOK
## If the last statement returns TRUE, all genes have passed the cuts. If not, is removed the offending genes and samples from the data.
if (!gsg$allOK)
{if (sum(!gsg$goodGenes)>0)
  printFlush(paste("Removing genes:", paste(names(athalData3)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(athalData3)[!gsg$goodSamples], collapse=", ")))
  athalData3= athalData3[gsg$goodSamples, gsg$goodGenes]
}

# Load meta-data for linking the dataset: TRAITMENTS BY SAMPLE
datTraits <- read.csv(here("data", "Traits_Infected_sample.csv"), header=TRUE,  sep=',')   
datTraits[1:5,2:7]
datTraits = subset(datTraits, select = -c(accesion, B_cinerea, C_higginsianum));
datTraits[1:5,]
#dim(datTraits)
# Form a data frame analogous to expression data that will hold the infection traits.
rownames(datTraits) = datTraits$ID
datTraits$ID = NULL                                  # Set ID sample column as index row
table(rownames(datTraits)==rownames(athalData3))     # should return TRUE if datasets align correctly, otherwise your names are out of order

# Load meta-data for linking the dataset: TRAITMENTS BY FUNGI INTERACTION AND HPI
datTraits2 <- read.csv(here("data", "Traits_Infected_hpi.csv"), header=TRUE,  sep=',')   
datTraits2 = subset(datTraits2, select = -c(accesion, B_cinerea, C_higginsianum));
datTraits2[1:5,]
dim(datTraits2)
## form a data frame analogous to expression data that will hold the infection traits.
rownames(datTraits2) = datTraits2$ID
datTraits2$ID = NULL                 
table(rownames(datTraits2)==rownames(athalData3)) 

#save(athalData3, datTraits, datTraits2, file="Athal_Infected_MatrixE_simulated5000.RData")    
#load("Athal_Infected_MatrixE_simulated5000.RData")

##########################################################################################################
#
#             ANALYSIS OF SCALE FREE TOPOLOGY FOR SOFT-THRESHOLDING
#        know more: https://rdrr.io/cran/WGCNA/man/pickSoftThreshold.html 
#        
##########################################################################################################

# Choosing a soft-threshold to fit a scale-free topology to the network
# pickSoftThreshold():
#                  1) Help to pick up an appropriate soft-thresholding power for network construction
#                  2) The function calculates the similarity of columns (genes) in datExpr by calling the function given in corFnc or distFnc,
#                     transforms the similarity according to type and raises it to power, resulting in a weighted network adjacency matrix.
powers = c(c(1:10), seq(from =10, to=20, by=2)) 
# Another example: 
# powers = c(c(1:10), seq(from =10, to=30, by=1)); 

sft = pickSoftThreshold(athalData3,
                        dataIsExpr = TRUE,
                        powerVector = powers,
                        corFnc = cor,                             # cor=Pearson correlation or any function returning values between -1 and 1 can be used. 
                        corOptions = list(use = 'p'),             # Almost all lists in R internally are Generic Vectors, whereas traditional dotted pair lists (as in LISP) remain available but rarely seen by users (except as formals of functions). 
                        verbose = 5,                              # from 0 to 5
                        networkType = "signed");                  # "unsigned", "signed", "signed hybrid"
# A signed ntw-preserve the natural continuity of the correlation (+ / -)
#          contrary to what is preserved in an unsigned ntw
#warnings()
#write.table(sft, file = "results/SFT_corr_Athal_Infected_MatrixE.txt", sep = ",", quote = FALSE, row.names = T)

# Plot the results
sizeGrWindow(5, 5)
par(mfrow = c(1,3));     # 1 panel in y / 3 in x  
cex1 = 1.5 # 0.9;
# Scale-free topology fit index as a function of the soft-threshold power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit, Signed R^2",
     type="n", main = paste("Scale independence for plant infected \n (Matrix E)"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
# Red line corresponds to using an R^2 cut-off
abline(h=0.60,col="red");  
abline(h=0.70,col="blue");  
abline(h=0.80,col="brown");  

# Mean connectivity as a function of the soft-threshold power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=160,col="blue");

# Median connectivity as a function of the soft-threshold power
plot(sft$fitIndices[,1], sft$fitIndices[,6],xlab="Soft Threshold (power)",ylab="Median Connectivity", type="n",main = paste("Median connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,6], labels=powers, cex=cex1,col="red")
abline(h=130 ,col="green");

softPower = 14; 

#############################################################################################
#
#   Generating adjacency and TOM diss/similarity matrices based on the selected soft-power
#
#############################################################################################

#calculate the adjacency matrix
adjacency = adjacency(athalData3,type = "signed",   
                     power = softPower);
dim(adjacency)
head(adjacency)
adjacency[1:5,1:5]
adjacency[2:10]
adjacency[3:10]

# translate the adjacency into TOM and calculate the corresponding dissimilarity
# this action minimize the effects of noise and spurious associations
TOM = TOMsimilarityFromExpr(adjacency,                         
                          TOMType = "signed", 
                          power = softPower);
TOM[1:5,1:5]
dim(TOM)
dissTOM = 1-TOM
dissTOM[1:5,1:5]

#save(dissTOM, softPower, file="Athal_Infected_DissTOM_MatrixE_simulated5000.RData") 
#load("Athal_Infected_DissTOM_MatrixE_simulated5000.RData")

#########################################################################################
##
##    Generate Modules 
##
#########################################################################################

## Generate a clustered gene tree with flashCLust or hclust
## hclust make a hierarchical cluster analysis on a set of dissimilarities
#geneTree = flashClust(as.dist(dissTOM),method="average")      # average, ward, single, complete, mcquitty, median or centroid
geneTree = hclust(as.dist(dissTOM), method="average")

## Plot the results
sizeGrWindow(9, 12)
#plot the resulting clustering tree (dendrogram)
plot(geneTree, xlab="Gene clusters", 
     main="Gene clustering on TOM−based dissimilarity for A.thaliana infected (Matrix E)", 
     cex.main = 1,
     cex.axis = 1,
     cex=0.3)   # sub="Signed-Ntw (average)" , labels= FALSE, hang=0.04

########################################################################################
##
##    Module identification using dynamic tree cut 
##
########################################################################################

## Set the minimum module size
minModuleSize = 20;
diff(geneTree$height)

## Method1: This chunck is the function for pruning of HClust with the method tree.
# dynamicMods = cutreeDynamic(dendro = geneTree,
#                            method="tree",
#                            deepSplit=TRUE,
#                            #cutHeight = 0.99,
#                            pamRespectsDendro= FALSE, 
#                            minClusterSize = minModuleSize);

## Method2: This chunck is the function for pruning of HClust with the method disstM
## We chose dissTOM for this analysis because all genes are asigned to some module *** ###
dynamicMods = cutreeDynamic(dendro= geneTree, 
                            distM = dissTOM,        ## Method "hybrid". The distance matrix used as input to hclust
                            deepSplit=2,            ## For method "hybrid",  range 0 to 4. If TRUE, method will favor sensitivity and produce more smaller clusters. When FALSE, there will be fewer bigger clusters.
                            pamRespectsDendro= FALSE,     ## La etapa PAM respetará el dendrograma en el sentido de que los pequeños clusters solo se asignarán a los clusters que pertenezcan a la misma rama 
                            minClusterSize = minModuleSize)

## when cutHeight not given, for method=="tree" it defaults to 0.99, for method=="hybrid" it defaults to 99% of the range between the 5th percentile and the maximum of the joining heights on the dendrogram

## gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)
# write.table(table(dynamicMods), file = "results/dynamicMods_Athal_Infected_MatrixE.txt", sep = ",", quote = FALSE, row.names = F)

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)      # colors code is assigned to plot
sort(table(dynamicColors), decreasing = TRUE)   # get the number of genes by color
# write.table(sort(table(dynamicColors), decreasing = TRUE), file = "results/dynamicColors_Athal_Infected_MatrixE.txt", sep = ",", quote = FALSE, row.names = F)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.05,
                    addGuide = TRUE, guideHang = 0.1,
                    cex.main = 1,
                    main = "Gene dendrogram and module colors for Athal infected \n Matrix E")

#save(dynamicMods, dynamicColors, geneTree, file="Athal_Infected_DynamicMods_MatrixE_simulated5000.RData")
#load("Athal_Infected_DynamicMods_MatrixE_simulated5000.RData")

##################################################################################################
#
##                      Plot the topological overlap matrix
#
##################################################################################################

## Raise the dissimilarity matrix to the power of 25 to bring out the module structure
sort(table(dynamicColors), decreasing = TRUE)
restGenes = (dynamicColors != "grey")                   #restGenes= (dynamicColors == "mediumpurple2")

## set the DIAGONAL of the dissimilarity to NA 
diag(dissTOM) = NA;
sizeGrWindow(7,7)
## Graphical representation of the TOM using a heatmap plot combined with the corresponding hclust dendrogram and module colors.
TOMplot(dissTOM, 
        geneTree, 
        as.character(dynamicColors[dynamicColors]),
        ColorsLeft = dynamicColors,
        terrainColors = TRUE,
        main = "TOM for A.thaliana infected / Matrix E");

##################################################################################################
#
#                          Merging of modules whose expression profiles are very similar
#
##################################################################################################

# Calculate eigengenes (PCA)
MEList= moduleEigengenes(athalData3, colors= dynamicColors,
                         excludeGrey = TRUE,
                         softPower = softPower)
# Here are the ME tables by color module
MEs = MEList$eigengenes
MEs[1:5,]
MEs$MEgreen4

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)          
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss),method= "average")
#save(MEs, MEDiss, METree, file= "Athal_Infected_Module_Identification_MatrixE_simulated5000.RData")
#load("Athal_Infected_Module_Identification_MatrixE_simulated5000.RData")

# Clustering of module eigengenes
sizeGrWindow(6,14)

## plots tree showing how the eigengenes cluster together
plot(METree, main= "Clustering of module eigengenes for A.thaliana infected (Matrix E)", 
     cex.main = 1.2,
     xlab= "",      #col.lab="red", cex.lab=1.5,
     ylab = "",
     cex = 0.5,     #labels = FALSE
     sub= "");   #cex.main = 1, cex.lab = 1, cex.axis = 1    

abline(h=0.50,col="red"); 

##########################################################################################
##
## set a threshold for merging modules. Set the threshold in the MEDissThres=0.0 variable
##
##########################################################################################

## We choose a height cut of 0.10, corresponding to correlation of 0.90, to merge

MEDissThres = 0.10
## Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
## Call an automatic merging function
merge = mergeCloseModules(athalData3, dynamicColors, 
                          cutHeight= MEDissThres, 
                          verbose =3)
merge$dendro
merge$oldDendro
merge$newMEs
dim(merge$newMEs)
merge$cutHeight

## The merged module colors
mergedColors = merge$colors
## Eigengenes of the new merged modules
mergedMEs = merge$newMEs
mergedMEs$MEblack[1:5]

length(mergedMEs)   ## number if merged MEs
sort(table(mergedColors), decreasing = TRUE)
write.table(sort(table(mergedColors), decreasing = TRUE), file = "results/mergedColors_5000.txt", sep = ",", quote = FALSE, row.names = F)

## INCLUDE THE NEXT LINE TO SAVE TO FILE
pdf(file="results/mergedColors_5000.pdf")

## plot dendrogram with module colors below it
plotDendroAndColors(geneTree, main= "Clustering of MEs corresponding to correlation of 0.90 merged \nA thaliana infected (Matrix E)", 
                    cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), 
                    dendroLabels = FALSE, 
                    cex.main = 1,
                    hang=0.03, addGuide= TRUE, 
                    guideHang=0.05)

# Reasign the new variables to the old one 
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));          # standardColors(50)
moduleLabels = match(mergedColors, colorOrder)-1     # moduleColors
MEs = mergedMEs
#INCLUE THE NEXT LINE TO SAVE TO FILE
#dev.off()

#save(MEs, moduleLabels, moduleColors, file= "Athal_Infected_MergedMods_MatrixE_simulated5000.RData")
#load("Athal_Infected_MergedMods_MatrixE_simulated5000.RData")

## Write modules no merged in the path location
length(dynamicColors)
module_colors = setdiff(unique(dynamicColors), "grey")    ## exclude the grey module reserved for not assigned genes
length(module_colors)
for (color in module_colors){
  module = SubGeneNames[which(dynamicColors==color)]
  write.table(module, paste("results/module_",color, ".txt",sep=""), 
              sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
}

## Write MEs merged in the path location
length(mergedColors)
module_colors = setdiff(unique(mergedColors), "grey")     ## exclude the grey module reserved for not assigned genes
length(module_colors)
#module_colors
for (color in module_colors){
  module=SubGeneNames[which(mergedColors==color)]
  write.table(module, paste("results/merged_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

##################################################################################################
#
##Calculation of the topological overlap matrix for merged modules 
#
##################################################################################################

sort(table(mergedColors), decreasing = TRUE)
##set the DIAGONAL of the dissimilarity to NA 
#diag(dissTOM) = NA;

TOMplot(dissTOM,
        geneTree,
        as.character(mergedColors[mergedColors]),
        ColorsLeft = mergedColors,
        terrainColors = TRUE,
        main = "TOM for A.thaliana infected / Matrix E");

##################################################################################################
##
##       Correlate traits (Quantifying module–trait associations)
##       Identify modules that are significantly associated with the measured biology traits
##       *** Here we use correlation by HPI 
##
##################################################################################################

# Define the number of genes and samples
nGenes = ncol(athalData3)
nSamples = nrow(athalData3)
## Recalculate MEs with color labels.   
MEs0 = moduleEigengenes(athalData3, moduleColors)$eigengenes  # Calculates module eigengenes (1st principal component) of modules in a given single dataset
MEs = orderMEs(MEs0)

# Link the external values to the correlation matrix
#moduleTraitCor2 = cor(MEs, datTraits2, use= "p")
moduleTraitCor2 = cor(MEs, datTraits, use= "p")

# Calculates Student asymptotic p-value for given correlations.
moduleTraitPvalue2 = corPvalueStudent(moduleTraitCor2, nSamples)   
# Asymptotic p-values are useful for large sample sizes when the calculation of an exact p-value is too computer-intensive.
# We color code each association by the correlation value

# Print correlation and p-value heatmap between modules and traits
#textMatrix2= paste(signif(moduleTraitCor2, 2), " (", signif(moduleTraitPvalue2, 1), ")", sep= "")
# Print just the correlation heatmap between modules and traits.
# textMatrix2= paste(signif(moduleTraitCor2, 2), " / (", signif(moduleTraitPvalue, 1), ")", sep= "")
textMatrix2= paste(signif(moduleTraitCor2, 2))
# Print "" if too many modules
#textMatrix= ""
dim(textMatrix2) = dim(moduleTraitCor2)
#par(mar= c(6, 8.5, 3, 3))
sizeGrWindow(9, 12)
par(mar= c(3.5, 10, 2, 1))    # margen inferior - y(izq) + margen superior - y(der)

# Display the corelation values with a heatmap
#INCLUE THE NEXT LINE TO SAVE TO FILE
#pdf(file="heatmap.pdf")
labeledHeatmap(Matrix= moduleTraitCor2,
               #xLabels= names(datTraits2),
               xLabels= names(datTraits),
               yLabels = names(MEs),
               yLabelsPosition = "left",
               yColorWidth = 0.1,
               cex.lab.y = 0.8,
               cex.lab.x = 1,
               ySymbols= names(MEs),
               #colorLabels= FALSE,
               colors= blueWhiteRed(50), #blueWhiteRed(50), greenBlackRed()
               textMatrix= (textMatrix2),
               setStdMargins= FALSE,
               cex.text= 0.6,
               #zlim= c(-1,1),
               main= paste("Module-trait by h.p.i"))

###################################################################################################
##                       
##                       Classical multidimensional scaling (MDS) plots
##                       
###################################################################################################                      

# Classical multidimensional scaling (MDS) of a data matrix, is also known as principal coordinates analysis (Gower, 1966)   
# We want to represent the distances among the objects in a parsimonious (and visual) way (i.e., a lower k-dimensional space). 
par(mfrow = c(1, 1))

cmd1 = cmdscale(as.dist(dissTOM), 3)
plot(cmd1, col = moduleColors, main = "MDS plot", xlab = "Scaling Dimension 1", 
     ylab = "Scaling Dimension 2")
 

