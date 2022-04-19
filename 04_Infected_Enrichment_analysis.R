#################################################################################################
###
### Goal: Enrichment analysis
###
### Made by: Cynthia Soto 
### Date: Feb 18, 2022
### Latest update: Apr 18, 2022
###
#################################################################################################

# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(athalData3, power = 9);
# Read in the annotation file
# annot = read.csv(file = "GeneAnnotation.csv");
annot = read.csv(here("data", "geneGO_metadata.csv"), header=TRUE, sep='\t')
annot[1:5,]
dim(annot)
names(annot)


# Select module
module = "darkmagenta";
# Select module genes
genes = names(athalData3);
inModule = (moduleColors == module);
modgenes = genes[inModule];

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modgenes, modgenes);
# exportNetworkToVisANT() export network data in format readable by VisANT can read
# VisANT software is available from http://www.visantnet.org/visantnet.html/
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$Genes, annot$GO_term));

vis[1:10,1:5];
dim(vis);

nTop = 30;
# softConnectivity() calculates connectivity of a weighted network
IMConn = softConnectivity(athalData3[, modgenes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$Genes, annot$GO_term) )
vis[1:10,1:5]
dim(vis)

##########################       Exporting Cytoscape     #######################################

# Recalculate topological overlap if needed
# Select modules
modules = c("brown", "red");
modules = c("darkmagenta");

# Select module genes
genes = names(athalData3)
inModule = is.finite(match(moduleColors, modules));
modgenes = genes[inModule];
modGenes = annot$GO_term[match(modgenes, annot$Genes)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];

dimnames(modTOM) = list(modgenes, modgenes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modgenes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])

