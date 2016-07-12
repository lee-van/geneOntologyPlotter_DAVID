###plot GO enrichment for transcripts with HAMR mods (GMUCT data)

library(RColorBrewer, quietly = T)
library(gplots, quietly = T, warn.conflicts = F)

args <- commandArgs(TRUE)
input = args[1]
upperbound = args[2]
lowerbound = args[3]
transform = args[4]
output = args[5]


#load parsed DAVID ouptut
summary = read.table(input, header = T, sep = "\t")

#if multiple GO IDs share the same GO term, collapse and and retain only the first instance. Most likely result from alternate ids for same term, so should share same stats
summary = summary[which(!duplicated(summary[,2])),]
rownames(summary) = summary$name
summary = summary[-c(1:2)]

#sort by desired level and catenate
if (upperbound != "NA"  & lowerbound != "NA") {
  lowerbound = as.integer(lowerbound)
  upperbound = as.integer(upperbound)
  summary_mf = summary[which(summary$level >= lowerbound & summary$level <= upperbound & summary$namespace == "molecular_function"), ,drop =F]
  summary_bp = summary[which(summary$level >= lowerbound & summary$level <= upperbound & summary$namespace == "biological_process"), ,drop =F]
  summary_cc = summary[which(summary$level >= lowerbound & summary$level <= upperbound & summary$namespace == "cellular_component"), ,drop =F]
} else {
  summary_mf = summary[which(summary$namespace == "molecular_function"), ,drop =F]
  summary_bp = summary[which(summary$namespace == "biological_process"), ,drop =F]
  summary_cc = summary[which(summary$namespace == "cellular_component"), ,drop =F]
}

#heatmaps
bp_matrix = summary_bp[,-(1:2)]
bp_matrix = data.matrix(bp_matrix)
bp_matrix[which(is.na(bp_matrix))] = 1
mf_matrix = summary_mf[,-(1:2)]
mf_matrix = data.matrix(mf_matrix)
mf_matrix[which(is.na(mf_matrix))] = 1
cc_matrix = summary_cc[,-(1:2)]
cc_matrix = data.matrix(cc_matrix)
cc_matrix[which(is.na(cc_matrix))] = 1

if (transform == "NA") {
} else if (transform == "neglog10") {
  bp_matrix = -log10(bp_matrix)
  mf_matrix = -log10(mf_matrix)
  cc_matrix = -log10(cc_matrix)
} else if (transform == "log10") {
  bp_matrix = log10(bp_matrix)
  mf_matrix = log10(mf_matrix)
  cc_matrix = log10(cc_matrix)
} else if (transform == "neglog2") {
  bp_matrix = -log2(bp_matrix)
  mf_matrix = -log2(mf_matrix)
  cc_matrix = -log2(cc_matrix)
} else if (transform == "log2") {
  bp_matrix = log2(bp_matrix)
  mf_matrix = log2(mf_matrix)
  cc_matrix = log2(cc_matrix)
} 

hmcol<-colorRampPalette(c("white","orange", "red"), bias = 1)(256)
pdf(file = output)
if (nrow(bp_matrix) > 1 & ncol(bp_matrix) > 1) {
  heatmap.2(bp_matrix, main = "                     Biological Process Term Enrichment\n Across Structure Deciles", 
            Colv = NA, dendrogram = 'row', colsep = seq(1, ncol(bp_matrix)), 
            rowsep = seq(1, nrow(bp_matrix)), sepwidth=c(.01,.01), col = hmcol, scale="none",
            margins=c(5,15), cexCol = .5, trace = "none", key = T, keysize = 1.5,
            cexRow = .5, density.info = "none")
}
if (nrow(mf_matrix) > 1 & ncol(mf_matrix) > 1) {
  heatmap.2(mf_matrix, dendrogram = 'row', main = "                     Molecular Function Term Enrichment\n Across Structure Deciles",
          Colv = NA, colsep = seq(1, ncol(mf_matrix)), 
          rowsep = seq(1, nrow(mf_matrix)), sepwidth=c(.01,.01), col = hmcol, scale="none", 
          margins=c(5,18), cexCol = .5, trace = "none", keysize = 1.2, cexRow = .5, density.info = "none")
}
if (nrow(cc_matrix) > 1 & ncol(cc_matrix) > 1) {
  heatmap.2(cc_matrix, dendrogram = 'row', main = "                     Cellular Component Term Enrichment\n Across Structure Deciles",
          Colv = NA, colsep = seq(1, ncol(cc_matrix)), 
          rowsep = seq(1, nrow(cc_matrix)), sepwidth=c(.01,.01), col = hmcol, scale="none", 
          margins=c(5,18), cexCol = .5, trace = "none", keysize = 1.2, cexRow = .5, density.info = "none")
}
invisible(dev.off())

