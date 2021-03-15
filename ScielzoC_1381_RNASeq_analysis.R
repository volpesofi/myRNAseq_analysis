#--------------- Prj set up ------------------ 
# set working diretory
prj = 'ScielzoC_1381_RNASeq' 
PI = 'Scielzio'
setwd(paste("/Users/tascini.annasofia/Dropbox (HSR Global)/WORKSPACE/",
            PI,"/",prj,
            "/7_bioinfo/",
            sep=''))

#--------------- load libraries --------------- 
suppressMessages(library("edgeR"))
suppressMessages(library(data.table))
suppressMessages(library("DESeq2"))
suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library("RColorBrewer"))
suppressMessages(library(enrichR))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library('VennDiagram'))
suppressMessages(library(venn))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))
suppressMessages(library("viridis"))
suppressMessages(library('pals'))
suppressMessages(library(RColorBrewer))
suppressMessages(library(wesanderson))
suppressMessages(library(patchwork))
suppressMessages(library(magick))
suppressMessages(library('philentropy'))
suppressMessages(library('IntClust'))
suppressMessages(library('pheatmap'))
suppressMessages(library(assertr))
suppressMessages(library("remotes"))
suppressMessages(library(GeneOverlap))
suppressMessages(library('stringr'))
suppressMessages(library(ggrepel))

#--------------- Counts and Nreplica ---------------
# Min number of replica in each group
Nreplica= 3

# path to the count file
# rsync -rlv tascini.annasofia@srhpclogin01.ihsr.dom:/beegfs/scratch/ric.cosr/ric.scielzo/ScielzoC_1381_RNASeq/samples/all/runs/all/fastq/merge-by-read/trimmed/trimmomatic/mapped/STAR/merged/featureCounts/merged/all.counts.gz .

filecount = "./all.counts.gz"

# import counts
annotation <- c('GeneID','Chr','Start','End','Strand','Length')
fCounts <- read.delim(file=filecount, header=TRUE, check.names = F)
fCountsData <- fCounts[
  , 
  -which(
    tolower(names(fCounts))
    %in% 
      tolower(annotation))]

fCountsAnnotation <- fCounts[
  , 
  which(
    tolower(names(fCounts))
    %in% 
      tolower(annotation))]

geneidColname <- 'Geneid'
geneidIdx <- which(tolower(annotation) %in% tolower(geneidColname))
rownames(fCountsData) <- fCounts[[geneidIdx]]

# ----------------- Metadata -------------------
metadata = data.frame(SampleID = sort(colnames(fCountsData)), 
                      condition = sort((str_sub(colnames(fCountsData),1,-2))))
row.names(metadata) <- metadata$SampleID
metadata
# Reordering counts matrix to have samples ordered as in metadata
# it assumes that the first column in metadata is sample name
fCountsData <- fCountsData[,match(metadata[,1], colnames(fCountsData))] 

#------------------------- biotypes ------------------------------------
#BiocManager::install("GEOquery")
source('./plot-biotypes.R')
biotypesFile = paste('/lustre1/workspace/Scielzo/ScielzoC_1205_2D_3D_cell_line/',
                     'gencode.v31.basic.annotation.BIOTYPES.DICTIONARY.gz', sep='')
# also saved in new cluster in 
# /beegfs/scratch/ric.cosr/ric.scielzo/ScielzoC_1205_2D_3D_cell_line
biotypes_function(
  countsFile    = filecount,
  biotypesFile  = biotypesFile,
  pngFolder     = 'Biotypes/',
  minSamples    = 3,
  filterExp     = TRUE,
  useRpkm       = TRUE,
  plotPie       = TRUE,
  sglSamplePlot = TRUE,
  writeTable    = TRUE,
  perc2plot     = 0.001,
  useGgplot     = TRUE)

#------------ save counts in excel file: rawcount, filtered, cmp, rpkm ------------------------
SAVE_variable <- list()
filename_xls <- paste('COUNTS_',prj,'.xlsx', sep='')
variable2save_names <- c('all_counts', 'expGenes_counts','expGenes_LogCPM', 'expGenes_LogRPKM')

# all_counts
y <- DGEList(counts=fCountsData, genes = fCountsAnnotation)
SAVE_variable[[variable2save_names[1]]] <- as.data.frame(y$counts)

# expGENES_counts
keep <- rowSums(cpm(y)>1)>=Nreplica
table(keep)
yf <- y[keep,]
SAVE_variable[[variable2save_names[2]]] <- as.data.frame(yf$counts)


#CPM
SAVE_variable[[variable2save_names[3]]] <- as.data.frame(cpm(y, log=T))[keep,]

#RPKM
SAVE_variable[[variable2save_names[4]]] <- as.data.frame(rpkm(y, log=T, gene.length =y$genes$Length))[keep,]

write.xlsx(SAVE_variable,
           file = filename_xls, 
           row.names = T,
           asTable = T, 
           sheetName =variable2save_names)

#------------ 500 most variant genes ------------
y <- DGEList(counts=fCountsData, genes = fCountsAnnotation)
# rpkm
fCountsRPKM = rpkm(y, log=T, gene.length =y$genes$Length)
# filter for expression
keep <- rowSums(cpm(y)>1)>=Nreplica
yf <- y[keep,]

# calculate 500 most variant genes and 
# save their counts in fCountsRPKMTOP
N=500
vary <- apply(fCountsRPKM[keep,],1,var)
vary_s <- sort(vary, decreasing = T)
TOP_N <- names(vary_s[1:N])
yTOP <-  y[TOP_N,]
fCountsRPKMTOP <- fCountsRPKM[TOP_N,]

######### control genes ############
gene_control_file = paste(
  '/lustre2/scratch/bioinfotree/common/bioinfotree/prj/',
  'MalteccaF_926_Corteccia_Ippocampo/local/share/data/',
  'standard.controls.hsapiens.txt', sep='')
gene_control = read.table(gene_control_file)
gene_control.r = as.character(gene_control$V1)
rpkm.control = fCountsRPKM[which(rownames(fCountsRPKM) %in% gene_control.r),]
#### heatmap control genes
colors <- colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)
pheatmap::pheatmap(rpkm.control,
                   cluster_rows = F,
                   cluster_cols = F,
                   main = 'Heatmap of housekeeping genes - RPKM',
                   show_rownames = T,
                   show_columnames = F,
                   fontsize = 12, fontsize_row = 10, fontsize_col = 14, 
                   display_numbers = F, 
                   col=colors,
                   filename = 'Heatmap_controlGenes.pdf',
                   width = 10, height = 7)

######### control genes ############
gene_control_file = './standard.controls.hsapiens.txt'
gene_control = read.table(gene_control_file)
gene_control.r = as.character(gene_control$V1)
rpkm.control = fCountsRPKM[which(rownames(fCountsRPKM) %in% gene_control.r),]
#### heatmap control genes
colors <- colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)
pheatmap::pheatmap(rpkm.control,
                   cluster_rows = F,
                   cluster_cols = F,
                   main = 'Housekeeping genes - RPKM',
                   show_rownames = T,
                   show_columnames = F,
                   cellheight = 20,
                   cellwidth = 20,
                   fontsize = 10, fontsize_row = 10, fontsize_col = 10, 
                   display_numbers = F, 
                   col=colors,
                   filename = 'Heatmap_controlGenes.pdf')

########### explorative PCA plots ##########
#PCA parameters
pcx = 1
pcy = 2
centering = TRUE
scaling = TRUE
# PCA
pca = prcomp(t(fCountsRPKMTOP), center=centering, scale=scaling)
var = round(matrix(((pca$sdev^2)/(sum(pca$sdev^2))), ncol=1)*100,1)
score = as.data.frame(pca$x)

# plot paramters
xlab = paste("PC", pcx, " (",var[pcx],"%)", sep="")
ylab = paste("PC", pcy, " (",var[pcy],"%)", sep="")
cum = var[pcx]+var[pcy]
names = rownames(pca$x)


#plot
pca_plot <- list()

i=0
for (col in 2:dim(metadata)[2]) {
  score$factor <-  metadata[,col]
  i=i+1
  pca_plot[[i]] <- ggplot(score, aes(x=score[,pcx], y=score[,pcy], color=factor))+
    geom_point(size=5)+
    labs(x=xlab, y=ylab, title=paste("PC",pcx," vs PC",pcy," scoreplot",sep="")) + 
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22, hjust =.5), 
          axis.title.x = element_text(face = "bold", color = "black", size = 24),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
          axis.title.y = element_text(face = "bold", color = "black", size = 24),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) 
}

pca_plot

                       
########### PCA plot for report ##########
score$condition = metadata$condition
score$sampleID = metadata$SampleID
#stringr::str_wrap(V1, 15)
#str_sub(x <- sampleID,7,-1)
pca<- ggplot(score, aes(x=score[,pcx], y=score[,pcy], 
                        color=condition))+
  geom_label_repel(data= score, aes(x=score[,pcx], y=score[,pcy], 
                                    color=condition, label = str_sub(x <- sampleID,10,-1)),
                   size = 6,  box.padding = unit(.5, "lines"), point.padding = unit(0.1, "lines"),
                   segment.color = 'grey50') +
  geom_point(size= 8)+
  labs(x=xlab, y=ylab, title=paste("PC",pcx," vs PC",pcy," scoreplot",sep="")) + 
  geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
  geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
  theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22, hjust =.5), 
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
        axis.title.y = element_text(face = "bold", color = "black", size = 24),
        legend.text = element_text(face = "bold", color = "black", size = 18),
        legend.position="right", 
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  scale_color_manual(values = brewer.pal(n = 3, name = "Set1"))  

pca
pdf('PCA_top500rpkm.pdf', width = 10, height = 8.5)
pca
dev.off()

#------------------- Heatmap for report ###################
annotation_column <- as.data.frame(metadata[,2:(dim(metadata)[2])])
colnames(annotation_column) = "condition" 

mycolors_c <- brewer.pal(n = 3, name = "Set1");     names(mycolors_c) = levels(annotation_column$condition)
ann_colors = list(
  condition = mycolors_c
)

row.names(annotation_column) <- metadata[,1]

crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)


HP <- pheatmap::pheatmap(fCountsRPKMTOP,
                         scale = 'row',
                         annotation_col = annotation_column,
                         annotation_colors = ann_colors, 
                         cluster_rows = T, 
                         cluster_cols = T,                       
                         #main = 'Heatmap: 500 most variable genes - RPKM',
                         show_rownames = F,
                         cutree_rows = 2,
                         cutree_cols = 3,
                         fontsize = 12, fontsize_row = 10, fontsize_col = 14, 
                         display_numbers = F, 
                         col=colors,
                         filename = 'Heatmap_500rpkm.pdf',
                         width = 10, height = 8.5)



#---------------- combine PCA and heatmap plots ----------------
A1 <- image_read_pdf('PCA_top500rpkm.pdf', density = 200)
A2 <- image_read_pdf('Heatmap_500rpkm.pdf', density = 200)
im = image_append(c(A1, A2))
im
image_write(im, path = "Clustering_composite.tiff", format = "tiff")

#------------------ DGE - DESeq object ------------------ 
dds <- DESeqDataSetFromMatrix(
  countData = fCountsData,
  colData  = metadata,
  design   = as.formula('~condition'))

filter <- rowSums(cpm(counts(dds)) >= 1) >= Nreplica
table(filter)
ddsFiltered <- dds[filter,]

ddsFiltered[['condition']] <- relevel(ddsFiltered[['condition']] , ref = 'PTsbasale')

dga <- DESeq(
  object = ddsFiltered,
  test = "Wald",
  fitType = "parametric",
  betaPrior = FALSE,
  minReplicatesForReplace = Inf)

plotDispEsts(dga)

#--------------- DGE - comparisons ------------------------------
contrasts = resultsNames(dga)[- which(resultsNames(dga) %in% 'Intercept')]
contrasts
alpha = 0.05

dgeResults <- list()
for (contrast in contrasts) {
  dgeResults[[contrast]] <- results(
    dga,
    name                 = contrast,
    cooksCutoff          = Inf,
    independentFiltering = TRUE,
    alpha                = alpha,
    pAdjustMethod        = "BH")
  print(summary(dgeResults[[contrast]]))
  # sorting gene list according to significance
  dgeResults[[contrast]] <- dgeResults[[contrast]][order(dgeResults[[contrast]]$pvalue, decreasing = F),]
}

# add additional comparison
v = "condition"
c1 = "PTs3Dday7"
c2 = "PTs2Dday7"
print(paste(c1,"_vs_",c2,sep=''))
dgeResults.tmp <- results(dga,
                          contrast             = c(v,c1,c2),
                          cooksCutoff          = Inf,
                          independentFiltering = TRUE,
                          alpha                = alpha,
                          pAdjustMethod        = "BH")
summary(dgeResults.tmp)
dgeResults[[paste(v,'_',c1,"_vs_",c2,sep='')]] <- dgeResults.tmp[order(dgeResults.tmp$pvalue, decreasing = F),] 


#------------------- DGE - Save results ---------------------------
# Save results
f="dgeResults"
dir.create(f, showWarnings=TRUE, recursive=TRUE)

lapply(
  names(dgeResults),
  function(x) write.table(
    data.table(
      data.frame(dgeResults[[x]]),
      keep.rownames=geneidColname),
    file.path(f, paste(x, ".tsv", sep="")),
    append=F,
    row.names=F,
    col.names=T,
    quote=F,
    sep="\t"))


dgeResults_table = list()    
dgeResults_table = lapply(
  names(dgeResults),
  function(x) 
    data.table(
      data.frame(dgeResults[[x]]),
      keep.rownames=geneidColname))

names(dgeResults_table) = names(dgeResults)


write.xlsx(dgeResults_table,
           file = paste(f,'/DGE_results.xlsx', sep=''), 
           row.names = F,
           asTable = T, 
           sheetName = str_sub(names(dgeResults),1,31)) 

#------------------ DGE - MAplot and Vulcano Plots ------------------
f="dgeResults"
n.label = 20
FDR = T
pvalue = 0.01
for (condition in names(dgeResults)) {
  results = as.data.frame(dgeResults[[condition]])
  results$DE = 'unm'
  if (!FDR) {
    if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]))>0) {
      results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]$DE = 'up'}
    if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]))>0) {
      results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < .1,]$DE = 'down' }
  } else {
    if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]))>0) {
      results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]$DE = 'up'}
    if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]))>0) {
      results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]$DE = 'down' }        
  }
  if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]))>0) {
    results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]$DE = 'SEQCup'}
  if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]))>0) {
    results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]$DE = 'SEQCdown' }
  if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]))>0) {
    results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]$DE = 'FDRup'}
  if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]))>0) {
    results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]$DE = 'FDRdown' }        
  
  results$DE <- factor(x = results$DE, levels = c("unm", "FDRdown","FDRup", 'SEQCdown','SEQCup'))
  mycolors = c('grey','dodgerblue4','darkred','dodgerblue2','coral'); names(mycolors) = levels(results$DE)
  results$DE2 = 'unm'; results[results$DE!='unm',]$DE2 = 'mod'
  results$DE2 <- factor(x = results$DE2, levels = c("mod","unm"))
  mysize = c(3,2); names(mysize) = levels(results$DE2)
  myalpha = c(1,0.2); names(mysize) = unique(results$DE2)
  
  # label N genes
  N = min(n.label, length(rownames(results[results$DE == 'FDRup',])))
  up_label = rownames(results[results$DE == 'FDRup',])[1:N]
  N = min(n.label, length(rownames(results[results$DE == 'FDRdown',])))
  down_label = rownames(results[results$DE == 'FDRdown',])[1:N]
  
  MAplot = ggplot(results) +
    geom_point(aes(x=baseMean, y=log2FoldChange, color = DE, alpha = DE2), size = 3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=baseMean, y=log2FoldChange, color = DE, alpha = DE2), size = 3) +
    xlim(c(0,1.e5)) +
    scale_x_continuous(trans='log10') +
    ggtitle(paste("MAPlot,", condition)) +
    scale_color_manual(values = mycolors) +
    #scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +  
    theme(plot.title = element_text(color="black", size=16, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "Mean expression", y = "log2 fold change")  
  
  #print(MAplot)
  #pdf(paste(f,'/','MAplot_',condition,'.pdf',sep=''),width=8, height=6.5)
  #plot(MAplot)
  #dev.off()
  
  # Vulcano plot 
  VP = ggplot(results) +
    geom_point(aes(x=log2FoldChange, y=-log10(padj), color = DE), size =2) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(padj), color = DE), size =2) +
    ggtitle(paste("Vulcano Plot,", condition)) +
    scale_color_manual(values = mycolors) +
    scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    geom_label_repel(data= results[c(up_label, down_label),], aes(x = log2FoldChange, y = -log10(padj), color = DE), 
                     label = row.names(results[c(up_label, down_label),]), size = 3,
                     box.padding = unit(0.35, "lines"), point.padding = unit(0.5, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=16, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 adjusted p-value")   
  
  print(VP)
  pdf(paste(f,'/','VulcanoPlot_',condition,'.pdf',sep=''),width=8, height=6.5)
  plot(VP)
  dev.off()
  
  options(repr.plot.width=14, repr.plot.height=6.5)
  p1 = MAplot ; p2 = VP;   
  print((p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
          (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt"))) +  
          plot_layout(guides = "collect"))
  
  pdf(paste(f,'/','MA_VP_',condition,'.pdf',sep=''),width=16, height=6.5)
  print((p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
          (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt"))) +  
          plot_layout(guides = "collect"))
  dev.off()
  
  A1 <- image_read_pdf(paste(f,'/','MA_VP_',condition,'.pdf',sep=''), density = 70)
  image_write(A1, path = paste(f,'/','MA_VP_',condition,'.tiff',sep=''), format = "tiff")
  
}

#------------------ FDR ------------------
# save FDR genes in a list
fdrUP = list()

alpha = 0.05
fdrUP = lapply(names(dgeResults), 
               function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$padj <= alpha & 
                                                        !is.na(dgeResults[[x]]$padj)&
                                                        dgeResults[[x]]$log2FoldChange > 0])
names(fdrUP)= names(dgeResults)       

fdrDW = list()
fdrDW = lapply(names(dgeResults), 
               function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$padj <= alpha & 
                                                        !is.na(dgeResults[[x]]$padj)&
                                                        dgeResults[[x]]$log2FoldChange < 0])
names(fdrDW)= names(dgeResults)

#------------------ SEQC ------------------
# save SEQC genes in a list
seqcUP = list()

pvalue = 0.01
seqcUP = lapply(names(dgeResults), 
                function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$pvalue <= pvalue & 
                                                         !is.na(dgeResults[[x]]$padj)&
                                                         dgeResults[[x]]$log2FoldChange > 1])
names(seqcUP)= names(dgeResults)               

seqcDW = list()
seqcDW = lapply(names(dgeResults), 
                function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$pvalue <= pvalue & 
                                                         !is.na(dgeResults[[x]]$padj)&
                                                         dgeResults[[x]]$log2FoldChange < -1])
names(seqcDW)= names(dgeResults)

print('FDRup');lengths(fdrUP); print('FDRdw');lengths(fdrDW)
print('SEQCup');lengths(seqcUP); print('SEQCdw');lengths(seqcDW)

print('FDR_tot');lengths(fdrUP) + lengths(fdrDW)
print('SEQC_tot'); lengths(seqcUP) + lengths(seqcDW)

#------------------ Heatmaps with FDR genes AS- vertical orientation ------------------
for (i in names(fdrUP)) {
  re = c(fdrUP[[i]],fdrDW[[i]])
  #a = unlist(str_split(i,'_'))
  #c1=a[2]; c2=a[4]
  #samples= as.character(metadata[metadata$condition %in% c(c1,c2),]$SampleID)
  HP <- pheatmap::pheatmap(fCountsRPKM[re,],
                           scale = 'row',
                           annotation_col = annotation_column,
                           annotation_colors = ann_colors, 
                           cluster_rows = T, 
                           cluster_cols = T, 
                           cutree_cols  = 3,
                           cellwidth = 25, cellheight = 0.03,
                           cutree_rows  = 2,
                           main = paste('FDR', i),
                           show_rownames = F,
                           fontsize = 12, fontsize_row = 10, fontsize_col = 14, 
                           display_numbers = F, 
                           col=colors,
                           filename = paste(f,'/allSample_Heatmap_v_',i,'.pdf',sep=''))
  print(HP)
}

############ Heatmaps with FDR genes SS - horizontal orientation ############
for (i in names(fdrUP)) {
  re = c(fdrUP[[i]],fdrDW[[i]])
  #a = unlist(str_split(i,'_'))
  #c1=a[2]; c2=a[4]
  #samples= as.character(metadata[metadata$condition %in% c(c1,c2),]$SampleID)
  HP <- pheatmap::pheatmap(t(fCountsRPKM[re,]),
                           scale = 'column',
                           annotation_row = annotation_column,
                           annotation_colors = ann_colors, 
                           cluster_rows = T, 
                           cluster_cols = T, 
                           cellwidth = 0.03, cellheight = 25,
                           cutree_rows = 2,
                           cutree_cols = 2,
                           main = paste('FDR', i),
                           show_rownames = T,
                           show_colnames = F,
                           fontsize = 12, fontsize_row = 9, fontsize_col = 14, 
                           display_numbers = F, 
                           col=colors,
                           filename =  paste(f,'/allSample_Heatmap_h_',i,'.pdf',sep=''))
  print(HP)
}

#------------------ Heatmaps with FDR genes SS - vertical orientation ------------------
for (i in names(fdrUP)) {
  re = c(fdrUP[[i]],fdrDW[[i]])
  a = unlist(str_split(i,'_'))
  c1=a[2]; c2=a[4]
  samples= as.character(metadata[metadata$condition %in% c(c1,c2),]$SampleID)
  HP <- pheatmap::pheatmap(fCountsRPKM[re,samples],
                           scale = 'row',
                           annotation_col = annotation_column,
                           annotation_colors = ann_colors, 
                           cluster_rows = T, 
                           cluster_cols = T, 
                           cutree_cols  = 2,
                           cellwidth = 25, cellheight = 0.03,
                           cutree_rows  = 2,
                           main = paste('FDR', i),
                           show_rownames = F,
                           fontsize = 12, fontsize_row = 10, fontsize_col = 14, 
                           display_numbers = F, 
                           col=colors,
                           filename = paste(f,'/Heatmap_v_',i,'.pdf',sep=''))
  print(HP)
}

#------------------ Heatmaps with FDR genes SS - horizontal orientation ------------------
for (i in names(fdrUP)) {
  re = c(fdrUP[[i]],fdrDW[[i]])
  a = unlist(str_split(i,'_'))
  c1=a[2]; c2=a[4]
  samples= as.character(metadata[metadata$condition %in% c(c1,c2),]$SampleID)
  HP <- pheatmap::pheatmap(t(fCountsRPKM[re,samples]),
                           scale = 'column',
                           annotation_row = annotation_column,
                           annotation_colors = ann_colors, 
                           cluster_rows = T, 
                           cluster_cols = T, 
                           cellwidth = 0.03, cellheight = 25,
                           cutree_rows = 2,
                           cutree_cols = 2,
                           main = paste('FDR', i),
                           show_rownames = T,
                           show_colnames = F,
                           fontsize = 12, fontsize_row = 9, fontsize_col = 14, 
                           display_numbers = F, 
                           col=colors,
                           filename =  paste(f,'/Heatmap_h_',i,'.pdf',sep=''))
  print(HP)
}

################ FDR gene intersection 3d vs 2d - venn::venn ################
list_Venn = fdrUP
names(list_Venn) = str_sub(names(list_Venn),11,-1)
a = venn(list_Venn,
         simplify=F, opacity = 0.3, box = FALSE, ilab=TRUE, plotsize = 200, 
         zcolor = rainbow(2), intersection = T, ilcs = 1.5, sncs = 1.5)
list_intersection = attr(x = a, "intersections")
names(list_intersection) <- str_replace_all(names(list_intersection), ':', '_AND_')
list_dataframe = list()

for (element in names(list_intersection)) {
  if (nchar(element) > 12) {
    entry = data.frame(Intersection = element, gene = list_intersection[[element]])
    list_dataframe[[element]] = entry
  }
}

write.xlsx(list_dataframe,
           file = 'VENNintersection_up.xlsx',  
           row.names = F,
           asTable = T,
           sheetName = str_sub(names(list_dataframe), 1,31))

list_Venn = fdrDW
names(list_Venn) = str_sub(names(list_Venn),11,-1)
a = venn(list_Venn,
         simplify=F, opacity = 0.3, box = FALSE, ilab=TRUE, plotsize = 200, 
         zcolor = rainbow(2), intersection = T, ilcs = 1.5, sncs = 1.5)
list_intersection = attr(x = a, "intersections")
names(list_intersection) <- str_replace_all(names(list_intersection), ':', '_AND_')
list_dataframe = list()
for (element in names(list_intersection)) {
  if (nchar(element) > 12) {
    entry = data.frame(Intersection = element, gene = list_intersection[[element]])
    list_dataframe[[element]] = entry
  }
}

write.xlsx(list_dataframe,
           file = 'VENNintersection_dw.xlsx',  
           row.names = F,
           asTable = T,
           sheetName = str_sub(names(list_dataframe), 1,31)) 

#################### FDR gene intersection 3d vs 2d - VennPlot ################
# intersections
colors = mycolors_c

# Venn plots
list_Venn = fdrUP
venn.plot <- venn.diagram(
  x = list_Venn,
  category.names = str_sub(names(fdrDW),11,-1),
  main = 'UP genes',
  filename = "Venn_genes_UP_3Dvs2D.tiff",
  scaled = TRUE,
  ext.text = TRUE,
  fill = mycolors_c,
  alpha = 0.6,
  cex = 2,
  cat.cex = 1,
  cat.pos = c(-20,20,180),
  cat.col = mycolors_c,
  main.cex = 2,
  sub.cex = 1)

list_Venn = fdrDW
venn.plot <- venn.diagram(
  x = list_Venn,
  category.names = str_sub(names(fdrDW),11,-1),
  main = 'DW genes',
  filename = "Venn_genes_DW_3Dvs2D.tiff",
  scaled = TRUE,
  ext.text = TRUE,
  fill = mycolors_c,
  alpha = 0.6,
  cex = 2,
  cat.cex = 1,
  cat.pos = c(-20,20,180),
  cat.col = mycolors_c,
  main.cex = 2,
  sub.cex = 1)

#BottleRocket2
Gstar = image_read("Venn_genes_UP_3Dvs2D.tiff")
Pstar = image_read("Venn_genes_DW_3Dvs2D.tiff")
im_star = image_append(c(Gstar, Pstar))
im_star
image_write(im_star, path = "Venn_UP_and_DW_3Dvs2D.tiff", format = "tiff")


#----------------  enrichment ----------------  

#' GSEA
#' Prepare lists of genes to run GSEA
outdir = 'GSEA/'
dir.create(outdir)
for (c in names(dgeResults)) {
  l_ranked = data.frame(GeneID= row.names(dgeResults[[c]]), LogFC = dgeResults[[c]][,'log2FoldChange'])
  l_ranked = l_ranked[order(l_ranked$LogFC, decreasing = T),]
  write.table(l_ranked, file = paste(outdir, c,'_ranked_list.rnk', sep =''), 
              quote = F, row.names= F, col.names = F, sep ='\t')
}

#' EnrichR
databases <- listEnrichrDbs()
dir.create('enrichR/', showWarnings=FALSE, recursive=TRUE)

#' enrichment Parameters
#' databases to make the enrichment of
enrich.databases <- c("GO_Biological_Process_2018",
                      "GO_Cellular_Component_2018",
                      "GO_Molecular_Function_2018",
                      "Reactome_2016",
                      "KEGG_2016",
                      "WikiPathways_2016",
                      "BioCarta_2016")
#' alpha used in DGE
padj.cutoff = alpha; 

#' Perform Enrichment
enrichr.list <- list()
for (i in 1:length(dgeResults)){
  print(names(dgeResults)[i])
  .res <- dgeResults[[i]]
  up.genes   <- fdrUP[[i]]
  down.genes <- fdrDW[[i]]
  both.genes <- c(up.genes, down.genes)
  write.table(up.genes, paste('./enrichR/FDRup_',names(dgeResults)[i],
                              '.txt', sep =''), quote = F, 
              row.names = F, col.names = F)
  write.table(down.genes, paste('./enrichR/FDRdw_',names(dgeResults)[i],
                                '.txt', sep =''), 
              quote = F, row.names = F, col.names = F)
  write.table(both.genes, paste('./enrichR/FDRboth_',names(dgeResults)[i],
                                '.txt', sep =''), quote = F, 
              row.names = F, col.names = F)
  
  
  enrichr.list[[i]] <- lapply(list(up.genes,down.genes,both.genes),function(x) {
    enrichR::enrichr(genes = x, databases = enrich.databases)
  })
  names(enrichr.list[[i]]) <-  c("fdr_up","fdr_down","fdr_both")
}
names(enrichr.list) <- names(dgeResults)

#' Write excels files
for (i in 1:length(dgeResults)){
  for (j in c("fdr_up","fdr_down","fdr_both")){
    filename = paste(
      file.path('enrichR', 
                names(dgeResults)[[i]]),
      j,
      ".xlsx",
      sep="_")
    write.xlsx(x = enrichr.list[[i]][[j]], file = filename)
  }
}


#----------- Jaccard distances -----------
#' Jaccard distances
dir = paste("./enrichR/JaccardPlots/", sep = '')
dir.create(dir)
#'WikiPathways_2016'
#database = c('KEGG_2016','Reactome_2016','GO_Biological_Process_2018')
p.value.thr = 0.05

breaksList = seq(0, 1, by = 0.001)

cutree_rows_values = c(5,12,3,10,5,1)

lf = list.files(paste("./enrichR/", sep=''), pattern=glob2rx("*.xlsx"))
for (file in lf) {    
  enrichR.file = paste("./enrichR/",file,sep='')
  s = openxlsx::getSheetNames(enrichR.file)
  Pathways.Table = data.frame()
  for (dat in s[1:length(s)]) {
    Table <- read.xlsx(xlsxFile = enrichR.file, 
                       sheet = dat, 
                       startRow = 1, 
                       colNames = TRUE,
                       rowNames = TRUE, 
                       detectDates = FALSE, 
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA", 
                       fillMergedCells = FALSE)
    
    Pathways.Table = rbind(Pathways.Table, Table)
  }
  
  gene.list = list()
  gene.all = character()
  pathways = unlist(row.names(Pathways.Table[Pathways.Table$Adjusted.P.value < p.value.thr,]))
  
  for (p in pathways) {
    gene.list[[p]] <- unlist(strsplit(Pathways.Table[p,]$Genes, ';')) 
    gene.all = c(gene.all, unlist(strsplit(Pathways.Table[p,]$Genes, ';')))
  }
  gene.all = unique(gene.all)
  
  if (length(gene.all) !=0) {
    # MAtrix
    M = matrix(0, nrow = length(pathways), ncol = length(gene.all))
    row.names(M) = pathways 
    colnames(M) = gene.all 
    
    for (pat in pathways) {
      for (gene in gene.all) {
        if (gene %in% gene.list[[pat]]) {
          M[pat,gene] <- 1 
        }            
      }    
    }
    
    if (length(pathways) >1) {
      # Jaccard dist
      Jacard.Matrix <- distance(M, method = "jaccard")
      if (length(pathways)==2) {
        Jacard.Matrix_new = as.matrix(rbind(c(0,Jacard.Matrix),c(Jacard.Matrix,0)))
        Jacard.Matrix = Jacard.Matrix_new
      }
      
      row.names(Jacard.Matrix) <- pathways
      colnames(Jacard.Matrix) <- pathways
      
      w=5; h=5; fs = 4; cutree_rows_N = 5
      if (file =="condition_PTs3Dday7_vs_PTsbasale_fdr_down_.xlsx") {cutree_rows_N = 2; w=25; h=25; fs= 22}
      myb = seq(0,1,by = 0.01)
      myc = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(length(myb))
      pheatmap(Jacard.Matrix,
               border_color = 'darkgrey',
               color = myc, 
               breaks = myb,
               cluster_rows = TRUE,
               cluster_cols = TRUE, 
               cellwidth = w, cellheight = h,
               cutree_rows = cutree_rows_N,
               show_colnames = FALSE,
               main = paste(file,'- Jaccard distance heatmap'),
               fontsize = 12,
               fontsize_row = fs,
               filename = paste(dir,file,'_JaccardDist.pdf', sep=''))
    }
  }
}

dev.off()

############## Plot top N pathways in a plot ################
fx <- function(x) eval(parse(text=enrichR.table[x,]$Overlap))
N=25
plotdir = paste('./enrichR/','TOP_',N,'_pathways/',sep='')
dir.create(plotdir, recursive = T)

path.sigs = list()

lf = list.files(paste("./enrichR/", sep=''), pattern=glob2rx("*.xlsx"))
for (file in lf) { 
  pathways.dataframe = data.frame()
  enrichR.file = paste("./enrichR/",file,sep='')
  s = openxlsx::getSheetNames(enrichR.file)
  enrichR.table = data.frame()
  for (dat in s[1:length(s)]) {
    Table <- read.xlsx(xlsxFile = enrichR.file, 
                       sheet = dat, 
                       startRow = 1, 
                       colNames = TRUE,
                       rowNames = TRUE, 
                       detectDates = FALSE, 
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA", 
                       fillMergedCells = FALSE)
    
    enrichR.table = rbind(enrichR.table, Table)
  }
  p = row.names(enrichR.table[enrichR.table$Adjusted.P.value < 0.005,])
  pathways.dataframe.entry = data.frame(Pathway = p,
                                        gene.ratio = sapply(p, fx),
                                        p.value = enrichR.table[p,]$P.value,
                                        p.value.adj = enrichR.table[p,]$Adjusted.P.value)
  pathways.dataframe <- rbind(pathways.dataframe, pathways.dataframe.entry) 
  
  path.sigs[[file]] = pathways.dataframe
  
  #up
  
  #pathways.dataframe = path.sigs[[file]]
  pathways.dataframe = pathways.dataframe[order(pathways.dataframe$p.value.adj),]
  pathways.dataframe$Pathway.num = dim(pathways.dataframe)[1]:1
  pathways.dataframe$Pathway.num = as.factor(pathways.dataframe$Pathway.num )
  
  
  pd = ggplot(pathways.dataframe[1:N,], aes(Pathway.num,-log10(p.value.adj), 
                                            color = -log10(p.value))) + 
    geom_point(aes(size = gene.ratio)) +
    scale_size_continuous(range = c(2,8), name = "Gene ratio") +
    scale_color_continuous(type = "viridis") + #, limits =c(1,35)) + 
    #scale_colour_viridis_c(begin = 0, end = 40) + 
    #ylim(c(10,40)) +
    coord_flip() +
    scale_x_discrete(breaks=pathways.dataframe$Pathway.num, 
                     labels=stringr::str_wrap(pathways.dataframe$Pathway,width = 85),
                     position = "top") +
    theme(plot.title = element_text(color="black", size=18, face="bold.italic"),
          plot.subtitle = element_text(color="black", size=16, face="italic"),
          axis.text.x = element_text(angle = 90,, color = 'black', size=14, hjust =1, family = "mono"), 
          axis.title.x = element_text(face = "bold", color = "black", size = 14),
          axis.text.y = element_text(angle = 0, color = 'black', size=14, family = "mono"),
          axis.title.y = element_text(face = "bold", color = "black", size = 14),
          legend.text = element_text(color = "black", size = 12),
          legend.title = element_text(face = "bold", color = "black", size = 14),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(title = paste("Top", N, 'pathways -', str_sub(file,1,-7)), 
         x = "Pathways", 
         y = "-Log10(p.adj.value)") 
  print(pd)
  pdf(paste(plotdir, str_sub(file,1,-7),'_TOP',N,'.pdf',sep=''), 
      width = 14, height = 8)
  print(pd)
  dev.off()
  
}    



