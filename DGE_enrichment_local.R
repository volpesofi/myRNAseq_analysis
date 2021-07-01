################ Prj set up ################ 
# set working diretory
prj = 'ScielzoC_1381_RNASeq' 
PI = 'Scielzio'
setwd(paste("/Users/tascini.annasofia/Dropbox (HSR Global)/WORKSPACE/",
            PI,"/",prj,
            "/7_bioinfo/",
            sep=''))

################ load libraries ################ 
suppressMessages(library("edgeR"))
suppressMessages(library(data.table))
suppressMessages(library("DESeq2"))
#suppressMessages(library(openxlsx))
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
suppressMessages(library('philentropy'))
suppressMessages(library('IntClust'))
suppressMessages(library('pheatmap'))
suppressMessages(library(assertr))
suppressMessages(library("remotes"))
suppressMessages(library(GeneOverlap))
suppressMessages(library('stringr'))
suppressMessages(library(ggrepel))
#BiocManager::install("magick",force=T)
suppressMessages(library(magick))
suppressMessages(library(pdftools))
suppressMessages(library(forcats))

library("msigdbr")
library("vroom")
library(enrichplot)
library(msigdbr)
#library(xlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(purrr)

################ DGE - Import results ######################## 

###import dgeResults list saved as Rdata
dgeResults <- list.files(pattern = "dgeResults") %>% map(readRDS)
dgeResults= unlist(dgeResults)


#  names(dgeResults),
# function(x) write.table(
#  data.table(
#   data.frame(dgeResults[[x]]),
#   keep.rownames=geneidColname),
# file.path(f, paste(x, ".tsv", sep="")),
# append=F,
# row.names=F,
# col.names=T,
#  quote=F,
# sep="\t"))


dgeResults_table = list()    
dgeResults_table = lapply(
  names(dgeResults),
  function(x) 
    data.table(
      data.frame(dgeResults[[x]]),
      keep.rownames=geneidColname))
names(dgeResults_table) = names(dgeResults)

library(openxlsx)

write.xlsx(dgeResults_table,
           file = paste('DGE_results.xlsx', sep=''), 
           row.names = F,
           asTable = T, 
           sheetName = str_sub(names(dgeResults),1,31)) 
detach("package:openxlsx", unload = TRUE)


names(dgeResults)

#######################creo una cartella per ogni confronto
for (f in names(dgeResults)){
  dir.create(f, showWarnings=TRUE, recursive=TRUE)
  setwd(f)
  file=as.data.frame(dgeResults[[f]])
  
  
  file$FDR = ifelse(file$padj <= alpha &!is.na(file$padj)&file$log2FoldChange > 0,1,(ifelse(file$padj <= alpha &!is.na(file$padj)&
                                                                                              file$log2FoldChange < 0,-1,0)))
  file$SEQC = ifelse(file$pvalue <= 0.01 &!is.na(file$pvalue)&file$log2FoldChange >1,1,(ifelse(file$pvalue <= 0.01 &!is.na(file$pvalue)&
                                                                                                 file$log2FoldChange < -1,-1,0)))
  write.table(file,file = paste(f,'DGE_results.tsv', sep='_'),sep="\t",col.names=NA,quote=F)
  setwd(dirname(getwd()))
}

####Import GSEA signatures for later

###Hallmarks
#m_t2g_H <- read.gmt("/Users/eugeniabezzecchi/beegfs/scratch/ric.cosr/ric.sitiagiovanni/SitiaG_1321_mRNASeq/dataset/DGE_1321/h.all.v7.2.symbols.gmt")
m_t2g_H <- read.gmt("./myRNAseq_analysis/h.all.v7.2.symbols.gmt")

###C5_GO
#m_t2g_C5=read.gmt("/Users/eugeniabezzecchi/beegfs/scratch/ric.cosr/ric.sitiagiovanni/SitiaG_1321_mRNASeq/dataset/DGE_1321/c5.go.v7.2.symbols.gmt")
m_t2g_C5=read.gmt("./myRNAseq_analysis/c5.go.v7.2.symbols.gmt")

###C2CP
#m_t2g_C2=read.gmt("/Users/eugeniabezzecchi/beegfs/scratch/ric.cosr/ric.sitiagiovanni/SitiaG_1321_mRNASeq/dataset/DGE_1321/c2.cp.v7.2.symbols.gmt")
m_t2g_C2=read.gmt("./myRNAseq_analysis/c2.cp.v7.2.symbols.gmt")


##########################################ciclo sulle cartelle
f=NULL

for (f in names(dgeResults)){
  setwd(paste("/Users/tascini.annasofia/Dropbox (HSR Global)/WORKSPACE/",
              PI,"/",prj,
              "/7_bioinfo/",
              sep=''))
  #f=names(dgeResults)[1]
  condition=f
  setwd(f)
  
  
  ################ DGE - MAplot and Vulcano Plots ######################## 
  n.label = 20
  FDR = T
  pvalue = 0.01
  files <- list.files(pattern = "DGE_results.tsv",full =F,recursive = F)
  results<-read.table(files[1],sep="\t",header = T,row.names=1)
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
  pdf(paste('VulcanoPlot_',condition,'.pdf',sep=''),width=8, height=6.5)
  plot(VP)
  dev.off()
  
  options(repr.plot.width=14, repr.plot.height=6.5)
  p1 = MAplot ; p2 = VP;   
  print((p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
          (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt"))) +  
          plot_layout(guides = "collect"))
  
  pdf(paste('MA_VP_',condition,'.pdf',sep=''),width=16, height=6.5)
  print((p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
          (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt"))) +  
          plot_layout(guides = "collect"))
  dev.off()
  
  A1 <- image_read_pdf(paste('MA_VP_',condition,'.pdf',sep=''), density = 70)
  image_write(A1, path = paste('MA_VP_',condition,'.tiff',sep=''), format = "tiff")
  
  ###nominal pval barplot
  p = ggplot(results, aes(x=pvalue)) + geom_histogram(binwidth=0.05, color="blue", fill="skyblue1") + theme_linedraw()
  p = p + labs(title="Nominal pvalue distribution", x="nominal pvalue", y="frequency")
  p = p + theme(axis.title=element_text(size=17), axis.text=element_text(size=15), title=element_text(size=17))
  
  grDevices::png("Nominal pvalue distribution.png", width=8, height=8, units="in", res=300)
  plot(p)
  dev.off()
  
  
  
  alpha = 0.05
  fdrUP = row.names(results)[results$padj <= alpha &!is.na(results$padj)&
                               results$log2FoldChange > 0]
  
  fdrDN = row.names(results)[results$padj <= alpha &!is.na(results$padj)&
                               results$log2FoldChange < 0]
  
  
  
  ####### SEQC ########
  
  
  pvalue = 0.01
  seqcUP = row.names(results)[results$pvalue <= pvalue &!is.na(results$padj)&
                                results$log2FoldChange > 1]
  
  seqcDN = row.names(results)[results$pvalue <= pvalue &!is.na(results$padj)&
                                results$log2FoldChange < -1]
  
  
  genes<- qpcR:::cbind.na(fdrUP,fdrDN,seqcUP,seqcDN)
  
  
  
  
  ###### ENRICHR ########
  databases <- listEnrichrDbs()
  
  # databases to make the enrichment of
  enrich.databases <- c("GO_Biological_Process_2018",
                        "GO_Cellular_Component_2018",
                        "GO_Molecular_Function_2018",
                        "Reactome_2016",
                        "KEGG_2016",
                        "WikiPathways_2016",
                        "BioCarta_2016")
  # alpha used in DGE
  padj.cutoff = alpha; 
  
  # -------------------------
  # Perform Enrichment
  # -------------------------
  
  
  enrichr.list<- lapply(list(fdrUP,fdrDN,seqcUP,seqcDN),function(x) {
    enrichR::enrichr(genes = x, databases = enrich.databases)
  })
  names(enrichr.list) <-  c("fdr_up","fdr_down","seqc_up","seqc_dn")
  
  
  # -----------------------------
  # Write excels files
  # -----------------------------
  dir = paste('./enrichR/',sep='')
  dir.create(dir, recursive = T)
  
  library("xlsx2dfs")
  
  for (j in c("fdr_up","fdr_down","seqc_up","seqc_dn")){
    filename =paste(dir,j,"_enrichR.xlsx",sep='')
    dfs2xlsx(enrichr.list[[j]],filename,rowNames = F)}
  
  detach("package:xlsx2dfs", unload = TRUE)
  
  ############## Plot top N pathways in a plot ################
  N=10
  plotdir = paste('./enrichR/','TOP_',N,'_pathways/',sep='')
  dir.create(plotdir, recursive = T)
  
  
  for (file in names(enrichr.list)){ 
    for (dat in 1:length(enrichr.list[[file]])) {
      y=enrichr.list[[file]][[dat]]
      DB=c("GO_BP","GO_CC","GO_MF","reactome","kegg","wiki","biocarta")
      db=DB[dat]
      y$logPadj<--log10(y$Adjusted.P.value)
      y$Overlap<-sapply(y$Overlap, function(x) eval(parse(text=x)))
      y<-y[c(1:N),]
      y=na.omit(y)
      
      if (nrow(y)==0) 
      {
        next 
      }
      
      #col=cut(y$logPadj,quantile(y$logPadj))
      pdf(paste(plotdir,file,"_",db,'_TOP',N,'.pdf',sep=''),15,10)
      p=ggplot(y, # you can replace the numbers to the row number of pathway of your interest
               aes(x = logPadj, y = fct_reorder(Term, logPadj,.desc=F))) + 
        geom_point(aes(size = Overlap, color = logPadj)) +
        theme(plot.title = element_text(color="black",hjust =0, vjust =1, size=18, face="bold.italic"),
              plot.subtitle = element_text(color="black", size=16, face="italic"),
              axis.text.x = element_text(angle = 90, color = 'black', size=14, hjust =1, family = "mono"), 
              axis.title.x = element_text(face = "bold", color = "black", size = 14),
              axis.text.y = element_text(angle = 0, color = 'black', size=12, family = "mono"),
              axis.title.y = element_text(face = "bold", color = "black", size = 14),
              legend.text = element_text(color = "black", size = 12),
              legend.title = element_text(face = "bold", color = "black", size = 14),
              panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
        labs(title = paste(file,"-",db), 
             x = "-Log10(p.adj.value)", 
             y = "Pathways") 
      plot(p)
      dev.off()
      
    } 
  }
  
  
  outdir = 'GSEA/'
  dir.create(outdir)
  setwd(outdir)
  
  results<-na.omit(results)
  results$gene<-row.names(results)
  results$gene<-toupper(results$gene)
  #results$entrez<-data.frame(mapIds(org.Hs.eg.db, results$gene,"ENTREZID","SYMBOL"))
  
  d=data.frame(results$gene,results$log2FoldChange)
  
  ## feature 1: numeric vector
  geneList = d[,2]
  ## feature 2: named vector
  names(geneList) = as.character(d[,1])
  ## feature 3: decreasing order
  geneList = sort(geneList, decreasing = TRUE)
  geneList=na.omit(geneList)
  gene=geneList
  
  ####C2CP
  
  library(xlsx)
  
  
  z<- GSEA(gene, TERM2GENE=m_t2g_C2, pvalueCutoff = 0.1)
  #z<- setReadable(z, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  Curated_GSEA=as.data.frame(z)
  
  
  if (nrow(Curated_GSEA)!=0){
    
    write.xlsx(Curated_GSEA,"GSEA_results.xlsx", sheetName="Curated_pathways", 
               col.names=TRUE, row.names=TRUE, append=FALSE)
    
    # Make plots.
    plot_list = list()
    
    dim=ifelse(nrow(Curated_GSEA)<150,nrow(Curated_GSEA),150)
    for (j in 1:dim) {
      p = gseaplot2(z, geneSetID =j, title = z$Description[j])
      plot_list[[j]] = p
    }
    
    # create pdf where each page is a separate plot.
    pdf("GSEA_Curated_pathways.pdf")
    for  (j in 1:dim) {
      print(plot_list[[j]])
    }
    dev.off()
    
  }else{
    no_enrichment=data.frame(a="no_enrichment")
    write.xlsx(no_enrichment,"GSEA_results.xlsx", sheetName="Curated_pathways", 
               col.names=F, row.names=F, append=T)
  }
  
  ####Hallmarks
  
  z<- GSEA(gene, TERM2GENE=m_t2g_H, pvalueCutoff = 0.1)
  Hallmark=as.data.frame(z)
  
  if (nrow(Hallmark)!=0){
    
    write.xlsx(Hallmark,"GSEA_results.xlsx", sheetName="Hallmarks", 
               col.names=T, row.names=T, append=T)
    
    # Make plots.
    plot_list = list()
    dim=ifelse(nrow(Hallmark)<50,nrow(Hallmark),50)
    for (j in 1:dim) {
      p = gseaplot2(z, geneSetID =j, title = z$Description[j])
      plot_list[[j]] = p
    }
    
    pdf("GSEA_Hallmarks.pdf")
    for  (j in 1:dim) {
      print(plot_list[[j]])
    }
    dev.off()
    
  }else{
    no_enrichment=data.frame(a="no_enrichment")
    write.xlsx(no_enrichment,"GSEA_results.xlsx", sheetName="Hallmarks", 
               col.names=F, row.names=F, append=T)
  }
  
  
  z<- GSEA(gene, TERM2GENE=m_t2g_C5, pvalueCutoff = 0.1)
  GO=as.data.frame(z)
  
  if (nrow(GO)!=0){
    
    write.xlsx(GO,"GSEA_results.xlsx", sheetName="GO_all", 
               col.names=TRUE, row.names=TRUE, append=T)
    
    
    
    # Make plots.
    plot_list = list()
    dim=ifelse(nrow(GO)<150,nrow(GO),150)
    for (j in 1:dim) {
      p = gseaplot2(z, geneSetID =j, title = z$Description[j])
      plot_list[[j]] = p
    }
    
    pdf("GSEA_GO_all.pdf",8,6)
    for  (j in 1:dim) {
      print(plot_list[[j]])
    }
    dev.off()
    
  }else{
    no_enrichment=data.frame(a="no_enrichment")
    write.xlsx(no_enrichment,"GSEA_results.xlsx", sheetName="GO_all", 
               col.names=F, row.names=F, append=T)
  }
  
}

#save.image(file = paste(prj,".RData", sep=''))




