
# Function for calculate differential expression
# matrix: matrix with samples in columns and genes in rows with symbol name
# metadata dataframe with metadata

differential.expression.voom <- function(matrix, metadata, design) {
  
  # Load matrix and normalization
  
  d0 <- DGEList(matrix)
 
  # Removing low cont genes
  
  keep <- filterByExpr(d0, design)
  d <- d0[keep,,keep.lib.sizes=FALSE]
  
  # Scale normalization
  d <- calcNormFactors(d)
  
  # Model matrix
  mm <- design

  # Voom
  y <- voom(d, mm)
  # Fittind the linear model
  fit <- lmFit(y, mm)
  return(fit)
}



differential.expression.trend<- function(matrix, metadata, design) {

  # Load matrix 
  d0 <- DGEList(matrix)
  
  # Removing low cont genes
  
  keep <- filterByExpr(d0, design)
  d <- d0[keep,,keep.lib.sizes=FALSE]
  
  # Scale normalization
  d <- calcNormFactors(d)
  
  # Convert to logCPM
  
  logCPM <- cpm(d, log=TRUE, prior.count=3)
  
  # Removing low cont genes
  

  fit <- lmFit(logCPM, design)
  fit <- eBayes(fit, trend=TRUE)
  
  return(fit)
}



Enrichment_Ananalysis_KEEG <- function(top.table){
  
  top.table <- top.table.DEA.Night.VS.Day.WT.voom
 
  
  # New annotation

  mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  
  IDs <- biomaRt::getBM(attributes = c("ensembl_gene_id","mgi_symbol"), 
                        values = row.names(top.table), mart = mart, useCache = F)
  
  colnames(IDs)<-c("sample", "mgi_symbol")
  
  top.table$mgi_symbol <- row.names(top.table)
  #top.table[,"mgi_symbol"] = toupper( top.table[,"mgi_symbol"])
  top.table <- merge(top.table, IDs, by="mgi_symbol")
  
  
  top.table <- top.table[!duplicated(top.table$sample), ]
  
  row.names(top.table) <- top.table$sample
  
  mask <- top.table$adj.P.Val < 0.05 
  
  deGenes <- rownames(top.table[mask, ])
  
  length(deGenes)
  
  geneUniverse <- rownames(top.table)
  length(geneUniverse)
  
  deGenes <- unlist(mget(deGenes, envir=org.Mm.egENSEMBL2EG,
                         ifnotfound = NA))
  
  geneUniverse <- unlist(mget(geneUniverse, envir=org.Mm.egENSEMBL2EG,
                              ifnotfound = NA))
  
  
  ## KEGG Enrichment
  ans.kegg <- enrichKEGG(gene = deGenes, 
                         organism = 'mmu', 
                         universe = geneUniverse, 
                         pvalueCutoff = 0.05,
                         use_internal_data=T,
                         pAdjustMethod = "none")
  
  
  
  geneList = top.table[mask, ][6]
  geneList$gene <- row.names(geneList)
  
  geneList <-   geneList[order(-geneList$adj.P.Val),]
  
  
  ids <- unlist(mget(row.names(geneList), envir=org.Mm.egENSEMBL2EG,
                         ifnotfound = NA))
  
  deGenes <- geneList$adj.P.Val
  names(deGenes) <- as.character(na.omit(ids))
  
  
  ans.kegg2 <- gseKEGG(gene = deGenes, 
                         organism = 'mouse', 
                         minGSSize    = 20,
                         pvalueCutoff = 0.05,
                         use_internal_data=T,
                         pAdjustMethod = "none")
  
  
  return(ans.kegg)
  
}