args <- commandArgs(trailingOnly=TRUE)
base <- args[1]
library(edgeR)
for(cohort in c('Discovery','Validation')) {
  counts <- as.matrix(read.delim(file.path(base,paste0(cohort,'_edger_counts.tsv')),row.names=1,check.names=FALSE))
  meta <- read.delim(file.path(base,paste0(cohort,'_edger_samples.tsv')),row.names=1,check.names=FALSE)
  stopifnot(identical(colnames(counts),rownames(meta)))
  meta$donor <- factor(meta$donor)
  ref <- if(cohort=='Discovery') 'PT Injured' else 'PT'
  case <- if(cohort=='Discovery') 'PT Inflammatory' else 'PTVCAM1'
  meta$state <- factor(meta$state,levels=c(ref,case))
  design <- model.matrix(~donor+state,data=meta)
  y <- DGEList(counts=counts,lib.size=meta$total_umi)
  # Preserve the all-feature UMI denominator even after filtering.
  nf <- calcNormFactors(y,method='TMM')$samples$norm.factors
  write.table(data.frame(sample=rownames(meta),norm_factor=nf,total_umi=meta$total_umi),file.path(base,paste0(cohort,'tmm_factors.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
  keep <- filterByExpr(y,design=design)
  write.table(data.frame(gene=rownames(y),retained=keep),file.path(base,paste0(cohort,'expression_filter.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
  y <- y[keep,,keep.lib.sizes=TRUE];y$samples$norm.factors <- nf
  y <- estimateDisp(y,design,robust=TRUE)
  fit <- glmQLFit(y,design,robust=TRUE)
  test <- glmQLFTest(fit,coef=ncol(design))
  tt <- topTags(test,n=Inf,sort.by='none')$table
  tt <- data.frame(gene=rownames(tt),tt,row.names=NULL)
  write.table(tt,file.path(base,paste0(cohort,'_edger_genomewide.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
  capture.output(list(cohort=cohort,reference=ref,case=case,n_samples=ncol(y),n_genes=nrow(y),design=design,common_dispersion=y$common.dispersion,session=sessionInfo()),file=file.path(base,paste0(cohort,'_edger_run.txt')))
  cat(cohort,nrow(y),'genes',sum(tt$FDR<.05),'FDR<0.05\n')
  print(tt[tt$gene %in% c('PPARGC1B','ESRRA','HNF4A','SDHAF1','MGME1'),])
}
