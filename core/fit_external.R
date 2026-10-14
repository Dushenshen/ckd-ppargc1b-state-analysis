library(edgeR)
b <- '.'
counts <- as.matrix(read.delim(file.path(b,'external_counts.tsv'),row.names=1,check.names=FALSE));counts <- rowsum(counts,group=rownames(counts),reorder=FALSE)
m <- read.delim(file.path(b,'external_samples.tsv'),row.names=1,check.names=FALSE)
eligible <- intersect(m$donor[m$state=='Adaptive PT' & m$n_cells>=20],m$donor[m$state=='Mature PT' & m$n_cells>=20])
m <- m[m$donor %in% eligible & m$state %in% c('Adaptive PT','Mature PT'),];counts<-counts[,rownames(m),drop=FALSE]
m$donor<-factor(m$donor);m$state<-factor(m$state,levels=c('Mature PT','Adaptive PT'));design<-model.matrix(~donor+state,m)
y<-DGEList(counts=counts,lib.size=m$total_umi);y<-calcNormFactors(y);keep<-filterByExpr(y,design=design)
write.table(data.frame(gene=rownames(y),retained=keep),file.path(b,'external_expression_filter.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
y<-y[keep,,keep.lib.sizes=TRUE];y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE);tt<-topTags(glmQLFTest(fit,coef=ncol(design)),n=Inf,sort.by='none')$table
write.table(data.frame(gene=rownames(tt),tt),file.path(b,'external_edger.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
capture.output(list(n_donors=length(eligible),n_genes=nrow(y),design=design,session=sessionInfo()),file=file.path(b,'external_edger_run.txt'))
print(tt[rownames(tt) %in% c('PPARGC1B','ESRRA','HNF4A','SDHAF1','MGME1'),])
