suppressMessages(library(WGCNA))
options(stringsAsFactors=FALSE);allowWGCNAThreads(nThreads=4);set.seed(20261021)
O<-"/path/to/storage/kleiber_birds_gemcat/results/growth_regulators_2026-10-21_v1"
inv<-read.delim(file.path(O,"config/full_inventory.tsv"),check.names=FALSE)
dir.create(file.path(O,"results/wgcna"),recursive=TRUE,showWarnings=FALSE)
writeLines(capture.output(sessionInfo()),file.path(O,"audit/R_full_sessionInfo.txt"))
receipt<-list()
for(panel in inv$panel){
 m<-as.matrix(read.delim(file.path(O,"inputs",paste0(panel,".log.tsv")),row.names=1,check.names=FALSE))
 for(mode in c("raw","sample_centered")){
  out<-file.path(O,"results/wgcna",paste(panel,mode,sep="__"));dir.create(out,showWarnings=FALSE)
  x<-t(m);if(mode=="sample_centered")x<-x-rowMeans(x)
  good<-goodSamplesGenes(x,verbose=0);x<-x[good$goodSamples,good$goodGenes,drop=FALSE]
  cat("START",panel,mode,nrow(x),ncol(x),"\n");flush.console()
  if(ncol(x)>5000){v<-apply(x,2,mad);anchor<-scan(file.path(O,"config/forced_anchors.txt"),what="character",quiet=TRUE);chosen<-unique(c(colnames(x)[order(-v,names(v))[1:5000]],intersect(anchor,colnames(x))));x<-x[,chosen,drop=FALSE]}
  powers<-c(1:10,12,14,16,18,20)
  sft<-pickSoftThreshold(x,powerVector=powers,networkType="signed",verbose=0)
  write.table(sft$fitIndices,file.path(out,"soft_threshold.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  eligible<-with(sft$fitIndices,SFT.R.sq>=.8 & slope<0)
  fallback<-if(nrow(x)<20)18 else if(nrow(x)<30)16 else if(nrow(x)<40)14 else 12
  power<-if(any(eligible))sft$fitIndices$Power[which(eligible)[1]] else fallback
  net<-blockwiseModules(x,power=power,networkType="signed",TOMType="signed",corType="pearson",minModuleSize=30,deepSplit=2,mergeCutHeight=.25,numericLabels=TRUE,pamRespectsDendro=FALSE,maxBlockSize=5000,randomSeed=20261021,verbose=0,saveTOMs=FALSE)
  color<-labels2colors(net$colors);means<-colMeans(x);sds<-apply(x,2,sd);z<-scale(x);members<-list();scores<-list();weights<-list();i<-0
  for(module in unique(color)){
   ii<-which(color==module);xx<-z[,ii,drop=FALSE];sv<-svd(xx,nu=1,nv=1);w<-sv$v[,1];sc<-drop(xx%*%w)
   if(cor(sc,rowMeans(xx))<0){w<--w;sc<--sc}
   kme<-drop(cor(xx,sc));i<-i+1
   members[[i]]<-data.frame(feature=colnames(x)[ii],module=module,kME=kme,mean_log=means[ii],sd_log=sds[ii],weight=w,n_features=length(ii),variance_PC1=sv$d[1]^2/sum(sv$d^2))
   scores[[i]]<-data.frame(sample=rownames(x),module=module,eigengene=sc)
  }
  write.table(do.call(rbind,members),file.path(out,"members.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  write.table(do.call(rbind,scores),file.path(out,"eigengenes.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  write.table(data.frame(panel=panel,mode=mode,power=power,n_samples=nrow(x),n_features=ncol(x),n_modules=length(setdiff(color,"grey")),scale_free_reached=any(eligible)),file.path(out,"run.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
  receipt[[length(receipt)+1]]<-data.frame(panel=panel,mode=mode,power=power,n_samples=nrow(x),n_features=ncol(x),n_modules=length(setdiff(color,"grey")),grey=sum(color=="grey"))
  cat("DONE",panel,mode,"modules",length(unique(color))-as.integer("grey"%in%color),"\n");flush.console();rm(net,z);gc()
 }
}
write.table(do.call(rbind,receipt),file.path(O,"results/wgcna_full_runs.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
