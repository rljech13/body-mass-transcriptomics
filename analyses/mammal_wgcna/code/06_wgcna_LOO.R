suppressMessages(library(WGCNA));options(stringsAsFactors=FALSE);allowWGCNAThreads(nThreads=4);set.seed(20261021)
O<-"/path/to/storage/kleiber_birds_gemcat/results/growth_regulators_2026-10-21_v1"
rows<-list();j<-0
for(panel in c("ma2016","bird_HG","bird_F12")){
 x<-t(as.matrix(read.delim(file.path(O,"inputs",paste0(panel,".log.tsv")),row.names=1,check.names=FALSE)))
 meta<-read.delim(file.path(O,"inputs",paste0(panel,".samples.tsv")),check.names=FALSE);names(meta)[1]<-"sample"
 full<-read.delim(file.path(O,"results/wgcna",paste0(panel,"__raw"),"members.tsv"))
 run<-read.delim(file.path(O,"results/wgcna",paste0(panel,"__raw"),"run.tsv"))
 for(sp in unique(meta$species)){
  xs<-x[meta$species[match(rownames(x),meta$sample)]!=sp,,drop=FALSE]
  good<-goodSamplesGenes(xs,verbose=0);xs<-xs[,good$goodGenes,drop=FALSE]
  net<-blockwiseModules(xs,power=run$power,networkType="signed",TOMType="signed",corType="pearson",minModuleSize=30,deepSplit=2,mergeCutHeight=.25,numericLabels=TRUE,pamRespectsDendro=FALSE,maxBlockSize=5000,randomSeed=20261021,verbose=0,saveTOMs=FALSE)
  cc<-labels2colors(net$colors);names(cc)<-colnames(xs)
  for(mod in setdiff(unique(full$module),"grey")){
   a<-full$feature[full$module==mod];best<-0;matchmod<-"none";intersection<-0;bestsize<-0
   for(mm in setdiff(unique(cc),"grey")){
    b<-names(cc)[cc==mm];k<-length(intersect(a,b));jac<-k/length(union(a,b))
    if(jac>best){best<-jac;matchmod<-mm;intersection<-k;bestsize<-length(b)}
   }
   j<-j+1;rows[[j]]<-data.frame(panel=panel,mode="raw",omitted_species=sp,module=mod,n_full=length(a),n_LOO=bestsize,best_module=matchmod,intersection=intersection,Jaccard=best)
  }
  cat("LOO DONE",panel,sp,"\n");flush.console();rm(net);gc()
  write.table(do.call(rbind,rows),file.path(O,"results/wgcna_LOO_membership.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
 }
}
