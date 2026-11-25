suppressPackageStartupMessages(library(jsonlite))
a<-commandArgs(trailingOnly=TRUE);w<-normalizePath(a[1]);rd<-function(n)read.delim(file.path(w,'data',n),check.names=FALSE,stringsAsFactors=FALSE)
e<-rd('saved_eigengenes.tsv');members<-rd('saved_members.tsv');species<-rd('species.tsv');old<-rd('saved_historical_mass_tests.tsv');hm<-rd('saved_Hallmark_enrichment.tsv')
rows<-list();scores<-list();i<-0
for(panel in unique(e$panel)){
 s<-species[species$panel==sub('_full$','',panel),];stopifnot(!anyDuplicated(s$species))
 for(module in sort(setdiff(unique(e$module[e$panel==panel]),'grey'))){
  d<-e[e$panel==panel&e$module==module,];k<-match(d$sample,s$species);stopifnot(!anyNA(k),nrow(d)==nrow(s),!anyDuplicated(d$sample))
  x<-log10(s$mass_g[k]);y<-d$eigengene;ct<-cor.test(x,y,method='pearson');md<-members[members$panel==panel&members$module==module,]
  h<-old[old$panel==panel&old$module==module,];stopifnot(nrow(h)==1,abs(ct$p.value-h$ols_p)<1e-8)
  enrichment<-hm[hm$panel==panel&hm$module==module,];enrichment<-enrichment[order(enrichment$q_all_tests,enrichment$p),]
  label<-if(nrow(enrichment)&&enrichment$q_all_tests[1]<.05)enrichment$pathway[1] else ''
  loo<-sapply(seq_along(x),function(j)cor(x[-j],y[-j]));i<-i+1
  rows[[i]]<-data.frame(panel=panel,module=module,n_species=nrow(d),n_genes=nrow(md),variance_PC1=md$variance_PC1[1],r=unname(ct$estimate),p=ct$p.value,
   historical_OLS_p=h$ols_p,historical_PGLS_p=h$p,historical_PGLS_q=h$q_joint_current_run,
   top_Hallmark=label,top_Hallmark_q=if(nrow(enrichment))enrichment$q_all_tests[1] else NA,
   loo_sign_fraction=mean(sign(loo)==sign(ct$estimate)),loo_min_abs_r=min(abs(loo)))
  scores[[i]]<-data.frame(panel=panel,module=module,species=d$sample,mass_g=s$mass_g[k],log10_mass=x,eigengene=y,z=as.numeric(scale(y)))
 }
}
t<-do.call(rbind,rows);t$q_within_study<-ave(t$p,t$panel,FUN=function(z)p.adjust(z,method='BH'));t$q_all3<-p.adjust(t$p,method='BH')
stopifnot(nrow(t)==123)
write.table(t,file.path(w,'data/module_mass_Pearson.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,scores),file.path(w,'data/module_relative_expression.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
sel<-t[t$panel=='ma2016_full'&t$q_within_study<.05,];sel<-sel[order(sel$r),]
write.table(sel,file.path(w,'data/displayed_modules.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
summary<-list(studies=lapply(split(t,t$panel),function(z)list(modules=nrow(z),BH05=sum(z$q_within_study<.05),BH_all3_05=sum(z$q_all3<.05),min_q=min(z$q_within_study))),selected_modules=sel,all_historical_OLS_p_reproduced=TRUE,network_rebuilt=FALSE,
 previous_historical_PGLS_mammal_tests_BH05=0,interpretation='New descriptive Pearson/BH presentation of the previously frozen modules; no new module discovery. Relative module expression, not absolute activation or flux.',human_review='not_completed',training_enabled=FALSE)
write_json(summary,file.path(w,'audit/summary.json'),auto_unbox=TRUE,pretty=TRUE,digits=12)
cat(toJSON(summary,auto_unbox=TRUE,pretty=TRUE,digits=6))
