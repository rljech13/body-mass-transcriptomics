suppressPackageStartupMessages(library(jsonlite))
a<-commandArgs(trailingOnly=TRUE);w<-normalizePath(a[1]);r<-file.path(w,'results/analysis');out<-list()
for(seed in c(20261023,20261024,20261025)){
 d<-read.delim(file.path(r,paste0('overlap_',seed,'.tsv')));p<-phyper(d$k_shared-1,d$K_targets_background,d$N_background-d$K_targets_background,d$n_module_covered,lower.tail=FALSE);q<-p.adjust(p,'BH');stopifnot(max(abs(p-d$p))<1e-10,max(abs(q-d$q))<1e-10)
 out[[as.character(seed)]]<-list(tests=nrow(d),p_max_error=max(abs(p-d$p)),q_max_error=max(abs(q-d$q)),hits=sum(q<.05))
}
x<-read.delim(file.path(r,'AUCell_species.tsv'),row.names=1,check.names=FALSE);traits<-read.delim(file.path(w,'inputs/species_traits.tsv'),row.names=1,check.names=FALSE);mass<-log10(traits[rownames(x),'mass_g']);tests<-read.delim(file.path(r,'regulon_mass_associations.tsv'),check.names=FALSE);ps<-sapply(tests$regulon,function(n)if(diff(range(x[[n]]))<=1e-12)1 else cor.test(x[[n]],mass,method='pearson')$p.value);stopifnot(max(abs(ps-tests$mass_p))<1e-9,max(abs(p.adjust(ps,'BH')-tests$mass_q))<1e-9)
write_json(list(hypergeometric_BH=out,mass_p_max_error=max(abs(ps-tests$mass_p)),mass_q_max_error=max(abs(p.adjust(ps,'BH')-tests$mass_q)),R=R.version.string,implementation='Independent R phyper, p.adjust and cor.test'),file.path(w,'audit/independent_R_verification.json'),auto_unbox=TRUE,pretty=TRUE,digits=14)
