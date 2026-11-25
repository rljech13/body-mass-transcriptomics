suppressPackageStartupMessages({library(ggplot2);library(ggrepel);library(jsonlite)})
a<-commandArgs(trailingOnly=TRUE);w<-normalizePath(a[1]);root<-normalizePath(file.path(w,'../..'))
.libPaths(c(file.path(root,'outputs/mammal_R_SVG_2026-10-23/R_library'),.libPaths()))
src<-file.path(root,'outputs/mammal_SCENIC_2026-10-23/selected')
auc<-read.delim(file.path(src,'results/analysis/AUCell_species.tsv'),check.names=FALSE)
traits<-read.delim(file.path(src,'inputs/species_traits.tsv'),check.names=FALSE)
checked<-read.delim(file.path(root,'outputs/mammal_WGCNA_GO_2026-10-23/data/ERG_FOXP4_verified.tsv'),check.names=FALSE)
stats<-checked[checked$seed==20261023,];stopifnot(nrow(stats)==2,nrow(auc)==15,!anyDuplicated(auc$species))
traits<-traits[match(auc$species,traits$species),];stopifnot(!anyNA(traits$mass_g),all(traits$mass_g>0))
rows<-list();fits<-list();cols<-c(ERG='#2166AC',FOXP4='#984EA3')
for(tf in c('ERG','FOXP4')) {
 n<-paste0(tf,'(+)');st<-stats[stats$regulon==n,];stopifnot(nrow(st)==1)
 d<-data.frame(TF=tf,species=auc$species,mass_g=traits$mass_g,log10_mass=log10(traits$mass_g),activity=auc[[n]])
 stopifnot(abs(cor(d$log10_mass,d$activity)-st$mass_r)<1e-10)
 fit<-lm(activity~log10_mass,data=d)
 fits[[tf]]<-data.frame(TF=tf,intercept=unname(coef(fit)[1]),slope=unname(coef(fit)[2]),r=st$mass_r,p=st$mass_p,q=st$mass_q,targets_nonself=st$targets_nonself)
 d$facet<-sprintf('%s · %d targets\nr = %.3f   q = %.4f',tf,st$targets_nonself,st$mass_r,st$mass_q)
 rows[[tf]]<-d
}
d<-do.call(rbind,rows);fits<-do.call(rbind,fits);d$facet<-factor(d$facet,levels=unique(d$facet))
xlim<-c(.85,4.5)
lines<-do.call(rbind,lapply(names(rows),function(tf){d<-rows[[tf]];f<-fits[fits$TF==tf,];data.frame(TF=tf,facet=d$facet[1],log10_mass=xlim,activity=f$intercept+f$slope*xlim)}))
p<-ggplot(d,aes(log10_mass,activity,color=TF))+
 geom_line(data=lines,linewidth=.8)+geom_point(size=2.6)+
 geom_text_repel(aes(label=species),seed=20261023,size=3.1,color='black',max.overlaps=Inf,box.padding=.4,point.padding=.3,min.segment.length=0,max.time=5,max.iter=100000,segment.color='#999999',segment.size=.3)+
 facet_wrap(~facet,nrow=1,scales='free_y')+
 scale_color_manual(values=cols,guide='none')+
 scale_x_continuous(limits=xlim,breaks=seq(1,4.5,.5),expand=expansion(mult=c(.01,.01)))+
 scale_y_continuous(breaks=function(x)pretty(x,n=6),expand=expansion(mult=c(.15,.18)))+
 labs(title='E',x='Log10 body mass (g)',y='Regulon activity (AUCell)')+
 theme_light(base_size=12,base_family='sans')+
 theme(plot.title=element_text(size=21,face='bold',hjust=0,margin=margin(b=9)),axis.title=element_text(size=12.5),axis.title.x=element_text(margin=margin(t=9)),axis.text=element_text(size=10.5,color='black'),panel.grid.minor=element_blank(),strip.background=element_blank(),strip.text=element_text(size=12,color='black',lineheight=1.25,margin=margin(b=12)),panel.spacing.x=grid::unit(2,'lines'),plot.margin=margin(12,20,10,10))
svglite::svglite(file.path(w,'Panel_E_ERG_FOXP4.svg'),width=12,height=6.4,bg='white');print(p);dev.off()
ggsave(file.path(w,'audit/preview.png'),p,width=12,height=6.4,dpi=160,bg='white')
write.table(d,file.path(w,'data/species_activity.tsv'),sep='\t',quote=TRUE,row.names=FALSE)
write.table(fits,file.path(w,'data/display_fits.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write_json(list(R=R.version.string,packages=setNames(lapply(c('ggplot2','ggrepel','svglite','jsonlite'),function(n)as.character(packageVersion(n))),c('ggplot2','ggrepel','svglite','jsonlite')),study='Ma2016',species=15,primary_seed=20261023,facets=c('ERG','FOXP4'),source_values='Previously independently verified primary species-mean AUCell and mass r/p/q; q across all 187 primary regulons.',new_SCENIC=FALSE,new_association_tests=FALSE,display_fits='Ordinary least squares AUCell against log10 mass; lines extend to entire displayed x range.',different_y_scales=TRUE,own_TF_in_signature=list(ERG=TRUE,FOXP4=FALSE),figure_approval='pending',presentation_modified=FALSE),file.path(w,'audit/render.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
