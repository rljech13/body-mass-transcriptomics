suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(jsonlite)})
a<-commandArgs(trailingOnly=TRUE);w<-normalizePath(a[1]);root<-normalizePath(file.path(w,'../..'))
.libPaths(c(file.path(root,'outputs/mammal_R_SVG_2026-10-23/R_library'),.libPaths()))
d<-read.delim(file.path(w,'data/species_gene_expression.tsv'),check.names=FALSE)
sel<-read.delim(file.path(w,'data/selected_targets.tsv'),check.names=FALSE)
sp<-read.delim(file.path(w,'data/species_order.tsv'),check.names=FALSE)
src<-file.path(root,'outputs/mammal_SCENIC_2026-10-23/selected/inputs')
x<-read.delim(gzfile(file.path(src,'expression.tsv.gz')),row.names=1,check.names=FALSE)
meta<-read.delim(file.path(src,'samples.tsv'),check.names=FALSE)
spp<-meta$species[match(rownames(x),meta$sample)]
stopifnot(nrow(d)==300,nrow(sel)==20,nrow(sp)==15,!anyDuplicated(sel$row_id),!anyNA(spp),!anyDuplicated(sel$gene),all(sel$gene %in% colnames(x)))
for(tf in c('ERG','FOXP4')) {
 e<-read.delim(file.path(root,'outputs/mammal_SCENIC_regulon_screen_2026-10-23/data',paste0(tf,'_all_targets.tsv')))
 e<-e[order(-e$importance,e$target),];s<-sel[sel$TF==tf,]
 stopifnot(identical(s$gene,e$target[1:10]),max(abs(s$importance-e$importance[1:10]))<1e-12)
}
mat<-sapply(sel$gene,function(g)vapply(sp$species,function(s)mean(x[spp==s,g]),numeric(1)))
z<-scale(mat,center=TRUE,scale=TRUE)
ix<-cbind(match(d$species,sp$species),match(d$gene,sel$gene))
mean_error<-max(abs(mat[ix]-d$species_mean_log2_CPM_plus1))
z_error<-max(abs(z[ix]-d$row_z))
stopifnot(mean_error<1e-12,z_error<1e-12,all(diff(sp$mass_g)>0))
write_json(list(implementation='Independent R reading original wide frozen expression, matching sample metadata, independently checking saved top-ten ranking, species means and sample-SD row z.',mean_max_abs_error=mean_error,z_max_abs_error=z_error,verified_cells=nrow(d),all_selected_weights_and_ranks_match=TRUE),file.path(w,'audit/independent_R_verification.json'),pretty=TRUE,auto_unbox=TRUE,digits=16)
d$species<-factor(d$species,levels=sp$species);sp$species<-factor(sp$species,levels=sp$species)
d$gene_row<-factor(d$row_id,levels=rev(sel$row_id));d$TF<-factor(d$TF,levels=c('ERG','FOXP4'))
ylabs<-setNames(sel$gene,sel$row_id)
lim<-ceiling(max(abs(d$row_z)))
hm<-ggplot(d,aes(species,gene_row,fill=row_z))+
 geom_tile(width=.98,height=.98)+
 facet_grid(TF~.,scales='free_y',space='free_y',switch='y',labeller=as_labeller(c(ERG='ERG targets',FOXP4='FOXP4 targets')))+
 scale_fill_gradient2(name='Expression\n(row z-score)',low='#2166AC',mid='white',high='#B2182B',midpoint=0,limits=c(-lim,lim),breaks=seq(-lim,lim,1),guide=guide_colorbar(barheight=grid::unit(3.1,'cm'),barwidth=grid::unit(.35,'cm')))+
 scale_y_discrete(labels=ylabs,expand=expansion(add=.05))+
 scale_x_discrete(expand=expansion(add=.05))+
 labs(x='Species ordered by increasing body mass →',y=NULL)+
 theme_light(base_size=11,base_family='sans')+
 theme(axis.text.x=element_text(angle=55,hjust=1,vjust=1,size=9.5,color='black'),axis.text.y=element_text(size=10.5,color='black',face='italic'),axis.title.x=element_text(size=11,margin=margin(t=10)),axis.ticks=element_blank(),panel.grid=element_blank(),panel.border=element_blank(),strip.placement='outside',strip.background=element_blank(),strip.text.y.left=element_text(angle=90,size=11,face='bold',color='black'),panel.spacing.y=grid::unit(.7,'lines'),legend.title=element_text(size=10),legend.text=element_text(size=9),legend.position='right',plot.margin=margin(0,10,8,8))
mass<-ggplot(sp,aes(species,1,fill=log10_mass))+
 geom_tile(width=.98,height=.75)+
 scale_fill_gradient(name='Log10 mass (g)',low='#EEEEEE',high='#333333',limits=c(1,4.5),breaks=c(1,2,3,4),guide=guide_colorbar(barwidth=grid::unit(.35,'cm'),barheight=grid::unit(1.8,'cm')))+
 scale_x_discrete(expand=expansion(add=.05))+scale_y_continuous(limits=c(.5,1.5),breaks=1,labels='Body mass')+
 labs(x=NULL,y=NULL)+theme_void(base_size=11,base_family='sans')+
 theme(axis.text.y=element_text(color='black',size=10),legend.position='right',legend.title=element_text(size=10),legend.text=element_text(size=9),plot.margin=margin(0,10,0,8))
p<-(mass/hm)+plot_layout(heights=c(.75,10))+
 plot_annotation(title='F',theme=theme(plot.title=element_text(size=21,face='bold',hjust=0,margin=margin(b=10)),plot.margin=margin(10,12,8,8)))
svglite::svglite(file.path(w,'Panel_F_regulon_targets.svg'),width=8.8,height=8.5,bg='white');print(p);dev.off()
ggsave(file.path(w,'audit/preview.png'),p,width=8.8,height=8.5,dpi=170,bg='white')
write_json(list(R=R.version.string,packages=setNames(lapply(c('ggplot2','patchwork','svglite','jsonlite'),function(n)as.character(packageVersion(n))),c('ggplot2','patchwork','svglite','jsonlite')),source='Frozen primary Ma2016 SCENIC input and nonself target weights, seed 20261023',rows=20,columns=15,full_color_range=c(-lim,lim),actual_range=range(d$row_z),clipped_cells=sum(abs(d$row_z)>lim),rows_clustered=FALSE,species_reordered_by_mass=TRUE,signs_flipped=FALSE,new_models=FALSE,new_association_tests=FALSE,figure_approval='pending',presentation_modified=FALSE),file.path(w,'audit/render.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
