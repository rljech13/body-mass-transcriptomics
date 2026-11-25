suppressPackageStartupMessages({library(ggplot2);library(jsonlite)})
a<-commandArgs(trailingOnly=TRUE);w<-normalizePath(a[1]);root<-normalizePath(file.path(w,'../..'));.libPaths(c(file.path(root,'outputs/mammal_R_SVG_2026-10-23/R_library'),.libPaths()))
d<-read.delim(file.path(w,'selected/results/analysis/module_regulon_results.tsv'),check.names=FALSE);mods<-c('orange','violet','cyan','skyblue3')
# Display the two leading primary candidates per module plus the preexisting E2F2 candidate.
chosen<-unique(unlist(lapply(mods,function(m){x<-d[d$module==m,];x<-x[order(x$q,x$p,x$TF),];head(x$TF,2)})))
if('E2F2'%in%d$TF)chosen<-unique(c(chosen,'E2F2'))
best<-sapply(chosen,function(tf)min(d$q[d$TF==tf]));chosen<-chosen[order(best,chosen)];seeds<-c(20261023,20261024,20261025);xx<-do.call(rbind,lapply(seq_along(seeds),function(i){a<-read.delim(file.path(w,'selected/results/analysis',paste0('overlap_',seeds[i],'.tsv')));a$run<-i;a[a$TF%in%chosen,]}));grid<-expand.grid(TF=chosen,module=mods,run=1:3,stringsAsFactors=FALSE);x<-merge(grid,xx,by=c('TF','module','run'),all.x=TRUE);x$TF<-factor(x$TF,levels=rev(chosen));x$module<-factor(x$module,levels=mods);x$run<-factor(x$run,levels=1:3);x$FDR<-factor(ifelse(x$q<.05,'< 0.05','≥ 0.05'),levels=c('< 0.05','≥ 0.05'));labels<-c(orange='orange\n77 genes',violet='violet\n56 genes',cyan='cyan\n116 genes',skyblue3='skyblue3\n44 genes')
p<-ggplot(x,aes(run,TF))+
 geom_point(data=x[!is.na(x$q)&x$k_shared>0,],aes(size=k_shared,fill=FDR),shape=21,color='#555555',stroke=.3)+
 geom_point(data=x[is.na(x$q),],shape=4,size=1.3,color='#AAAAAA',stroke=.4)+
 facet_grid(.~module,labeller=labeller(module=labels),drop=FALSE)+
 scale_size_area(name='Shared genes',max_size=5.8,breaks=c(1,5,10,20))+
 scale_fill_manual(name='FDR q',values=c('< 0.05'='#0072B2','≥ 0.05'='#D9D9D9'),drop=FALSE)+
 scale_x_discrete(drop=FALSE)+
 labs(title='D',x='SCENIC run',y=NULL)+
 theme_light(base_size=12,base_family='sans')+
 theme(plot.title=element_text(size=21,face='bold',hjust=0,margin=margin(b=8)),axis.title=element_text(size=12.5),axis.text=element_text(size=10.5,color='black'),strip.background=element_rect(fill='white',color='#CCCCCC'),strip.text=element_text(size=10.5,color='black'),panel.grid.minor=element_blank(),panel.grid.major=element_line(color='#EEEEEE',linewidth=.3),legend.position='bottom',legend.box='vertical',legend.title=element_text(size=10.5),legend.text=element_text(size=10),plot.margin=margin(12,12,10,8))+
 guides(size=guide_legend(order=1,nrow=1),fill=guide_legend(order=2,nrow=1,override.aes=list(size=4)))
dir.create(file.path(w,'figures'),showWarnings=FALSE);svglite::svglite(file.path(w,'figures/Panel_D_SCENIC.svg'),width=6.5,height=6.5,bg='white');print(p);dev.off();ggsave(file.path(w,'audit/Panel_D_SCENIC_preview.png'),p,width=6.5,height=6.5,dpi=170,bg='white')
write.table(x,file.path(w,'data/Panel_D_SCENIC_points.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write_json(list(R=R.version.string,packages=sapply(c('ggplot2','svglite'),function(x)as.character(packageVersion(x))),source='full pySCENIC three predeclared seeds on original Ma2016 profiles; existing hg38 cisTarget database',displayed_TFs=chosen,selection='two leading primary module-regulon tests per fixed module plus prespecified E2F2; all primary187x4 tests remain in BH family; repeat families use all195x4 and172x4 tests',point_size='number of shared genes',point_fill='within-run module overlap BH q<0.05; three runs shown separately; x means TF regulon absent in that run',figure_approval='pending',presentation_modified=FALSE),file.path(w,'audit/render.json'),auto_unbox=TRUE,pretty=TRUE,digits=14)
