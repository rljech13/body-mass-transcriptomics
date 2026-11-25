fit <- function(run_dir=NULL, argv=character()) {
  suppressMessages(library(WGCNA))
  options(stringsAsFactors=FALSE);allowWGCNAThreads(nThreads=4);set.seed(20261021)
  O <- run_dir
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
}

loo <- function(run_dir=NULL, argv=character()) {
  suppressMessages(library(WGCNA));options(stringsAsFactors=FALSE);allowWGCNAThreads(nThreads=4);set.seed(20261021)
  O <- run_dir
  rows<-list();j<-0
  for(panel in c("ma2016_full")){
   x<-t(as.matrix(read.delim(file.path(O,"inputs",paste0(panel,".log.tsv")),row.names=1,check.names=FALSE)))
   meta<-read.delim(file.path(O,"inputs",paste0(panel,".samples.tsv")),check.names=FALSE);names(meta)[1]<-"sample"
   full<-read.delim(file.path(O,"results/wgcna",paste0(panel,"__raw"),"members.tsv"))
   run<-read.delim(file.path(O,"results/wgcna",paste0(panel,"__raw"),"run.tsv"))
   for(sp in unique(meta$species)){
    xs<-x[meta$species[match(rownames(x),meta$sample)]!=sp,,drop=FALSE]
    good<-goodSamplesGenes(xs,verbose=0);xs<-xs[,good$goodGenes,drop=FALSE]
  if(ncol(xs)>5000){v<-apply(xs,2,mad);anchor<-scan(file.path(O,"config/forced_anchors.txt"),what="character",quiet=TRUE);chosen<-unique(c(colnames(xs)[order(-v,names(v))[1:5000]],intersect(anchor,colnames(xs))));xs<-xs[,chosen,drop=FALSE]}
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
    write.table(do.call(rbind,rows),file.path(O,"results/wgcna_full_LOO_membership.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
   }
  }
}

avian <- function(run_dir=NULL, argv=character()) {
  #!/usr/bin/env Rscript
  # WGCNA on the rebuilt matrices of script 68 (common_only, window 5000 up / 1000 down).
  #
  # Same networks as script 64, written next to them rather than over them: the
  # old run on 2890 genes stays as the comparison. Per-species networks are skipped
  # here — they were a sensitivity on 7-11 samples, below the package floor of 15.

  suppressMessages(library(WGCNA))
  options(stringsAsFactors = FALSE)
  enableWGCNAThreads(nThreads = 16)

  IN  <- file.path(run_dir, "results/oxphos_species_thread/wgcna_input_v2")
  OUT <- file.path(run_dir, "results/oxphos_species_thread/wgcna_out_v2")
  dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
  TAG <- "common_only_u5000d1000"

  read_m <- function(f) {
    m <- read.table(file.path(IN, f), sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
    as.matrix(m)
  }
  samples <- read.table(file.path(IN, sprintf("samples_%s.tsv", TAG)), sep = "\t", header = TRUE)

  pick_power <- function(datExpr, label) {
    powers <- c(1:10, seq(12, 30, 2))
    sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed hybrid",
                             verbose = 0)
    ok <- sft$fitIndices$SFT.R.sq >= 0.8
    n <- nrow(datExpr)
    fallback <- if (n < 20) 18 else if (n < 30) 16 else if (n < 40) 14 else 12
    p <- if (any(ok)) sft$fitIndices$Power[which(ok)[1]] else fallback
    write.table(cbind(label = label, sft$fitIndices), file.path(OUT, "soft_threshold.tsv"),
                sep = "\t", row.names = FALSE, append = file.exists(file.path(OUT, "soft_threshold.tsv")),
                col.names = !file.exists(file.path(OUT, "soft_threshold.tsv")), quote = FALSE)
    cat(sprintf("[%s] n=%d генов=%d   softPower=%d%s\n", label, n, ncol(datExpr), p,
                if (any(ok)) "" else "  (порог R2=0.8 не достигнут, значение из таблицы WGCNA)"))
    p
  }

  run_net <- function(m, label, minModuleSize = 30, center_by = NULL) {
    if (!is.null(center_by)) {
      for (lv in unique(center_by)) {
        j <- which(center_by == lv)
        m[, j] <- m[, j] - rowMeans(m[, j])
      }
    }
    datExpr <- t(m)
    gsg <- goodSamplesGenes(datExpr, verbose = 0)
    if (!gsg$allOK) {
      datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
      cat(sprintf("[%s] отброшено генов %d, образцов %d\n", label,
                  sum(!gsg$goodGenes), sum(!gsg$goodSamples)))
    }
    p <- pick_power(datExpr, label)
    net <- blockwiseModules(datExpr, power = p, networkType = "signed hybrid",
                            TOMType = "signed", minModuleSize = minModuleSize,
                            deepSplit = 2, mergeCutHeight = 0.25,
                            numericLabels = TRUE, pamRespectsDendro = FALSE,
                            maxBlockSize = 20000, verbose = 0)
    colors <- labels2colors(net$colors)
    me <- moduleEigengenes(datExpr, colors)$eigengenes
    kme <- signedKME(datExpr, me)
    gm <- data.frame(gene = colnames(datExpr), module = colors, label_num = net$colors,
                     stringsAsFactors = FALSE)
    gm$kME_own <- sapply(seq_len(nrow(gm)), function(i) {
      cn <- paste0("kME", gm$module[i])
      if (cn %in% colnames(kme)) kme[i, cn] else NA_real_
    })
    write.table(gm, file.path(OUT, sprintf("gene_modules_%s.tsv", label)), sep = "\t",
                row.names = FALSE, quote = FALSE)
    me2 <- data.frame(sample = rownames(datExpr), me, check.names = FALSE)
    write.table(me2, file.path(OUT, sprintf("module_eigengenes_%s.tsv", label)), sep = "\t",
                row.names = FALSE, quote = FALSE)
    sz <- as.data.frame(table(colors), stringsAsFactors = FALSE)
    names(sz) <- c("module", "n_genes")
    write.table(sz, file.path(OUT, sprintf("module_sizes_%s.tsv", label)), sep = "\t",
                row.names = FALSE, quote = FALSE)
    cat(sprintf("[%s] модулей %d (без grey %d), grey %d генов\n", label,
                nrow(sz), sum(sz$module != "grey"),
                ifelse("grey" %in% sz$module, sz$n_genes[sz$module == "grey"], 0)))
    list(datExpr = datExpr, colors = colors, power = p)
  }

  cat("=== сети по тканям через все виды (панель v2) ===\n")
  Mfib <- read_m(sprintf("M_fib_%s.tsv", TAG))
  Mepi <- read_m(sprintf("M_epi_%s.tsv", TAG))
  nf <- run_net(Mfib, "fib")
  ne <- run_net(Mepi, "epi")

  cat("\n=== обе ткани, центрирование внутри ткани ===\n")
  Mall <- read_m(sprintf("M_all_%s.tsv", TAG))
  tis <- samples$tissue[match(colnames(Mall), samples$sample)]
  na <- run_net(Mall, "all_tc", center_by = tis)

  cat("\n=== consensus по семи видам как по семи sets, внутри фибробластов ===\n")
  spp <- unique(samples$species)
  multi <- vector("list", length(spp))
  common <- rownames(Mfib)
  for (i in seq_along(spp)) {
    cols <- grep(paste0("^", spp[i], "\\|fib\\|"), colnames(Mfib))
    multi[[i]] <- list(data = as.data.frame(t(Mfib[common, cols, drop = FALSE])))
  }
  names(multi) <- spp
  multi <- structure(multi, class = "multiData")
  pc <- min(sapply(multi, function(x) nrow(x$data)))
  cat(sprintf("наборов %d, минимальный размер набора %d образцов\n", length(multi), pc))
  cons <- try(blockwiseConsensusModules(multi, power = nf$power, networkType = "signed hybrid",
                                        TOMType = "signed", minModuleSize = 20, deepSplit = 2,
                                        mergeCutHeight = 0.25, numericLabels = TRUE,
                                        maxBlockSize = 20000, verbose = 0), silent = TRUE)
  if (!inherits(cons, "try-error")) {
    cc <- labels2colors(cons$colors)
    write.table(data.frame(gene = common, module = cc), file.path(OUT, "consensus_modules_fib.tsv"),
                sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("consensus: модулей %d (без grey %d)\n", length(unique(cc)),
                length(setdiff(unique(cc), "grey"))))
  } else {
    cat("consensus не сошёлся:", attr(cons, "condition")$message, "\n")
  }

  cat("\n=== сохранность модулей: фибробласты <-> эпибласт ===\n")
  gcommon <- intersect(colnames(nf$datExpr), colnames(ne$datExpr))
  cat(sprintf("общих генов %d\n", length(gcommon)))
  mp_input <- list(fib = list(data = as.data.frame(nf$datExpr[, gcommon])),
                   epi = list(data = as.data.frame(ne$datExpr[, gcommon])))
  col_fib <- nf$colors[match(gcommon, colnames(nf$datExpr))]
  col_epi <- ne$colors[match(gcommon, colnames(ne$datExpr))]
  mp <- try(modulePreservation(mp_input, list(fib = col_fib, epi = col_epi),
                               referenceNetworks = c(1, 2), nPermutations = 100,
                               randomSeed = 20260921, quickCor = 0, verbose = 0), silent = TRUE)
  if (!inherits(mp, "try-error")) {
    for (ref in c(1, 2)) {
      test <- if (ref == 1) 2 else 1
      z <- mp$preservation$Z[[ref]][[test]]
      obs <- mp$preservation$observed[[ref]][[test]]
      out <- data.frame(module = rownames(z), moduleSize = obs$moduleSize,
                        Zsummary = z$Zsummary.pres, medianRank = obs$medianRank.pres)
      f <- file.path(OUT, sprintf("preservation_ref_%s.tsv", c("fib", "epi")[ref]))
      write.table(out, f, sep = "\t", row.names = FALSE, quote = FALSE)
      cat(sprintf("ссылка %s -> проверка %s: модулей %d, Zsummary>10 у %d, <2 у %d\n",
                  c("fib", "epi")[ref], c("fib", "epi")[test], nrow(out),
                  sum(out$Zsummary > 10, na.rm = TRUE), sum(out$Zsummary < 2, na.rm = TRUE)))
    }
  } else {
    cat("modulePreservation не сошлась:", attr(mp, "condition")$message, "\n")
  }
  cat("\nготово, вывод в", OUT, "\n")
}

summary <- function(run_dir=NULL, argv=character()) {
  suppressPackageStartupMessages(library(jsonlite))
  w<-normalizePath(run_dir);rd<-function(n)read.delim(file.path(w,'data',n),check.names=FALSE,stringsAsFactors=FALSE)
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
}

if (sys.nframe() == 0) { a<-commandArgs(TRUE); if(length(a)<2)stop("operation and run directory required"); do.call(get(a[1]),list(run_dir=a[2],argv=a[-c(1,2)])) }
