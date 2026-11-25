args <- commandArgs(trailingOnly=TRUE)
x <- readRDS(args[1])
x <- x[grepl('^GO:', x[['gs_subcollection']]), ]
k <- paste(x[['gs_subcollection']], x[['gs_name']], sep='|')
gg <- split(x[['db_gene_symbol']],k)
out <- data.frame(collection=sub('[|].*','',names(gg)),term=sub('^[^|]*[|]','',names(gg)),genes=vapply(gg,function(y)paste(sort(unique(y)),collapse=';'),character(1)))
write.table(out,args[2],sep='\t',row.names=FALSE,quote=FALSE)
print(table(out$collection))
