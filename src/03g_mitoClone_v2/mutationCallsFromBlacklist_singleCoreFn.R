# Edit the mutationCallsFromBlacklist function to remove multiple cores, which
# caused errors on my local computer
mutationCallsFromBlacklist <- function(BaseCounts,lim.cov=20, min.af=0.2, 
                                       min.num.samples=0.01*length(BaseCounts), 
                                       min.af.universal =min.af, 
                                       universal.var.cells=0.95*length(BaseCounts), 
                                       blacklists.use = blacklists, max.var.na = 0.5, 
                                       max.cell.na = 0.95, ...) {
  varaf <- parallel::mclapply(BaseCounts,function(x){
    ## focus on A,G,C,T
    x <- x[,1:4]
    ## find cell that have less than 20 cov over agct at a given pos
    zeroes <- rowSums(x) < lim.cov
    ## af calc
    #x.af <- x/rowSums(x)
    x.af <- x / (x+apply(x,1,max))
    x.af <- reshape2::melt(x.af)
    colnames(x.af) <- c('pos','nt','af')
    ## remove reference af's
    x.af <- x.af[!(mito.dna[x.af$pos] == x.af$nt),]
    ## remove N site
    x.af <- x.af[!(mito.dna[x.af$pos] == 'N'),]
    x.af$name <- paste0(x.af$pos,' ',mito.dna[x.af$pos],'>',x.af$nt)
    ## find dominant NT
    x.af$af[x.af$pos %in% which(zeroes)] <- NA
    x <- x.af$af
    names(x) <- x.af$name
    return(x)
  }, mc.cores=1) ####################### remove parallelism here ###############
  varaf <- do.call(cbind, varaf)
  ## you could allow for only sites with coverage! currently you filter at a rate of 10% cells dropping out max
  ##varaf <- varaf[rowSums(is.na(varaf))/length(mc.out) < max.fraction.na,]
  varaf <- varaf[rowSums(varaf > min.af, na.rm=TRUE) >= min.num.samples,]
  
  is.names <- sapply(blacklists.use, function(x) typeof(x) == "character")
  #part 2 - filter based on the blacklist
  if(sum(is.names) > 0){
    removal.names.list <- unique(unlist(blacklists.use[is.names]))
    varaf <- varaf[!row.names(varaf) %in% removal.names.list,]
  }
  if(sum(!is.names) > 0){
    removal.ranges.list <- unique(unlist(GenomicRanges::GRangesList(blacklists.use[!is.names])))
    varaf <- varaf[-c(S4Vectors::queryHits(GenomicRanges::findOverlaps(mut2gr(row.names(varaf)),removal.ranges.list))),]
  }
  #if(drop.empty){
  varaf <- varaf[rowSums(varaf,na.rm=T) > 0,] #colSums(varaf,na.rm=T) > 0
  #}
  
  varaf <- varaf[!rowSums(varaf >= min.af.universal,na.rm=TRUE) >= universal.var.cells,]
  ## vars must have less than X % NA's
  varaf <- varaf[rowSums(is.na(varaf)) < max.var.na*NCOL(varaf),]
  ## cells must have less than X % NA's
  varaf <- varaf[,colSums(is.na(varaf)) < max.cell.na*NROW(varaf)]
  
  MN <- pullcounts.vars(BaseCounts, rownames(varaf), colnames(varaf))
  mutationCallsFromMatrix(t(MN$M), t(MN$N), ...)
}