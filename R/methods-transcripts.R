## This is where I will put methods to overload things like
## transcripts() and exons()...

## new argument: columns here can be any legit value for columns. (not just
## tx_id and tx_name etc.)

## filter will just pass through to the internal transcripts call.

## columns arg is just for b/c support and will just pass through to
## the internal transcripts call


## ## For consistency, the helper columns just wraps around cols method...
## setMethod("columns", "MultiDb", function(x){.cols(x)})


## .getTxDb <- function(x){
##     ## trick: there will *always* be a TXID
##     res <- .lookupDbFromKeytype(x, "TXID")
##     if(!is.null(res)){
##         return(res)
##     }else{
##         return(NA)
##     }
## }

.getTxDb <- function(x){
    res <- x@txdbSlot
    if(!is.null(res)){
        return(res)
    }else{
        return(NA)
    }
}


## expose method for gettting A TxDb (if there is one)
setMethod("getTxDbIfAvailable", "MultiDb", function(x, ...){.getTxDb(x)})


## And actually, just make a setter / getter for TxDbs on OrganismDb objects
## In this case, these methods are exclusive to OrganismDb objects.
## getter
setMethod("TxDb", "OrganismDb", function(x, ...){.getTxDb(x)})


## TxDb setter method
## .updateTxDb() helper makes graphInfo from mdb, modifies it to use new
## txdb info and then calls MultiDb() constructor...
.updateTxDb <- function(x, value){
    ## Here I need to work out what needs an update and update it...
    ## I need to find the TxDb in the object and replace it with
    ## the one in value
    if(is(value, 'TxDb')) stop('Replacement value must be a TxDb object.')
    
    ## 1st get the current TxDbs name
    txDbName <- OrganismDbi:::.lookupDbNameFromKeytype(x, 'TXID')
    ## we will use a generated name for internals when user does this.
    newTxDbName <- makePackageName(value)
    
    ## To modify the TxDb value rebuild the MultiDb
    ## 1) Extract/modify the keys/graphData
    gd <- x@keys
    gd[gd %in% txDbName] <- newTxDbName
    ## 2) Extract/modify the resources
    resources <- x@resources
    resources[names(resources) %in% txDbName] <- dbfile(value)
    names(resources)[names(resources) %in% txDbName] <- newTxDbName
    ## 3) rebuild (which should populate any slots etc.)
    graphInfo <- list(graphData=gd, resources=resources)
    x <- OrganismDbi:::OrganismDb(graphInfo=graphInfo)
    ## then return the new MultiDb object.
    x
}
setReplaceMethod("TxDb", "OrganismDb", function(x, value) .updateTxDb(x, value))

## test for setter:
## library(OrganismDbi); example(makeOrganismDbFromTxDb); odb;
## saveDb(txdbMouse, file='myTxDb.sqlite')

## library(TxDb.Mmusculus.UCSC.mm9.knownGene); txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene;

## library(OrganismDbi); txdbMouse <- loadDb('myTxDb.sqlite'); odb <- makeOrganismDbFromTxDb(txdb=txdbMouse)
## debug(OrganismDbi:::.updateTxDb) ## works!
## debug(OrganismDbi:::.getTxDb)
## TxDb(odb) <- txdb;       odb; odb@resources


## OK.  It looks like setter is working, but for this particular
## example, my getter is still grabbing from global scope sometimes?
## (like in this example) - all because the names are being used too
## much by things like getters for the object.  So in this example, I
## am at the mercy of the last thing that was loaded...  :/


## Actually, what I want to do is to make sure that any time I get DB
## resource, I first try to do loadDb on the matching thing from the
## resources slot AND THEN try to call get() on an existing name (but
## as a backup plan!).  - and this seems to work, but there is a bad fail in one unit test that I have to check into now...


## And if I can't solve the problem in that way:
## Then this can also be made more robust if instead I put an actual
## TxDb and OrgDb object into custom slots (with custom getter methods
## for each that act only for OrganismDb objects).  That way I can
## insulate my object from external interference from name clashes.










########################################################################
## TODO: .compressMetadata() might be useful to move into IRanges, as
## a complement to expand() methods?  - Discuss this with Val (who
## apparently may have similar issues in vcf...

## .compressMetadata() processes data.frame data into a DataFrame with
## compressed chars

## It does so by taking a special factor (f) and then applying it to
## ALL of the columns in a data.frame (meta) except for the one that
## was the basis for the special factor (avoidID)
.compressMetadata <- function(f, meta, avoidID){
    meta <- meta[, !colnames(meta) %in% avoidID, drop=FALSE]
    meta <- lapply(meta, function(column, f) unique(splitAsList(column, f)), f)
    DataFrame(meta, row.names=NULL)
}

## This helper does book keeping that is relevant to my situation here.
.combineMetadata <- function(rngs, meta, avoidID, joinID, columns){
    ## compress the metadata by splitting according to f
    joinValue <- as.character(mcols(rngs)[[joinID]])
    f <- factor(meta[[avoidID]], levels=unique(joinValue))
    if (anyNA(f))
        stop("not all annotations have matching ranges")

    if (avoidID %in% columns)           # don't avoid the avoidID
        avoidID <- NULL

    res <- .compressMetadata(f, meta, avoidID)
    ridx <- match(joinValue, levels(f))
    res <- c(mcols(rngs), res[ridx, , drop=FALSE])
    res[!(colnames(res) %in% c("tx_id","exon_id","cds_id","gene_id"))]
}



## How will we merge the results from select() and transcripts()?  We
## will join on tx_id (for transcripts)
.transcripts <- function(x, columns, filter){
    ## 1st get the TxDb object.
    txdb <- .getTxDb(x)
    ## call transcripts method (on the TxDb)
    txs <- transcripts(txdb, columns="tx_id", filter=filter)  
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=as.character(mcols(txs)$tx_id), columns, "TXID")    
    ## assemble it all together.
    mcols(txs) <- .combineMetadata(txs,meta,avoidID="TXID",joinID="tx_id",
                                   columns=columns) 
    txs
}

setMethod("transcripts", "MultiDb",
          function(x, columns=c("TXID", "TXNAME"), filter=NULL){
              .transcripts(x, columns, filter)} )


## test usage:
## library(Homo.sapiens); h = Homo.sapiens; columns = c("TXNAME","SYMBOL")
## transcripts(h, columns)


## How will we merge the results from select() and transcripts()?  We
## will join on tx_id (for transcripts)
.exons <- function(x, columns, filter){
    ## 1st get the TxDb object.
    txdb <- .getTxDb(x)
    
    ## call transcripts method (on the TxDb)
    exs <- exons(txdb, columns="exon_id", filter=filter)
    
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=as.character(mcols(exs)$exon_id), columns, "EXONID")
    
    ## assemble it all together.
    mcols(exs) <- .combineMetadata(exs,meta,avoidID="EXONID",joinID="exon_id",
                                   columns=columns)
    exs
}

setMethod("exons", "MultiDb",
          function(x, columns="EXONID", filter=NULL){
              .exons(x, columns, filter)})


## test usage:
## library(Homo.sapiens); h = Homo.sapiens; columns = c("CHR","REFSEQ")
## exons(h, columns)


## How will we merge the results from select() and transcripts()?  We
## will join on tx_id (for transcripts)
.cds <- function(x, columns, filter){
    ## 1st get the TxDb object.
    txdb <- .getTxDb(x)
    
    ## call transcripts method (on the TxDb)
    cds <- cds(txdb, columns="cds_id", filter=filter)
    
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=as.character(mcols(cds)$cds_id), columns, "CDSID")
    
    ## assemble it all together.
    mcols(cds) <- .combineMetadata(cds,meta,avoidID="CDSID",joinID="cds_id",
                                   columns=columns)
    cds
}

setMethod("cds", "MultiDb",
          function(x, columns="CDSID", filter=NULL){
              .cds(x, columns, filter)})


## test usage:
## library(Homo.sapiens); h = Homo.sapiens; columns = c("GENENAME","SYMBOL")
## cds(h, columns)





## How will we merge the results from select() and transcripts()?  We
## will join on tx_id (for transcripts)
.genes <- function(x, columns, filter){
    ## 1st get the TxDb object.
    txdb <- .getTxDb(x)
    
    ## call transcripts method (on the TxDb)
    genes <- genes(txdb, columns="gene_id", filter=filter)
    
    ## call select on the rest and use tx_id as keys 
    meta <- select(x, keys=as.character(mcols(genes)$gene_id), columns,
                   "GENEID")
    
    ## assemble it all together.
    mcols(genes) <- .combineMetadata(genes,meta,avoidID="GENEID",
                                     joinID="gene_id",
                                     columns=columns)
    genes
}

setMethod("genes", "MultiDb",
          function(x, columns="GENEID", filter=NULL){
              .genes(x, columns, filter)})


## test usage:
## library(Homo.sapiens); h = Homo.sapiens; columns = c("GENENAME","SYMBOL")
## genes(h, columns)




########################################################################
########################################################################
##                       The "By" methods
########################################################################
########################################################################

## "By" methods will just cram the same metadata into the INTERNAL
## metadata slot so that it appears with the show method.
## No attempt will be made to manage the insanity of knowing which
## metadata types belong in which spot...

.byToKeytype <- function(by){
    switch(by,
           'gene'='GENEID',
           'exon'='EXONID',
           'cds'='CDSID',
           'tx'='TXID')
}

.transcriptsBy <- function(x, by, columns, use.names, outerMcols){
    ## 1st get the TxDb object.
    txdb <- .getTxDb(x)
    txby <- transcriptsBy(txdb, by=by, use.names=use.names)

    if(length(columns) >= 1){ 
        ## get the tx_ids from the transcripts
        ## AND I need to one from the internal slot.
        gr <- txby@unlistData
        k  <- as.character(mcols(gr)$tx_id)
 
        ## call select on the rest and use tx_id as keys 
        meta <- select(x, keys=k, columns, "TXID") 
        ## assemble it all together.
        mcols(gr) <- .combineMetadata(gr, meta, avoidID="TXID",
                                      joinID="tx_id", columns=columns) 
        ## now cram it back in there.
        txby@unlistData <- gr
        ## AND ALSO put the metadata in for the 'outer' mcols...
        if(outerMcols==TRUE){
            k2 <- names(txby)
            keytype <- .byToKeytype(by)
            meta2 <- select(x, keys=k2, columns, keytype)
            ## Step here needed to make meta2 from data.frame into DataFrame
            f <- factor(meta2[[keytype]], levels=unique(as.character(k2)))
            metaC <- .compressMetadata(f, meta2, avoidID=keytype)
            mcols(txby) <- metaC
        }
    }
    txby
}

setMethod("transcriptsBy", "MultiDb",
          function(x, by="gene", columns=character(), use.names=FALSE,
                   outerMcols=FALSE){
              if(missing(by) || !any(by %in% c("gene","exon","cds")) ||
                 length(by) !=1){
                  stop("You must provide a valid argument for by")}
              .transcriptsBy(x, by, columns, use.names=use.names,
                             outerMcols=outerMcols)})


.exonsBy <- function(x, by, columns, use.names, outerMcols){
    ## 1st get the TxDb object.
    txdb <- .getTxDb(x)
    exby <- exonsBy(txdb, by=by, use.names=use.names)

    if(length(columns) >= 1){ 
        ## get the tx_ids from the transcripts
        ## AND I need to one from the internal slot.
        gr <- exby@unlistData
        k  <- as.character(mcols(gr)$exon_id)
    
        ## call select on the rest and use tx_id as keys 
        meta <- select(x, keys=k, columns, "EXONID")    
        ## assemble it all together.
        mcols(gr) <- .combineMetadata(gr, meta, avoidID="EXONID",
                                      joinID="exon_id", columns=columns) 
        ## now cram it back in there.
        exby@unlistData <- gr
        ## AND ALSO put the metadata in for the 'outer' mcols...
        if(outerMcols==TRUE){
            k2 <- names(exby)
            keytype <- .byToKeytype(by)
            meta2 <- select(x, keys=k2, columns, keytype)
            ## Step here needed to make meta2 from data.frame into DataFrame
            f <- factor(meta2[[keytype]], levels=unique(as.character(k2)))
            metaC <- .compressMetadata(f, meta2, avoidID=keytype)
            mcols(exby) <- metaC
        }
    }
    exby
}

setMethod("exonsBy", "MultiDb",
          function(x, by="tx", columns=character(), use.names=FALSE,
                   outerMcols=FALSE){
              if(missing(by) || !any(by %in% c("tx", "gene")) ||
                 length(by) !=1){
                  stop("You must provide a valid argument for by")}
              .exonsBy(x, by, columns, use.names=use.names,
                       outerMcols=outerMcols)})

.cdsBy <- function(x, by, columns, use.names, outerMcols){
    ## 1st get the TxDb object.
    txdb <- .getTxDb(x)
    cdsby <- cdsBy(txdb, by=by, use.names=use.names)

    if(length(columns) >= 1){ 
        ## get the tx_ids from the transcripts
        ## AND I need to one from the internal slot.
        gr <- cdsby@unlistData
        k  <- as.character(mcols(gr)$cds_id)
 
        ## call select on the rest and use tx_id as keys 
        meta <- select(x, keys=k, columns, "CDSID") 
        ## assemble it all together.
        mcols(gr) <- .combineMetadata(gr, meta, avoidID="CDSID",
                                      joinID="cds_id", columns=columns) 
        ## now cram it back in there.
        cdsby@unlistData <- gr
        ## AND ALSO put the metadata in for the 'outer' mcols...
        if(outerMcols==TRUE){
            k2 <- names(cdsby)
            keytype <- .byToKeytype(by)
            meta2 <- select(x, keys=k2, columns, keytype)
            ## Step here needed to make meta2 from data.frame into DataFrame
            f <- factor(meta2[[keytype]], levels=unique(as.character(k2)))
            metaC <- .compressMetadata(f, meta2, avoidID=keytype)
            mcols(cdsby) <- metaC
        }
    }
    cdsby
}

setMethod("cdsBy", "MultiDb",
          function(x, by="tx", columns=character(), use.names=FALSE,
                   outerMcols=FALSE){
              if(missing(by) || !any(by %in% c("tx", "gene")) ||
                 length(by) !=1){
                  stop("You must provide a valid argument for by")}
              .cdsBy(x, by, columns, use.names=use.names,
                     outerMcols=outerMcols)})

## TODO: (known issues)
## 1) columns don't come back in same order that the went in

## 2) some values (tx_id and tx_name come to mind) are not relabeled
## in a pretty way and may not have been requested (to solve this we
## have to adress issue #3)

## 3) I now have a columns AND a columns argument for the transcripts()
## family of methods.  This is totally redundant.  Proposed fix:
## rename arguments base method to be columns (maybe this is also an
## opportunity to rename columns everywhere), but rename it so that it's
## consistent, and then here, just only have one argument...

## 4) exonsBy and cdsBy may have some extra issues that I am missing...


###############################################################################
## Wrapper methods for TxDb methods (that are the same)



setMethod(promoters, 'MultiDb',
          function(x, upstream=2000, downstream=200, ...){
              promoters(getTxDbIfAvailable(x), upstream,
                        downstream, ...)})

setMethod(microRNAs, 'MultiDb',
          function(x){microRNAs(getTxDbIfAvailable(x))})

setMethod(tRNAs, 'MultiDb',
          function(x){tRNAs(getTxDbIfAvailable(x))})

setMethod(intronsByTranscript, 'MultiDb',
          function(x, use.names=FALSE){
              intronsByTranscript(getTxDbIfAvailable(x),
                                  use.names=use.names)})

setMethod(fiveUTRsByTranscript, 'MultiDb',
          function(x, use.names=FALSE){
              fiveUTRsByTranscript(getTxDbIfAvailable(x),
                                   use.names=use.names)})

setMethod(threeUTRsByTranscript, 'MultiDb',
          function(x, use.names=FALSE){
              threeUTRsByTranscript(getTxDbIfAvailable(x),
                                    use.names=use.names)})

setMethod(extractUpstreamSeqs, 'MultiDb',
          function(x, genes, width=1000, exclude.seqlevels=NULL){
              extractUpstreamSeqs(x, getTxDbIfAvailable(genes),
                                  width=width,
                                  exclude.seqlevels=exclude.seqlevels)})

## problem: no way for this to dispatch correctly...
## So we are basically trampling the original method definition
## here. (which is just NOT elegant)
## But: we plan to move OrganismDbi down into TxDbs later on so this
## is temporary.
## setMethod(extractTranscriptSeqs, 'BSgenome',
##           function(x, transcripts, strand="+"){
##               if(class(transcripts)=='MultiDb'){
##                   transcripts<-getTxDbIfAvailable(transcripts)
##               }
##               extractTranscriptSeqs(x, transcripts, strand=strand)})

## This works now without the need for an overload...
## library(Homo.sapiens);library(BSgenome.Hsapiens.UCSC.hg19);genome <- BSgenome.Hsapiens.UCSC.hg19;debug(GenomicFeatures:::.normarg_transcripts);tx_seqs <- extractTranscriptSeqs(genome, Homo.sapiens)


setMethod(isActiveSeq, 'MultiDb',
          function(x){isActiveSeq(getTxDbIfAvailable(x))})

.updateTxDbSeqMultiDb <-function(x, value){
    ## This will change the val in 'x' as well...
    txdb <- getTxDbIfAvailable(x)
    if (is(txdb, "TxDb")) { ## will be NA if there isn't one.
        isActiveSeq(txdb) <- value
    } else {
        stop('This object does not contain a TxDb object')
    }
    x
}
setReplaceMethod('isActiveSeq', 'MultiDb',
          function(x, value){.updateTxDbSeqMultiDb(x, value)})

setMethod(asBED, 'MultiDb', function(x){asBED(getTxDbIfAvailable(x))})
setMethod(asGFF, 'MultiDb', function(x){asGFF(getTxDbIfAvailable(x))})

## These ones dispatch on compound types (not just TxDbs):
setMethod(distance, c('GenomicRanges','MultiDb'),
          function(x, y, ignore.strand=FALSE, ..., id, 
                   type=c("gene", "tx", "exon", "cds")){
              distance(x, getTxDbIfAvailable(y), ignore.strand=ignore.strand,
                       ..., id=id, type=type)})

setMethod(mapToTranscripts, c('ANY', 'MultiDb'),
          function(x, transcripts, ignore.strand=TRUE, 
                   extractor.fun = GenomicFeatures::transcripts, ...){
              mapToTranscripts(x, transcripts=getTxDbIfAvailable(transcripts),
                               ignore.strand=ignore.strand, 
                               extractor.fun = extractor.fun, ... )})

