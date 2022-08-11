## New helper to lookup which org object or package should be used
## based on the taxonomy ID.  It takes a tax ID and returns an appropriate OrgDb object.
.packageTaxIds <- function(){
    c('180454'='org.Ag.eg.db',
      '3702'='org.At.tair.db',
      '9913'='org.Bt.eg.db',
      '9615'='org.Cf.eg.db',
      '9031'='org.Gg.eg.db',
      '9598'='org.Pt.eg.db',
      '511145'='org.EcK12.eg.db',
      '386585'='org.EcSakai.eg.db',
      '7227'='org.Dm.eg.db',
      '9606'='org.Hs.eg.db',
      '10090'='org.Mm.eg.db',
      '9823'='org.Ss.eg.db',
      '10116'='org.Rn.eg.db',
      '9544'='org.Mmu.eg.db',
      '6239'='org.Ce.eg.db',
      '8355'='org.Xl.eg.db',
      '559292'='org.Sc.sgd.db',
      '7955'='org.Dr.eg.db')
}

.taxIdToOrgDb <- function(taxid) {
    ## These are the packaged TaxIds
    packageTaxIds <- .packageTaxIds()
    if (taxid %in% names(packageTaxIds)) {
        pkg <- packageTaxIds[names(packageTaxIds) %in% taxid]
        nmspc <- loadNamespace(pkg)
        res <- get(pkg, nmspc)
    } else {
        ## If we don't have a package, then lets get the taxIds and AHIds
        ## for the hub objects
        loadNamespace("AnnotationHub")
        ah <- AnnotationHub::AnnotationHub()
        ah <- subset(ah, ah$rdataclass=='OrgDb')
        mc <- mcols(ah)[,'taxonomyid', drop=FALSE]
        ## Then just get the object
        AHID <- rownames(mc[mc$taxonomyid==taxid,,drop=FALSE])
        if (!length(AHID))
            stop("no OrgDb package found for taxid ", taxid)
        else res <- ah[[AHID]]
    }
    res
}
## examples:
## .taxIdToOrgDb(9606)
## .taxIdToOrgDb(9986)


## Need another helper to get us from taxID to the OrgDbName...
.taxIdToOrgDbName <- function(taxid) {
    packageTaxIds <- .packageTaxIds()
    if (taxid %in% names(packageTaxIds)) {
        pkg <- packageTaxIds[names(packageTaxIds) %in% taxid]
        nmspc <- loadNamespace(pkg)
        path <- dbfile(get(pkg, nmspc))
        res <- sub("sqlite$", "db", basename(path))
    } else {
        ## If we don't have a package, then lets get the taxIds and AHIds
        ## for the hub objects
        loadNamespace("AnnotationHub")
        ah <- AnnotationHub::AnnotationHub()
        ah <- subset(ah, ah$rdataclass=='OrgDb')
        mc <- mcols(ah)[,c('taxonomyid','title'), drop=FALSE]
        ## Then just get the object
        data <- mc[mc$taxonomyid==taxid,,drop=FALSE]
        res <- sub("sqlite","db", data$title)
    }
    res
}
## examples
## .taxIdToOrgDbName(9606)
## .taxIdToOrgDbName(9986)


#############################################################################

## simplify DB retrievals from metadata table
.getMetaDataValue <- function(db, name){
  con <- AnnotationDbi::dbconn(db)
  res <- dbGetQuery(con,
    paste0("SELECT value FROM metadata WHERE name='", name,"'"))[[1]]
  if(!is.character(res))
      stop("missing metadata table value for: ", name)
  res
}

## function to allow us to convert a list into the inernally preferred form...
.mungeGraphData <- function(graphData){
  ##pkgs <- sapply(graphData, names) ## This is no good.   :(
  pkgs <- matrix(unlist(lapply(graphData, names)), ncol=2, byrow=TRUE)
  keys <- matrix(unlist(graphData), ncol=2, byrow=TRUE)
  graphData <- cbind(pkgs, keys)
  colnames(graphData) <- c("xDbs","yDbs","xKeys","yKeys")
  graphData
}

## early sanity checks for graphData
.testGraphData <- function(graphData){
  if(ncol(graphData) !=4L)
      stop("'graphData' must contain exactly 4 columns.")
}

## test keys for graphData BEFORE we make any objects (which just
## means that we are not going to use data from the @resources slot
## for this function)
.testKeys <- function(fkeys, env){
  instore <- ls(envir = env)
  ss <- split(fkeys, names(fkeys) %in% instore)
  pkgs <- lapply(setNames(nm = names(ss[['TRUE']])), get, envir = env)
  pkgs <- c(pkgs, lapply(setNames(nm = names(ss[['FALSE']])), get))
  pkgs <- pkgs[names(fkeys)]
  res <- unlist(Map(function(x, y) { x %in% columns(y) }, x = fkeys, y = pkgs))
  if(!all(res))
    stop("some foreign keys are not present in their associated databases")
}

## helper to list bioc Annot packages (for filling in things like
## suggests fields)
.biocAnnPackages <- function(){
    availAnns <- as.data.frame(available.packages(
                               contrib.url(repositories()[["BioCann"]],
                               "source")))
    as.character(availAnns[["Package"]])
}


## Helper to set up to just load packages that need loading.
.extractDbFiles <- function(gd, deps){
    pkgs <- unique(names(.extractPkgsAndCols(gd)))
    ## Before we can proceed, we may need to call library on the deps...
    lapply(deps, library, character.only = TRUE)
    files <- unlist(lapply(pkgs, function(x){dbfile(get(x))}))
    setNames(files, pkgs)
}


## We want makeOrganismPackage to be self contained (have all it needs)
## IOW we want to store .sqlite files in a local inst/extdata when
## they are not known packages.
makeOrganismPackage <- function(pkgname,
                                graphData,
                                organism,
                                version,
                                maintainer,
                                author,
                                destDir,
                                license="Artistic-2.0"){
   ## there should only be one template
   template_path <- system.file("OrgPkg-template",package="OrganismDbi")
   ## We need to get a list of dependencies:
   ## 1st convert graphData into a data.frame
   gd <- .mungeGraphData(graphData)
   allDeps <- unique(as.vector(gd[,1:2]))
   ## Filter dependencies to make sure they are really package names
   biocPkgNames <- .biocAnnPackages()
   deps <- allDeps[allDeps %in% biocPkgNames]
   depsStr <- paste(deps,collapse=", ")
   ## We need to define some symbols in order to have the
   ## template filled out correctly.
   symvals <- list(
    PKGTITLE=paste("Annotation package for the",pkgname,"object"),
    PKGDESCRIPTION=paste("Contains the",pkgname,"object",
      "to access data from several related annotation packages."),
    PKGVERSION=version,
    AUTHOR=author,
    MAINTAINER=maintainer,
    LIC=license,
    ORGANISM=organism,
    ORGANISMBIOCVIEW=gsub(" ","_",organism),
    PKGNAME=pkgname,
    DEPENDENCIES=depsStr
   )
   ## Check the graphData object and rename if needed
   .testGraphData(gd)
## Try to call require on all the supporting packages.
## pkgs <- unique(names(.extractPkgsAndCols(gd)))
## for (pkg in pkgs)
##     .initPkg(pkg)

   ## #########################################################################
   ## Extract the dbFile information from each object and store that
   ## into resources

   ## EXCEPT Don't: (I can't really do this for 'packages' as the data
   ## is system specific)

   ## TODO: change this so that it isn't getting a bunch of dbFiles
   ## and then throwing them aways (or so that it's optional with a
   ## parameter or whatever seems appropriate for this function)

   resources <- .extractDbFiles(gd, deps)
   resources <- unlist(lapply(resources, function(x){return("")}))
   ## #########################################################################

   ## Also check that the fkeys are really columns for the graphData
   fkeys <- .extractPkgsAndCols(gd)
   .testKeys(fkeys)

   ## Should never have duplicates
   if (any(duplicated(names(symvals))))
       stop("'symvals' contains duplicated symbols")
   ## All symvals should by single strings (non-NA)
   is_OK <- sapply(symvals, isSingleString)
   if (!all(is_OK)) {
       bad_syms <- paste(names(is_OK)[!is_OK], collapse="', '")
       stop("values for symbols '", bad_syms, "' are not single strings")
   }
   createPackage(pkgname=pkgname,
                 destinationDir=destDir,
                 originDir=template_path,
                 symbolValues=symvals)

   ## Now just do work to save the data.frame (pass that in instead of the
   ## other stuff) in /data as a serialized R file.
   ## There will already be a /data dir in the template
   ## So just save to it:
   graphData <- gd
   graphInfo <- list(graphData=graphData, resources=resources)
   ## create data dir (because R CMD build removes empty dirs)
   ## And then save the data there.
   dir.create(file.path(destDir,pkgname,"data"))
   save(graphInfo, file=file.path(destDir,pkgname,"data","graphInfo.rda"))

   ## Get and other things that need to be saved and stash them into
   ## /inst/extdata
   otherDeps <- allDeps[!allDeps %in% biocPkgNames]
   .saveFromStr <- function(x, file){saveDb(x=get(x), file=file)}
   if(length(otherDeps)>0){
       ## Then we have to save stuff locally
       extPath <- file.path(destDir,pkgname,"inst","extdata")
       dir.create(extPath, recursive=TRUE)
       mapply(.saveFromStr, x=otherDeps, file=file.path(extPath,
                                           paste0(otherDeps,".sqlite")))
   }
}






################################################################################
## Now for some create functions that are more specialized:
## To create OrgansimDb objects from UCSC or from biomaRt
## the initial versions of these will just create the object (and not
## a package)

## Also: there are some pre-agreed upon expectations for this
## function.  It will make a specific kind of OrganismDb object.  One
## that contains certain specific elements (GO, OrgDb and TxDb - to
## start).  And the expectation is that we will always have those
## things when this function finishes.  So the contract is less
## general than before.


########
## Helper to set up to just load packages that need loading.  BUT this
## version of this function is less agressive and doesn't try to load
## a file if it already exists.  This is a different behavior than we
## want for packaging where things should be more strict.
.gentlyExtractDbFiles <- function(gd, deps, env){
    pkgs <- unique(names(.extractPkgsAndCols(gd)))
    not_pkg <- setdiff(ls(envir = env), deps)
    pkgs <- pkgs[!pkgs %in% not_pkg]
    ## Before we can proceed, we may need to call library on the deps...
    .library <- function(dep){
        if(!exists(dep)){
            library(dep, character.only = TRUE)
        }
    }
    lapply(deps, .library)
    files <- unlist(lapply(pkgs, function(x){dbfile(get(x))}))
    if (length(not_pkg))
        files <- c(files, rep("", times = length(not_pkg)))
    setNames(files, c(pkgs, not_pkg))
}


## from TxDb
## There are some issues with this as it's currently implemented:
## Because its not using the txdb that is passed in but is instead
## making one up and then assigning it to the global namespace.  This
## is what we want for makeOrgansimDbFromXXX when XXX is UCSC or
## BiomaRt, but its *not* what we want for a simple exsting txdb.
## It's straightforward to do something different when the txdb
## exists().  But: how to I look up it's original name in the
## .GlobalEnv ???  I need to know it's name for the graphData list
## below.
## The other big problem is that if my TxDb objects don't have a
## dbfile() then they can't be saved and re-loaded later.  But this
## problem is not "new" for this object.
## Still TODO: 1) either put the assigned functions in a special env that is
## accesible to these objects OR else make a subclass that can hold those
## objects. 2) make a save method that complains if any of the objects has not
## had saveDb called on it before calling a constructor.  And also add a warning
## to the constructor when somethings dbfile() is an empty string.

### Private environment for storing TxDb objects as needed (failed strategy)
## To work I would have to make the env public.  - It was never a
## great idea fo r this use case.
## TxDbObjs <- new.env(hash=TRUE, parent=emptyenv())
## I may still want to go with the other option of stashing this data
## into a subclass...  For one thing, that option can't have name
## clashes...
## Also: the name this local TxDb gets assigned to cannot be the same as is used
## by a package.  Otherwise a shortened 'custom' TxDb can be overwritten by a
## name clash with a package name...  This could end up being true even if I
## store the TxDb locally inside of a named sub-class.
## Also also: the name should not be made 'special' in the case where
## makeOrganismDbFromTxDb is called as a helper function from within
## makeOrganismDbFromUCSC or makeOrganismDbFromBiomart.

makeOrganismDbFromTxDb <- function(txdb, keytype=NA, orgdb=NA){
    if (!is(txdb, "TxDb"))
        stop("'txdb' must be a TxDb object")
    if(!is(orgdb, "OrgDb") && !is.na(orgdb))
        stop("'orgdb' must be an OrgDb object or NA")
    if (!isSingleStringOrNA(keytype))
        stop("'keytype' must be a single string or NA")

    txdb_store <- new.env(parent = emptyenv())
    ## Then assign that object value to the appropriate name:
    txdbName <- makePackageName(txdb)
    ## We temp assign to global scope
    ## (b/c you need it there if you 'generated' it)
    ## After we can remove it? (will be stored in the object)
    
    assign(txdbName, txdb, envir = txdb_store)
    ## Then get the tax ID:
    taxId <- taxonomyId(txdb_store[[txdbName]])

    ## Then get the name and valued for the OrgDb object
    if (!is(orgdb, "OrgDb") && is.na(orgdb)){
        orgdbName <- OrganismDbi:::.taxIdToOrgDbName(taxId)
        orgdb <- OrganismDbi:::.taxIdToOrgDb(taxId)
        assign(orgdbName, orgdb, envir = txdb_store)
    }else{
        org <- metadata(orgdb)[metadata(orgdb)$name=='ORGANISM',2]
        org <- sub(" ", "_", org)
        orgdbName <- paste0('org.',org,'.db')
        orgdb <- orgdb
        assign(orgdbName, orgdb, envir = txdb_store)
    }
    ## get the primary key for the OrgDb object:
    if(is.na(keytype)){
        geneKeyType <- chooseCentralOrgPkgSymbol(txdb_store[[orgdbName]])
    }else{
        geneKeyType <- keytype
    }

    graphData <- list(join1 = setNames(object=c('GOID', 'GO'),
                                       nm=c('GO.db', orgdbName)),
                      join2 = setNames(object=c(geneKeyType, 'GENEID'),
                                       nm=c(orgdbName, txdbName)))

    ## get the organism
    organism <- organism(txdb_store[[txdbName]])

    #############################################################
    ## Process and then test the graph Data
    gd <- OrganismDbi:::.mungeGraphData(graphData)
    OrganismDbi:::.testGraphData(gd)
    allDeps <- unique(as.vector(gd[,1:2]))
    biocPkgNames <- OrganismDbi:::.biocAnnPackages()
    deps <- allDeps[allDeps %in% biocPkgNames]
    resources <- OrganismDbi:::.gentlyExtractDbFiles(gd, deps, txdb_store)
    ## Check that the fkeys are really columns for the graphData
    fkeys <- OrganismDbi:::.extractPkgsAndCols(gd)
    OrganismDbi:::.testKeys(fkeys, env = txdb_store)
    ## Then make the object:
    graphInfo <- list(graphData=gd, resources=resources)
    OrganismDbi:::OrganismDb(graphInfo=graphInfo)
}


## from UCSC
makeOrganismDbFromUCSC <- function(genome="hg19",
                                   tablename="knownGene",
                                   transcript_ids=NULL,
                                   circ_seqs=NULL,
                                   url="http://genome.ucsc.edu/cgi-bin/",
                               goldenPath.url=getOption("UCSC.goldenPath.url"),
                                   miRBaseBuild=NA){

    ## So call the function to make that TxDb
    txdb <- makeTxDbFromUCSC(genome=genome,
                             tablename=tablename,
                             transcript_ids=transcript_ids,
                             circ_seqs=circ_seqs,
                             url=url,
                             goldenPath.url=goldenPath.url,
                             miRBaseBuild=miRBaseBuild)
    makeOrganismDbFromTxDb(txdb)
}

makeOrganismDbFromBiomart <- function(biomart="ENSEMBL_MART_ENSEMBL",
                                      dataset="hsapiens_gene_ensembl",
                                      transcript_ids=NULL,
                                      circ_seqs=NULL,
                                      filter="",
                                      id_prefix="ensembl_",
                                      host="https://www.ensembl.org",
                                      port,
                                      miRBaseBuild=NA,
                                      keytype="ENSEMBL",
                                      orgdb = NA){

    if (!missing(port))
        warning("The 'port' argument is deprecated and will be ignored.")

    ## So call the function to make that TxDb
    txdb <- makeTxDbFromBiomart(biomart=biomart,
                                dataset=dataset,
                                transcript_ids=transcript_ids,
                                circ_seqs=circ_seqs,
                                filter=filter,
                                id_prefix=id_prefix,
                                host=host,
                                miRBaseBuild=miRBaseBuild)
    makeOrganismDbFromTxDb(txdb, keytype=keytype, orgdb=orgdb)
}

## PROBLEM: OrganismDbi:::.extractDbFiles(gd, deps) requires (strictly) that all objects be available as files somewhere (no exceptions allowed)
## This means that when I get to this stage, with biomaRt, it fails because there is not a TxDb on disc...






















################################################################################
################################################################################


.getMaxEns <- function(srcUrl){
    release <- sub("^ftp://ftp.ensembl.org/pub/release-","",srcUrl)
    release <- unique(sub("/gtf/.*gtf.gz$","", release))
    max(release)
}

## I need a function that will list the GTFs that end users can use to
## make into TxDbs (will probably not overlap perfectly with available
## OrgDbs)
available.GTFsForTxDbs <- function() {
    loadNamespace("AnnotationHub")
    ah <-  AnnotationHub::AnnotationHub()
    ## get OrgDb species
    aho <- subset(ah, ah$rdataclass=='OrgDb')
    oTaxids <-unique(aho$taxonomyid)

    ## get GTF species from ensembl
    ahg <- subset(ah, grepl('gtf.gz$',ah$sourceurl))
    ahg <- subset(ahg, ahg$dataprovider=='Ensembl')
    max <- .getMaxEns(ahg$sourceurl)
    maxStr <- paste0("ftp://ftp.ensembl.org/pub/release-",max)
    ahg <- subset(ahg, grepl(maxStr,ahg$sourceurl))
    gTaxids <-unique(ahg$taxonomyid)

    ## intersect of taxIds
    taxInt <- intersect(oTaxids, gTaxids)
    ## subset down to just the ahgs that can work...
    ahg <- subset(ahg, ahg$taxonomyid %in% taxInt)
    ## for now return the subsetted annotationHubObject
    ahg
}


## And I need a function that will make the transformation for them
makeHubGTFIntoTxDb <- function(ahg){
    if(length(ahg) > 1){
        stop('This function expects only one hub object at a time.')}
    ## get the available GTFs.
    ahgs <- available.GTFsForTxDbs()
    ## Is this one of those? If so, then make it happen
    if(names(ahg) %in%  names(ahgs)){
        txMeta <- data.frame(name='Data source', value='Ensembl GTF')
        txdb <- makeTxDbFromGRanges(ahg[[1]],
                                    metadata= txMeta,
                                    taxonomyId=ahg$taxonomyid)
        require(OrganismDbi)
        ## requires using the 'ENSEMBL' keytype (for these TxDbs)
        ## odb <- makeOrganismDbFromTxDb(txdb, keytype='ENSEMBL')
        odb <- makeOrganismDbFromTxDb(txdb)
    }else{
        stop('No OrgDb information for ', ahg$species)
    }
    odb
}


## testing:
## ahgs <- OrganismDbi:::available.GTFsForTxDbs()
## ahg <-  ahgs[1]
## odb <- OrganismDbi:::makeHubGTFIntoTxDb(ahg)

## debug(OrganismDbi:::makeHubGTFIntoTxDb)
## debug(OrganismDbi:::makeOrganismDbFromTxDb)
## debug(OrganismDbi:::.gentlyExtractDbFiles)


## 1st bad problem is that it basically can't find this:
## org.Ailuropoda_melanoleuca.eg.db.  This is bad since it is one more
## thing that I have to put into the global namespace (temp hack).
## This will need to be done a different way!  Should I use a special
## package environment?  Will it help if I put a slot in that can hold
## an OrdDb?  Or should I find a way to make the resources vector more
## generic by having it only apply to hub items (IOW combine a generic
## character vector with hubcache information). - It is clear that the
## resources vector cannot work for all kinds of paths and situations
## though.  Since people can put objects ANYWHERE.

## 2nd bad problem (even worse if you can believe it), is that these
## orgDbs don't have ensembl IDs in them.  This is bad because that
## basically means that none of these organisms can be supported for
## these TxDbs (which will use ensembl based gene identifiers.  I can
## update the OrgDbs, it's just a crucial missing piece that has to be
## retrieved before I can proceed.

## 3rd problem this uncovered is that I had to use AH cars for the
## object just to get easy access to the metadata.  It seems like we
## could benefit from a universal application of that metadata to the
## contents that come out of these cars...  Then it would be
## straightforward to write a method like I wanted to here.


