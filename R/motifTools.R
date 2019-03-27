#' PWM hits and motif scores as an RLElist
#'
#' Creates rlelist of pwm hits. 
#' 
#'
#'
#' @docType methods
#' @name pwmToCoverage
#' @rdname pwmToCoverage
#' 
#' @author Thomas Carroll
#'
#' @param pwm A PWM matrix object.
#' @param genome A BSgenome object
#' @param min pwm score (as percentage of maximum score) cutoff
#' @param removeRand Remove contigs with rand string
#' @param chrsOfInterest Chromosomes to use
#' @return  A RLElist of motif density per base pair to be used as input to main soggi function.
#' @examples 
#' data(pwmCov)
#' data(singleGRange)
#' 
#' 
#' @export
pwmToCoverage <- function(pwm,genome,min="70%",removeRand=FALSE,chrsOfInterest=NULL){
  
  allchrs <- seqnames(genome)
  if(!is.null(allchrs)){
    allchrs <- allchrs[allchrs %in% chrsOfInterest]
  }
  if(removeRand){
    allchrs <- allchrs[!grepl("random",allchrs,ignore.case=TRUE)]
  }
  intergerList <- lapply(allchrs,function(x)pwmHitAsCoverage(pwm,genome,min,x))
  myrle <- RleList(intergerList,compress=FALSE)
  names(myrle) <- allchrs
  myrle
}  

pwmHitAsCoverage <- function(pwm,genome,min,chrofinterest){
  posMotifs <- matchPWM(pwm,genome[[chrofinterest]],min.score=min)
  negMotifs <- matchPWM(reverseComplement(pwm),genome[[chrofinterest]],min.score=min)
  if(length(posMotifs) > 0){
    rleMotifHitPos <- coverage(GRanges(seqnames=chrofinterest,ranges(posMotifs),strand="+"),width=length(genome[[chrofinterest]]))
  }else{
    rleMotifHitPos <- RleList(rep(0,length(genome[[chrofinterest]])))
  }
  if(length(negMotifs) > 0){
    rleMotifHitNeg <- coverage(GRanges(seqnames=chrofinterest,ranges(negMotifs),strand="-"),width=length(genome[[chrofinterest]]))
  }else{
    rleMotifHitNeg <- RleList(rep(0,length(genome[[chrofinterest]])))    
  }
  rleTotal <- rleMotifHitPos+rleMotifHitNeg
  return(rleTotal[[1]])
}

pwmHitAsGRanges <- function(pwm,genome,min,chrofinterest){
  posMotifs <- matchPWM(pwm,genome[[chrofinterest]],min.score=min,with.score=TRUE)
  negMotifs <- matchPWM(reverseComplement(pwm),genome[[chrofinterest]],min.score=min,with.score=TRUE)
  posMotifs <- GRanges(seqnames=rep(chrofinterest,length(posMotifs)),
                       ranges=ranges(posMotifs),
                       strand=rep("+",length(posMotifs))                      
                       )
  negMotifs <- GRanges(seqnames=rep(chrofinterest,length(negMotifs)),
                       ranges=ranges(negMotifs),
                       strand=rep("-",length(negMotifs))                      
                        )
  #strand(negMotifs) <- rep("-",length(negMotifs))
  GRangesTotal <- c(posMotifs,negMotifs)
  return(GRangesTotal)
}

pwmToGranges <- function(pwm,genome,min,chrs=NULL){
  if(is.null(chrs)){
    chrs <- seqnames(genome)
  }
  chrs <- unique(chrs)
  Res <- lapply(chrs,function(x)pwmHitAsGRanges(pwm,genome,min,x))
  pwmHitsAsGRanges <- unlist(GRangesList(unlist(Res)))
}

makeGRangesWithSummary <- function(GRangesSet,scoreBy){
  temp <- viewSums(Views(scoreBy,
                 ranges(GRangesSet)
            )
  )  
  mcols(GRangesSet) <- data.frame(score=temp)
  return(GRangesSet)
}
  

rleFromScoresInGRanges <- function(GRangesSet,scoreBy,chrs){
  chrs <- chrs[chrs %in% names(scoreBy)]
  chrs <- chrs[chrs %in% names(seqlengths(GRangesSet))]  
  res <- lapply(chrs,function(x)
                makeGRangesWithSummary(
                  GRangesSet[
                    seqnames(GRangesSet) %in% x,],
                  scoreBy[[x]]        
                  ))
  return(unlist(GRangesList(unlist(res))))
}

motifCov <- function(genome,regions,pwm,chrOfInterest,atCentre=FALSE){
  reducedregions <- reduce(regions[seqnames(regions) %in% chrOfInterest])
  regionViews <- Views(genome[[chrOfInterest]],ranges(reducedregions))
  trial <- matchPWM(pwm,regionViews,min.score = 0,with.score = TRUE)
  if(atCentre==TRUE){
    theRanges <- resize(as(trial,"IRanges"),1,"center")
  }
  if(atCentre==FALSE){
    theRanges <- as(trial,"IRanges")
  }
  if(length(theRanges) > 0){
    motifCov <- unlist(coverage(GRanges(chrOfInterest,theRanges,"*",mcols(trial)$score),weight="mcols.trial..score"))
  }else{
    motifCov <- unlist(RleList(rep(0,length(genome[[chrOfInterest]]))))
  }
  return(motifCov)
}

#' Motif score as an RLElist
#'
#' @name makeMotifScoreRle
#' @rdname pwmToCoverage
#' 
#'
#' @param regions GRanges object to include in pwm rlelist
#' @param extend bps to extend regions by
#' @param strandScore Method for averaging strand. Options are max, mean, sum, bothstrands
#' @param atCentre TRUE/FALSE. TRUE assigns score onto 1bp position at centre of motif.
#' FALSE assigns every basepair the sum of scores of all overlapping motifs. 
#' @export
makeMotifScoreRle <- function(pwm,regions,genome,extend,removeRand=FALSE,strandScore="mean",atCentre=FALSE){
  regions <- GRanges(seqnames(regions),IRanges(start(regions)-extend,end(regions)+extend),strand=strand(regions),mcols(regions))
  lengths <- seqlengths(genome)
  ## Filter testRanges to those contained within chromosomes.
  message("Filtering regions which extend outside of genome boundaries...",appendLF = FALSE)
  testRangeNames <- unique(seqnames(regions))
  temptestranges <- GRanges()
  maxDistance <- extend
  for(i in 1:length(testRangeNames)){
    perchrRanges <- regions[seqnames(regions) %in% as.vector(testRangeNames[i])]
    temptestranges <- c(temptestranges,perchrRanges[end(perchrRanges)+maxDistance < lengths[names(lengths) %in% testRangeNames[i]]
                                                    & start(perchrRanges)-maxDistance > 0 ])
    #print(i)
  }
  
  message("..Done")
  message("Filtered ",length(regions)-length(temptestranges)," of ",length(regions)," regions")
  
  regions <- temptestranges
  
  allchrs <- as.vector(unique(seqnames(regions)))
  
  if(removeRand){
    allchrs <- allchrs[!grepl("random",allchrs,ignore.case=TRUE)]
  }
  forMotif <- list()
  revMotif <- list()
  message("Scoring motifs on positive strand...",appendLF = FALSE)  
  #motifScoreRLE <- lapply(allchrs,function(x)motifCov(genome,regions,pwm,x))
  for(k in 1:length(allchrs)){
    message(allchrs[k])
    forMotif[[k]] <- unlist(motifCov(genome,regions,pwm,allchrs[k],atCentre))  
  }
  message("..done")
  
  message("Scoring motifs on negative strand...",appendLF = FALSE)  
  #motifScoreRLE <- lapply(allchrs,function(x)motifCov(genome,regions,pwm,x))
  for(k in 1:length(allchrs)){
    message(allchrs[k])    
    revMotif[[k]] <- unlist(motifCov(genome,regions,reverseComplement(pwm),allchrs[k]))    
  }
  message("..done")
  
  #motifScoreRLEComplement <- lapply(allchrs,function(x)motifCov(genome,regions,reverseComplement(pwm),x))
  #myrle <- RleList(motifScoreRLE,compress=F)
  # myrleComplement <- RleList(motifScoreRLEComplement,compress=F)
  revMotif <- RleList(revMotif,compress=FALSE)
  forMotif <- RleList(forMotif,compress=FALSE)
  if(strandScore=="sum"){
    MotifScore <- revMotif+forMotif
    names(MotifScore) <- allchrs    
  }
  if(strandScore=="mean"){  
    MotifScore <- (revMotif+forMotif)/2
    names(MotifScore) <- allchrs    
  }
  if(strandScore=="max"){  
    MotifScore <- pmax(revMotif,forMotif)
    names(MotifScore) <- allchrs
  }  
  if(strandScore=="bothstrands"){  
    MotifScore <- list(revMotif,forMotif)
    names(MotifScore[[1]]) <- allchrs
    names(MotifScore[[2]]) <- allchrs    
  }  
  return(MotifScore)
}


plotCuts <- function(BAM,PWMlist,Genome,selection=c(0,300),chromosomes=NULL,distanceAround=500,cutoff=NULL,ntop=20000,verbose=TRUE,ext="",name=NULL){
  # require(GenomicAlignments)
  # require(ggplot2)
  if(is.null(cutoff)) cutoff <- rep("80%",length(PWMlist))
  
  
  if(file.exists(BAM) & is.na(index(BamFile(BAM)))){
    message("Creating index for ",BAM)
    indexBam(BAM)
    message("..done")
  }
  
  seqTableDF <- seqinfo(BamFile(BAM))
  seqTableGR <- GRanges(seqTableDF)
  #Genome <- Genome[seqnames(Genome) %in% seqnames(seqTableGR)]
  if(!is.null(chromosomes)){
    seqTableGR <- seqTableGR[seqnames(seqTableGR) %in% chromosomes]
    #Genome <- Genome[seqnames(Genome) %in% chromosomes]
    chromosomes <- chromosomes[chromosomes %in% seqnames(seqTableGR)]
  }else{
    chromosomes <- unique(seqnames(seqTableGR))
  }
  if(verbose) message("Reading BAM file..",appendLF=FALSE)
  atacReads_Open <- readGAlignmentPairs(BAM,
                                        param=ScanBamParam(which = seqTableGR,
                                                           what=c("qname","mapq","isize")))
  if(verbose) message("..done",appendLF=TRUE)
  total <- length(atacReads_Open)
  if(verbose) message("Read ",total," fragments",appendLF=TRUE)              
  if(verbose) message("Filtering fragments by insert size..",appendLF=FALSE)
  insertSizes <- abs(elementMetadata(GenomicAlignments::first(atacReads_Open))$isize)
  atacReads_Open <- atacReads_Open[insertSizes > selection[1] & insertSizes < selection[2]]
  totalLeft <- length(atacReads_Open)
  if(verbose) message("..done",appendLF=TRUE)
  if(verbose) message("Filter ",total-totalLeft," fragments outside insert size range of ",selection[1]," to ",selection[2],appendLF=TRUE)   
  if(verbose) message("After filtering, ",totalLeft," fragments",appendLF=TRUE)       
  if(verbose) message("Creating cuts..",appendLF=FALSE)
  
  read1 <- GenomicAlignments::first(atacReads_Open)
  read2 <- GenomicAlignments::second(atacReads_Open)
  
  Firsts <- resize(granges(read1),fix="start",1)
  First_Pos_toCut <- shift(granges(Firsts[strand(read1) == "+"]),
                           4)
  First_Neg_toCut <- shift(granges(Firsts[strand(read1) == "-"]),
                           -5)
  
  Seconds <- resize(granges(read2),fix="start",1)
  Second_Pos_toCut <- shift(granges(Seconds[strand(read2) == "+"]),
                            4)
  Second_Neg_toCut <- shift(granges(Seconds[strand(read2) == "-"]),
                            -5)
  
  test_toCut <- c(First_Pos_toCut,First_Neg_toCut,
                  Second_Pos_toCut,Second_Neg_toCut)
  if(verbose) message("..done",appendLF=TRUE)
  if(verbose) message("Creating coverage from cuts..",appendLF=FALSE)
  cutsCoverage <- coverage(test_toCut)
  if(verbose) message("..done",appendLF=TRUE)
  if(verbose) message("Scannning for motifs",appendLF=TRUE)
  emptyMGR <- GRanges()
  for(i in 1:length(PWMlist)){
    PWM <- PWMlist[[i]]
    if(verbose) message("Checking positive strand..",appendLF=FALSE)
    myAllMatch <- lapply(chromosomes,
                         function(x)matchPWM(PWM,Genome[[x]],min.score = cutoff[i],with.score = TRUE))
    names(myAllMatch) <- chromosomes
    myAllMatchPos <- unlist(GRangesList(lapply(names(myAllMatch),
                                               function(x)GRanges(x,ranges(myAllMatch[[x]]),score=mcols(myAllMatch[[x]])$score))
    ))
    if(verbose) message("..done",appendLF=TRUE)
    if(verbose) message("Found ",length(myAllMatchPos)," ",names(PWMlist)[i]," motifs on positive strand",appendLF=TRUE)       
    
    if(verbose) message("Checking negative strand..",appendLF=FALSE)
    
    myAllMatch <- lapply(chromosomes,
                         function(x)matchPWM(reverseComplement(PWM),Genome[[x]],min.score = cutoff[i],with.score = TRUE))
    names(myAllMatch) <- chromosomes
    myAllMatchNeg <- unlist(GRangesList(lapply(names(myAllMatch),
                                               function(x)GRanges(x,ranges(myAllMatch[[x]]),score=mcols(myAllMatch[[x]])$score))
    ))
    if(verbose) message("..done",appendLF=TRUE)
    if(verbose) message("Found ",length(myAllMatchNeg)," ",names(PWMlist)[i]," motifs on negative strand",appendLF=TRUE)       
    
    myAllMatch <- c(myAllMatchPos,myAllMatchNeg)
    myAllMatch$Motif <- names(PWMlist)[i]
    if(verbose) message("Found ",length(myAllMatch)," ",names(PWMlist)[i]," total motifs",appendLF=TRUE) 
    myAllMatch <- myAllMatch[order(myAllMatch$score,decreasing = TRUE),]
    message(max(ntop,length(myAllMatch)))
    # myAllMatch <- myAllMatch[1:max(ntop,length(myAllMatch)),]
    if(verbose) message("Using ",min(ntop,length(myAllMatch))," of total motifs",appendLF=TRUE) 
    emptyMGR <- c(emptyMGR,myAllMatch)
  }
  suppressPackageStartupMessages(require(soGGi))
  if(verbose) message("Creating cut profile around motifs..",appendLF=FALSE)       
  
  Motif_Cuts <- regionPlot(cutsCoverage,
                           testRanges = emptyMGR,
                           style = "point",
                           format="rlelist",distanceAround = distanceAround,verbose=verbose)
  if(verbose) message("..done",appendLF=TRUE)
  myP <- plotRegion(Motif_Cuts,summariseBy="Motif",freeScale=TRUE)
  if(is.null(name)){
    outName <- paste0(dirname(BAM),
                      paste0(gsub("\\.bam","",basename(BAM)),"_",ext,"_MotifCuts.png"))
  }else{
    outName <- name
  }
  ggsave(myP,
         file=outName)
  return(Motif_Cuts)
}

