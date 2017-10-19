
#################################################################################################################
# Big-LD <input>
#' @title Estimation of LD block regions
#' @name Big_LD
#'
#' @description \code{Big_LD} returns the estimation of LD block regions of given data.
#'
#' @param geno A data frame or matrix of additive genotype data, each column is additive genotype of each SNP.
#' @param SNPinfo A data frame or matrix of SNPs information.  1st column is rsID and 2nd column is bp position.
#' @param CLQcut A numeric value of threshold for the correlation value |r|, between 0 to 1.
#' @param clstgap  An integer value to specifying the threshold of physical distance (bp) between two consecutive SNPs
#' which do not belong to the same clique, i.e., if a physical distance between two consecutive SNPs in a clique
#' greater than \code{clstgap}, then the algorithm split the cliques satisfying each
#' clique do not contain such consecutive SNPs
#' @param leng  An integer value to specify the number of SNPs in a preceding and a following region
#' of each sub-region boundary, every SNP in a preceding and every SNP in a following region need to be in weak LD.
#' @param MAFcut An numeric value to specifying the MAF threshold. 
#' @param subSegmSize  An integer value to specify the upper bound of the number of SNPs in a one-take sub-region.
#' @param appendRare If \code{appendRare = TRUE}, the algorithm append rare SNPs (MAF<MAFcut) to the constructed LD blocks or add a new LD blocks
#' @param checkLargest If \code{checkLargest = TRUE}, the algorithm use heuristic procedure to reduce runtime of CLQ-D execution
#'
# <output>
#' @return  A data frame of block estimation result.
#' Each row of data frame shows the starting SNP and end SNP of each estimated LD block.
#' 
#' @author Sun-Ah Kim <sunny03@snu.ac.kr>, Yun Joo Yoo <yyoo@snu.ac.kr>
#' @seealso \code{\link{CLQD}}, \code{\link{LDblockHeatmap}}
#'
#' @examples
#'
#' data(geno)
#' data(SNPinfo)
#' Big_LD(geno, SNPinfo)
#' Big_LD(geno, SNPinfo, CLQcut = 0.5, clstgap = 40000, leng = 200, subSegmSize = 1500)
#'
# sub-Functions 1. CLQD < built - in > 2. cutsequence.modi, 3.intervalCliqueList, 4. find.maximum.indept, 5. constructLDblock,
#' @export
#' 
Big_LD <- function(geno, SNPinfo, CLQcut = 0.5, clstgap = 40000, leng = 200, subSegmSize = 1500, MAFcut = 0.05, 
                   appendrare = TRUE, checkLargest = FALSE) {
  # packages
  # library(igraph)
  #######################################################################################################
  # sub-Functions 1. cutsequence.modi, 2.intervalCliqueList, 3. find.maximum.indept, 4. constructLDblock, 5. CLQ
  
  cutsequence.modi = function(geno, leng, subSegmSize) {
    print("split whole sequence into subsegments")
    modeNum <- 1
    lastnum <- 0 
    ## region length<=3000
    if (dim(geno)[2] <= subSegmSize) {
      print("cutting sequence, done")
      print("there is only one sub-region!")
      return(list(dim(geno)[2], NULL))
    } else {
      cutpoints <- NULL
      i = leng
      while (i <= (dim(geno)[2] - leng)) {
        if((i-lastnum) > 5*subSegmSize){
          modeNum <- 2
          break; 
        }
        # tick size = leng * (1/10)
        tick <- as.integer(leng * 1/10)
        nowcm <- cor(geno[, (i - tick + 1):(i + tick)],  use="pairwise.complete.obs")
        nowr2 <- nowcm^2
        nowr2[which(nowr2 < 0.5)] = 0
        diag(nowr2) <- 0
        if (length(which(nowr2[1:tick, (tick + 1):(2 * tick)] > 0)) > 0) {
          # print(i)
          i = i + 1
          next
        }
        # tick size = leng
        tick <- leng
        nowcm <- cor(geno[, (i - tick + 1):(i + tick)],  use="pairwise.complete.obs")
        nowr2 <- nowcm^2
        nowr2[which(nowr2 < 0.5)] = 0
        diag(nowr2) <- 0
        if (length(which(nowr2[1:tick, (tick + 1):(2 * tick)] > 0)) > 0) {
          i = i + 1
          next
        } else {
          print(c(i))
          cutpoints <- c(cutpoints, i)
          i <- i + (leng/2)
          lastnum = i
        }
      }##end while
      if(modeNum == 1){
        cutpoints <- c(cutpoints, dim(geno)[2])
        # separate too big regions candi.cutpoints return(cutpoints,candi.cutpoints)
        atfcut <- NULL
        while (max(diff(cutpoints)) > subSegmSize) {
          diffseq = diff(cutpoints)
          recutpoint <- which(diffseq > subSegmSize)
          nowmaxsize = max(diff(cutpoints))
          # print(nowmaxsize) print(recutpoint)
          tt <- cbind((cutpoints[recutpoint] + 1), cutpoints[recutpoint + 1])
          numvec = NULL
          for (i in 1:dim(tt)[1]) {
            st <- tt[i, 1]
            ed <- tt[i, 2]
            if (ed > (dim(geno)[2] - leng)) {
              ed <- dim(geno)[2] - leng
            }
            weakcount <- sapply(c((st + leng):(ed - leng)), function(x) {
              tick <- as.integer(leng/5)
              nowCM <- cor(geno[, (x - tick + 1):(x + tick)],  use="pairwise.complete.obs")
              nowr2 <- nowCM^2
              diag(nowr2) <- 0
              length(which(nowr2[(1:tick), (tick + 1):(2 * tick)] > 0.5))
            })
            weakcount.s <- sort(weakcount)
            weaks <- weakcount.s[10]
            weakpoint <- which(weakcount <= weaks)
            weakpoint <- weakpoint + st + leng - 1
            nearcenter = sapply(weakpoint, function(x) abs((ed - x) - (x - st)), simplify = TRUE)
            addcut <- weakpoint[which(nearcenter == min(nearcenter))][1]
            numvec <- c(numvec, addcut)
            atfcut <- c(atfcut, addcut)
          }  ##end for
          cutpoints <- sort(c(cutpoints, numvec))
          newcandi = which(diff(cutpoints) > subSegmSize)
          remaxsize = max(diff(cutpoints))
          # print(remaxsize) print(newcandi)
          if (length(newcandi) == 0) {
            break
          }
        }
      }
      ##end while
      if(modeNum == 2){
        cutpoints = seq(subSegmSize, dim(geno)[2], subSegmSize/2)
        atfcut = cutpoints
        if(max(cutpoints) == dim(geno)[2]){
          atfcut = atfcut[-(length(atfcut))]
        }else{
          cutpoints = c(cutpoints, dim(geno)[2])
        }
        
      }
    }
    print("cutting sequence, done")
    return(list(cutpoints, atfcut))
  }
  intervalCliqueList = function(clstlist, allsnps, onlybp) {
    bp.clstlist <- lapply(clstlist, function(x) onlybp[x])  ###
    bp.allsnps <- lapply(allsnps, function(x) onlybp[x])
    
    IMsize <- length(bp.clstlist)  ## adjacency matrix of intervals in interval graph
    adjacencyM <- matrix(0, IMsize, IMsize)
    for (i in 1:IMsize) {
      for (j in 1:IMsize) {
        adjacencyM[i, j] <- length(intersect(bp.allsnps[[i]], bp.allsnps[[j]]))
      }
    }
    diag(adjacencyM) <- 0
    interval.graph <- graph.adjacency(adjacencyM, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL)
    # print(paste("max coreness", max(coreness(interval.graph))))
    # print(paste("ecount", ecount(interval.graph), "vertex*5 ", 5*IMsize))
    if(max(coreness(interval.graph))>10){ #ecount(interval.graph)> 5*IMsize|
      interval.cliques <- maximal.cliques(interval.graph, min = 1)  
    }else{
      interval.cliques <- cliques(interval.graph, min = 1)  
    }
    interval.cliques <- interval.cliques[order(sapply(interval.cliques, min))]
    
    intervals <- lapply(interval.cliques, function(x) unlist(bp.clstlist[x]))
    intervals <- lapply(intervals, sort)
    intervals <- lapply(intervals, unique)
    weight.itv <- sapply(intervals, length)
    
    intervals.range <- t(sapply(intervals, range))
    unique.intervals.range <- unique(intervals.range)
    
    rangeinfo <- cbind(intervals.range, weight.itv)
    
    interval.info <- apply(unique.intervals.range, 1, function(x) {
      info1 = which(rangeinfo[, 1] == x[1])
      info2 = which(rangeinfo[, 2] == x[2])
      sameitv <- intersect(info1, info2)
      maxweight <- max(rangeinfo[sameitv, 3])
      sameitv[which(maxweight == rangeinfo[sameitv, 3])][1]
    })
    
    final.intervals <- intervals[interval.info]
    final.intervals.w <- rangeinfo[interval.info, 3]
    return(list(final.intervals, final.intervals.w))
  }
  # find maximum weight independent set input: clique interval list, clique weights
  find.maximum.indept = function(sample.itv, sample.weight) {
    n.of.sample <- length(sample.itv)
    interval.range <- t(sapply(sample.itv, range))
    pre.range <- as.list(rep(NA, n.of.sample))  #pre.range : range of predecessor
    ## pre.range : n by 2, i row (x,y) : possible predecessors of i interval are from x interval to y interval
    for (i in 1:n.of.sample) {
      nowstart <- interval.range[i, 1]
      if (length(which(interval.range[, 2] < nowstart)) > 0) {
        pre.range[[i]] <- which(interval.range[, 2] < nowstart)
      }
    }
    sources <- c(1:n.of.sample)[(sapply(pre.range, function(x) all(is.na(x)) == TRUE))]
    
    ## source of comparability graph of complement of Interval graph
    if (length(sources) < n.of.sample) {
      not.s <- setdiff(c(1:n.of.sample), sources)
      for (i in not.s) {
        pre.pre <- sort(unique(unlist(pre.range[pre.range[[i]]])))
        pre.range[[i]] <- setdiff(pre.range[[i]], pre.pre)
      }
      names(pre.range) <- sample.weight
      n.interval <- c(1:n.of.sample)
      route.weights <- rep(0, n.of.sample)  ##cumulative weights
      route.weights[sources] <- sample.weight[sources]
      pointers <- rep(0, n.of.sample)  ## predecessor of current interval
      pointers[sources] <- NA
      explored <- rep(0, n.of.sample)
      explored[sources] <- 1
      info <- cbind(n.interval, route.weights, pointers, explored)
      
      for (i in not.s) {
        maybe.pred <- pre.range[[i]]
        now.info <- info[maybe.pred, , drop = FALSE]
        max.info <- now.info[which(now.info[, 2] == max(now.info[, 2])), , drop = FALSE]
        if (dim(max.info)[1] > 1)
          max.info <- max.info[1, , drop = FALSE]
        info[i, 2] <- sample.weight[i] + max.info[2]
        info[i, 3] <- max.info[1]
        info[i, 4] <- 1
      }
      
      #### trace maximum independent set
      start.itv <- which(info[, 2] == max(info[, 2]))[1]
      predecessor <- info[start.itv, 3]
      indept.set <- c(predecessor, start.itv)
      while (!is.na(predecessor)) {
        predecessor <- info[predecessor, 3]
        indept.set <- c(predecessor, indept.set)
      }
      
      indept.set <- as.vector(indept.set)
      indept.set <- indept.set[-which(is.na(indept.set))]
      indept.set.weight <- max(info[, 2])
    } else {
      indept.set = which(sample.weight == max(sample.weight))
      indept.set.weight = max(sample.weight)
    }
    
    final.result <- list(indept.set = indept.set, indept.set.weight = indept.set.weight)
    return(final.result)
  }
  constructLDblock = function(clstlist, subSNPinfo) {
    # subfunction: intervalCliqueList, find.maximum.indept
    Totalblocks = NULL
    while (length(clstlist) > 0) {
      allsnps <- lapply(clstlist, function(x) c(min(x):max(x)))
      onlybp <- subSNPinfo[, 2]
      candi.interval <- intervalCliqueList(clstlist, allsnps, onlybp)
      intervals <- candi.interval[[1]]  ## list of SNPs in each cliques
      weight.itv <- candi.interval[[2]]  ## weights of each cliques
      MWIS <- find.maximum.indept(intervals, weight.itv)  ##find independent set
      indept.set <- intervals[MWIS[[1]]]
      LDintervals <- lapply(indept.set, function(x) match(x, subSNPinfo[, 2]))
      subLDblocks <- t(sapply(LDintervals, range))
      Totalblocks <- rbind(Totalblocks, subLDblocks)
      takenSNPs <- apply(Totalblocks, 1, function(x) c(min(x):max(x)))
      takenSNPs <- as.vector(unlist(takenSNPs))
      clstlist <- lapply(clstlist, function(x) setdiff(x, takenSNPs))
      clstlist <- clstlist[sapply(clstlist, function(x) length(x) > 1)]
      if (length(clstlist) == 0) break
      addinglist <- NULL
      for (n in 1:length(clstlist)) {
        nowbin <- clstlist[[n]]
        intersection <- intersect(c(min(nowbin):max(nowbin)), takenSNPs)
        if (length(intersection) > 0) {
          part1 <- nowbin[which(nowbin < min(intersection))]
          part2 <- setdiff(nowbin, c(min(part1):max(intersection)))
          clstlist[[n]] <- part1
          addinglist <- c(addinglist, list(part2))
        }
      }
      clstlist <- c(clstlist, addinglist)
      clstlist <- clstlist[sapply(clstlist, function(x) length(x) > 1)]
    }
    return(Totalblocks)
  }
  subBigLD = function(subgeno, subSNPinfo,  CLQcut, clstgap, checkLargest){
    subbinvec <- CLQD(subgeno, subSNPinfo, CLQcut, clstgap, CLQmode = "Density",codechange = FALSE, checkLargest)
    # print('CLQ done!')
    bins <- c(1:max(subbinvec[which(!is.na(subbinvec))]))
    clstlist <- sapply(bins, function(x) which(subbinvec == x), simplify = FALSE)
    clstlist <- lapply(clstlist, sort)  ###
    clstlist <- clstlist[order(sapply(clstlist, min))]  ###
    nowLDblocks <- constructLDblock(clstlist, subSNPinfo)
    # print('constructLDblock done!')
    nowLDblocks <- nowLDblocks[order(nowLDblocks[, 1]), , drop = FALSE]
    return(nowLDblocks)
    
  }
  appendSGTs = function(LDblocks, Ogeno, OSNPinfo, CLQcut, clstgap, checkLargest){
    expandB = NULL
    # failB = NULL
    snp1 = which(OSNPinfo[,2]<LDblocks[1,5])
    if(length(snp1)>2){
      OSNPs = 1:max(snp1)
      firstB = LDblocks[1,]  
      secondSNPs = which(OSNPinfo[,2]>=firstB$start.bp & OSNPinfo[,2] <= firstB$end.bp)
      cor2 = suppressWarnings(cor(Ogeno[,c(secondSNPs, OSNPs),drop=FALSE],  use="pairwise.complete.obs")^2)
      cor2 = cor2[1:length(secondSNPs), -(1:length(secondSNPs)), drop=FALSE]
      cor2num = apply(cor2, 2, function(x) {
        sum(x>CLQcut^2)
      })
      cor2ratio = cor2num/(dim(cor2)[1])
      # grid.draw(LDr2map(genoDp(Ogeno[,min(firstSNPs):max(secondSNPs)]),
      # c(0, length(firstSNPs), length(firstSNPs)+length(OSNPs),length(firstSNPs)+length(OSNPs)+length(secondSNPs)),1))
      cor2numT = cor2ratio>0.6
      cor2numT = c(cor2numT, 1)
      points2 = min(which(cor2numT>0))
      NsecondSNPs = points2:max(secondSNPs)
      reOSNPs = setdiff(c(1:max(NsecondSNPs)), NsecondSNPs)
      if(length(reOSNPs)>1){
        subgeno = Ogeno[, reOSNPs]
        subSNPinfo = OSNPinfo[reOSNPs,]
        subBlocks = subBigLD(subgeno, subSNPinfo,  CLQcut, clstgap, checkLargest)
        subBlocks = subBlocks+min(reOSNPs)-1
        expandB = rbind(expandB, subBlocks)
      }
      firstSNPs=NsecondSNPs
    }else{
      firstB = LDblocks[1,]  
      firstSNPs = which(OSNPinfo[,2]>=firstB$start.bp & OSNPinfo[,2] <= firstB$end.bp)
    }
    if(dim(LDblocks)[1]>1){
      for(i in 1:(dim(LDblocks)[1]-1)){
        secondB = LDblocks[(i+1),]
        secondSNPs = which(OSNPinfo[,2]>=secondB$start.bp & OSNPinfo[,2]<= secondB$end.bp)
        OSNPs = setdiff(max(firstSNPs):min(secondSNPs), c(max(firstSNPs), min(secondSNPs)))
        if(length(OSNPs)==0){
          expandB = rbind(expandB, range(firstSNPs))
          firstSNPs = secondSNPs
        }else{
          cor1 = suppressWarnings(cor(Ogeno[,c(firstSNPs, OSNPs),drop=FALSE],  use="pairwise.complete.obs")^2)
          cor1 = cor1[1:length(firstSNPs), -(1:length(firstSNPs)), drop=FALSE]
          cor1num = apply(cor1, 2, function(x) {
            sum(x>CLQcut^2)
          })
          cor1ratio = cor1num/(dim(cor1)[1])
          # cor1num = apply(cor1r2, 2, function(x) length(which(x>CLQcut^2))/length(firstSNPs))
          cor2 = suppressWarnings(cor(Ogeno[,c(secondSNPs, OSNPs),drop=FALSE],  use="pairwise.complete.obs")^2)
          cor2 = cor2[1:length(secondSNPs), -(1:length(secondSNPs)), drop=FALSE]
          cor2num = apply(cor2, 2, function(x) {
            sum(x>CLQcut^2)
          })
          cor2ratio = cor2num/(dim(cor2)[1])
          # grid.draw(LDr2map(cor(Ogeno[,min(firstSNPs):max(secondSNPs)]),
          # c(0, length(firstSNPs), length(firstSNPs)+length(OSNPs),length(firstSNPs)+length(OSNPs)+length(secondSNPs)),1))
          cor1numT = cor1ratio>0.6
          cor2numT = cor2ratio>0.6
          cor1numT = c(1, cor1numT, 0)
          cor2numT = c(0, cor2numT, 1)
          points1 = max(firstSNPs)+max(which(cor1numT>0))-1
          NfirstSNPs = min(firstSNPs):points1
          points2 = max(firstSNPs)+max(which(cor2numT>0))-1
          NsecondSNPs = points2:max(secondSNPs)
          if(max(NfirstSNPs)<min(NsecondSNPs)){
            expandB = rbind(expandB, range(NfirstSNPs))
            reOSNPs = setdiff(c(min(NfirstSNPs):max(NsecondSNPs)), c(NfirstSNPs, NsecondSNPs))
            if(length(reOSNPs)>1){
              subgeno = Ogeno[, reOSNPs]
              subSNPinfo = OSNPinfo[reOSNPs,]
              subBlocks = subBigLD(subgeno, subSNPinfo,  CLQcut, clstgap, checkLargest)
              subBlocks = subBlocks+min(reOSNPs)-1
              expandB = rbind(expandB, subBlocks)
            }
            firstSNPs=NsecondSNPs
          }else{
            #merge two blocks
            subgeno = Ogeno[, c(min(firstSNPs):max(secondSNPs))]
            subSNPinfo = OSNPinfo[c(min(firstSNPs):max(secondSNPs)),]
            subBlocks = subBigLD(subgeno, subSNPinfo,  CLQcut, clstgap, checkLargest)
            subBlocks = subBlocks+min(firstSNPs)-1
            if(dim(subBlocks)[1]==1) {
              firstSNPs = subBlocks[1,1]:subBlocks[1,2]
            }else{
              expandB = rbind(expandB, subBlocks[-(dim(subBlocks)[1]),]) 
              firstSNPs = subBlocks[(dim(subBlocks)[1]),1]:subBlocks[(dim(subBlocks)[1]),2]
            }
          }
          print(c(i, dim(LDblocks)[1]))
          # print(tail(expandB))
          # if(i >= 30) break
        }
      }
    }
    # firstSNPs
    if(max(firstSNPs)<(dim(Ogeno)[2]-1)){
      OSNPs = (max(firstSNPs)+1):(dim(Ogeno)[2])
      cor1 = suppressWarnings(cor(Ogeno[,c(firstSNPs, OSNPs),drop=FALSE],  use="pairwise.complete.obs")^2)
      cor1 = cor1[1:length(firstSNPs), -(1:length(firstSNPs)), drop=FALSE]
      cor1num = apply(cor1, 2, function(x) {
        sum(x>CLQcut^2)
      })
      cor1ratio = cor1num/(dim(cor1)[1])
      cor1numT = cor1ratio>0.6
      cor1numT = c(1, cor1numT, 0)
      points1 = max(firstSNPs)+max(which(cor1numT>0))-1
      NfirstSNPs = min(firstSNPs):points1
      expandB = rbind(expandB, range(NfirstSNPs))
      reOSNPs = setdiff(c(min(NfirstSNPs):dim(Ogeno)[2]), c(NfirstSNPs))
      if(length(reOSNPs)>1){
        subgeno = Ogeno[, reOSNPs]
        subSNPinfo = OSNPinfo[reOSNPs,]
        subBlocks = subBigLD(subgeno, subSNPinfo,  CLQcut, clstgap, checkLargest)
        subBlocks = subBlocks+min(reOSNPs)-1
        expandB = rbind(expandB, subBlocks)
      }
    }else{
      expandB = rbind(expandB, range(firstSNPs))
    }
    # LDblocks = expandB 
    expandB = expandB[(expandB[,1]!=expandB[,2]),,drop=FALSE]
    start.bp <- OSNPinfo[, 2][expandB[, 1]]
    end.bp <- OSNPinfo[, 2][expandB[, 2]]
    start.rsID <- as.character(OSNPinfo[, 1][expandB[, 1]])
    end.rsID <- as.character(OSNPinfo[, 1][expandB[, 2]])
    TexpandB <- data.frame(expandB, start.rsID, end.rsID, start.bp, end.bp)
    colnames(TexpandB) <- c("start", "end", "start.rsID", "end.rsID", "start.bp", "end.bp")
    return(TexpandB)
  }
  #######################################################################################################
  # Main part input data check!!!!!!!!!!!!!!!!!
  Ogeno = geno
  OSNPinfo = SNPinfo
  if (dim(Ogeno)[2] != dim(OSNPinfo)[1]) {
    stop("N of SNPs in geno data and N of SNPs in SNPinfo data Do Not Match!!")
    
  } else if (dim(OSNPinfo)[2] != 2) {
    stop("SNPinfo data Must Contain 2 columns!!")
  }
  Omono = apply(Ogeno, 2, function(x) {
    y<- x[!is.na(x)]
    length(unique(y))!=1
  })
  Ogeno <- Ogeno[,Omono]
  monoSNPs = OSNPinfo[!Omono,]
  OSNPinfo <- OSNPinfo[Omono,]
  maf = apply(Ogeno, 2, function(x) mean(x,na.rm=TRUE)/2)
  maf_ok=ifelse(maf>=0.5,1-maf,maf)
  maf=maf_ok
  mafprun <- which(maf >= MAFcut)
  geno <- Ogeno[,mafprun]
  SNPinfo <- OSNPinfo[mafprun,]
  # print("split whole sequence into subsegments")
  cutpoints.all <- cutsequence.modi(geno, leng, subSegmSize)
  cutpoints <- cutpoints.all[[1]]
  atfcut <- (cutpoints.all[[2]])
  if (!is.null(atfcut)){
    atfcut <- sort(atfcut)
  }
  cutblock <- cbind(c(1, cutpoints + 1), c(cutpoints, dim(geno)[2]))
  cutblock <- cutblock[-(dim(cutblock)[1]), , drop = FALSE]
  LDblocks <- matrix(NA, dim(SNPinfo)[1], 2)
  # partition each segment into LD blocks
  for (i in 1:dim(cutblock)[1]) {
    nowst <- cutblock[i, 1]
    nowed <- cutblock[i, 2]
    subgeno <- geno[, nowst:nowed]
    subSNPinfo <- SNPinfo[nowst:nowed, ]
    # subbinvec <- CLQD(subgeno, subSNPinfo, CLQcut, clstgap, CLQmode = "Density", codechange = FALSE)
    subbinvec <- CLQD(subgeno, subSNPinfo, CLQcut, clstgap, CLQmode = "Density",codechange = FALSE, checkLargest)
    print('CLQ done!')
    bins <- c(1:max(subbinvec[which(!is.na(subbinvec))]))
    clstlist <- sapply(bins, function(x) which(subbinvec == x), simplify = FALSE)
    clstlist <- lapply(clstlist, sort)  ###
    clstlist <- clstlist[order(sapply(clstlist, min))]  ###
    nowLDblocks <- constructLDblock(clstlist, subSNPinfo)
    # print('constructLDblock done!')
    nowLDblocks <- nowLDblocks + (cutblock[i, 1] - 1)
    nowLDblocks <- nowLDblocks[order(nowLDblocks[, 1]), , drop = FALSE]
    preleng1 <- length(which(!is.na(LDblocks[, 1])))
    LDblocks[(preleng1 + 1):(preleng1 + dim(nowLDblocks)[1]), ] <- nowLDblocks
    print(c(i, dim(cutblock)[1]))
    print(Sys.time())
  }
  doneLDblocks <- LDblocks[which(!is.na(LDblocks[, 1])), , drop = FALSE]
  if (length(atfcut) != 0) {
    newLDblocks <- matrix(NA, dim(SNPinfo)[1], 2)
    for(i in 1:(dim(doneLDblocks)[1]-1)){
      # if(i==1080) break;
      print(paste(i, dim(doneLDblocks)[1]))
      #if(i == 1){
      endblock = doneLDblocks[i,]
      #}
      # if(nowblock[1]<end) next;
      nextblock = doneLDblocks[(i+1),]
      gap = c(endblock[2]:nextblock[1])
      if(length(intersect(gap, atfcut))>0){
        ## merge 
        nowatfcut = intersect(gap, atfcut)
        newbigblock = range(c(endblock, nextblock))
        newbigblocksize = diff(newbigblock)+1
        nowst = newbigblock[1]
        nowed = newbigblock[2]
        subgeno <- geno[, nowst:nowed]
        subSNPinfo <- SNPinfo[nowst:nowed, ]
        subbinvec <- CLQD(subgeno, subSNPinfo, CLQcut, clstgap, CLQmode = "Density", codechange = FALSE, checkLargest)
        # print('CLQ done!')
        bins <- c(1:max(subbinvec[which(!is.na(subbinvec))]))
        clstlist <- sapply(bins, function(x) which(subbinvec == x), simplify = FALSE)
        clstlist <- lapply(clstlist, sort)  ###
        clstlist <- clstlist[order(sapply(clstlist, min))]  ###
        nowLDblocks <- constructLDblock(clstlist, SNPinfo[nowst:nowed, ])
        # print('constructLDblock done!')
        nowLDblocks <- nowLDblocks + (newbigblock[1] - 1)
        nowLDblocks <- nowLDblocks[order(nowLDblocks[, 1]), , drop = FALSE]
        preleng1 <- length(which(!is.na(newLDblocks[, 1])))
        nowLDbleng = dim(nowLDblocks)[1]
        #if(nowLDbleng != 1){
        newLDblocks[(preleng1 + 1):(preleng1 + dim(nowLDblocks)[1]), ] <- nowLDblocks
        #}
        #endblock <- nowLDblocks[nowLDbleng, ,drop = FALSE]
        # end <- max(nowLDblocks)
        # if(diff(newbigblock)+1 < subSegmSize)
        print(Sys.time())
      }else{
        addlinei = min(which(is.na(newLDblocks[,1])==TRUE))
        newLDblocks[addlinei,] <-endblock
        #endblock <- nextblock
        if(i == (dim(doneLDblocks)[1]-1)){
          addlinei = min(which(is.na(newLDblocks[,1])==TRUE))
          newLDblocks[addlinei,] <-nextblock
        }
      }
    }
    LDblocks = newLDblocks
  }else{
    LDblocks = doneLDblocks
  }
  LDblocks <- LDblocks[which(!is.na(LDblocks[, 1])), , drop = FALSE]
  LDblocks <- LDblocks[order(LDblocks[, 1]), , drop = FALSE]
  #overlapping LD block merging
  i = 1
  while(i <dim(LDblocks)[1]){
    nowb = LDblocks[i,]
    nextb = LDblocks[(i+1),]
    if(nowb[2]<nextb[1]){
      i = i+1
      next
    }else{
      newb = c(min(c(nowb, nextb)), max(c(nowb, nextb)))
      LDblocks[i,]<-newb
      LDblocks[(i+1),]<-c(NA, NA)
      LDblocks<-LDblocks[!is.na(LDblocks[,1]),]
    }
  }
  start.bp <- SNPinfo[, 2][LDblocks[, 1]]
  end.bp <- SNPinfo[, 2][LDblocks[, 2]]
  start.rsID <- as.character(SNPinfo[, 1][LDblocks[, 1]])
  end.rsID <- as.character(SNPinfo[, 1][LDblocks[, 2]])
  LDblocks <- data.frame(LDblocks, start.rsID, end.rsID, start.bp, end.bp)
  colnames(LDblocks) <- c("start", "end", "start.rsID", "end.rsID", "start.bp", "end.bp")
  if(appendrare==TRUE){
    LDblocks<-appendSGTs(LDblocks, Ogeno, OSNPinfo, CLQcut=CLQcut, clstgap = clstgap, checkLargest)
  }
  allSNPbp = sort(c(monoSNPs[,2], OSNPinfo[,2]))
  LDblocks$start <- match(LDblocks$start.bp, allSNPbp)
  LDblocks$end <- match(LDblocks$end.bp, allSNPbp)
  return(LDblocks)
}
