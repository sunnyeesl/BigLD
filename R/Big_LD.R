
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
#' @param clstgap  An integer value to specifing the threshold of physical distance (bp) between two consecutive SNPs
#' which do not belong to the same clique, i.e., if a physical distance between two consecutive SNPs in a clique
#' greater than \code{clstgap}, then the algorithm split the cliques satisfying each
#' clique do not contain such consecutive SNPs
#' @param leng  An integer value to specify the number of SNPs in a preceding and a following region
#' of each sub-region boundary, every SNP in a preceding and every SNP in a following region need to be in weak LD.
#'
#' @param subSegmSize  An integer value to specify the upper bound of the number of SNPs in a one-take sub-region.
#'
#'
#'
# <output>
#' @return  A data frame of block estimation result.
#' Each row of data frame shows the starting SNP and end SNP of each estimated LD block
#'
#' @author Sun-Ah Kim <sunny03@snu.ac.kr>
#' @seealso \code{\link{CLQD}}, \code{\link{LDblockHeatmap}}
#'
#' @examples
#'
#' data(geno)
#' data(SNPinfo)
#' Big_LD(geno, SNPinfo)
#' Big_LD(geno, SNPinfo, CLQcut = 0.5, clstgap = 40000, leng = 200, subSegmSize = 1500)
#' @export
#'
# sub-Functions 1. CLQD < built - in > 2. cutsequence.modi, 3.intervalCliqueList, 4. find.maximum.indept, 5. constructLDblock,
Big_LD = function(geno, SNPinfo, CLQcut = 0.5, clstgap = 40000, leng = 200, subSegmSize = 1500) {
    # packages
    # library(igraph)
    #######################################################################################################
    # sub-Functions 1. cutsequence.modi, 2.intervalCliqueList, 3. find.maximum.indept, 4. constructLDblock, 5. CLQ
    cutsequence.modi = function(geno, leng, CLQcut, subSegmSize) {
        ## region length<=3000
        if (dim(geno)[2] <= subSegmSize) {
            return(list(dim(geno)[2], NULL))
        } else {
            cutpoints <- NULL
            i = leng
            while (i <= (dim(geno)[2] - leng)) {
                # tick size = leng * (1/10)
                tick <- as.integer(leng * 1/10)
                nowcm <- cor(geno[, (i - tick + 1):(i + tick)])
                nowr2 <- nowcm^2
                nowr2[which(nowr2 < CLQcut)] = 0
                diag(nowr2) <- 0
                if (length(which(nowr2[1:tick, (tick + 1):(2 * tick)] > 0)) > 0) {
                  # print(i)
                  i = i + 1
                  next
                }
                # tick size = leng
                tick <- leng
                nowcm <- cor(geno[, (i - tick + 1):(i + tick)])
                nowr2 <- nowcm^2
                nowr2[which(nowr2 < CLQcut)] = 0
                diag(nowr2) <- 0
                if (length(which(nowr2[1:tick, (tick + 1):(2 * tick)] > 0)) > 0) {
                  i = i + 1
                  next
                } else {
                  print(c(i))
                  cutpoints <- c(cutpoints, i)
                  i <- i + (leng/2)
                }
            }  ##end while
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
                    nowCM <- cor(geno[, (x - tick + 1):(x + tick)])
                    nowr2 <- nowCM^2
                    diag(nowr2) <- 0
                    length(which(nowr2[(1:tick), (tick + 1):(2 * tick)] > CLQcut))
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
            }  ##end while
            print("cutting sequence, done")
            Nregion = length(cutpoints) + length(atfcut)
            print(paste("N of sub-region:", length(cutpoints), "+", length(atfcut)))
            return(list(cutpoints, atfcut))
        }
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
        interval.graph <- graph.adjacency(adjacencyM, mode = "undirected", weighted = TRUE, diag = NULL, add.colnames = NULL)
        interval.cliques <- cliques(interval.graph)
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
                if (dim(max.info)[1] != 1)
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
            if (length(clstlist) == 0)
                break
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
    #######################################################################################################
    # Main part input data check!!!!!!!!!!!!!!!!!
    if (dim(geno)[2] != dim(SNPinfo)[1]) {
        stop("N of SNPs in geno data and N of SNPs in SNPinfo data Do Not Match!!")

    } else if (dim(SNPinfo)[2] != 2) {
        stop("SNPinfo data Must Contain 2 columns!!")
    }
    # split all sequence into subsegments
    cutpoints.all <- cutsequence.modi(geno, leng, CLQcut, subSegmSize)
    cutpoints <- cutpoints.all[[1]]
    atfcut <- (cutpoints.all[[2]])
    if (!is.null(atfcut))
        atfcut <- sort(atfcut)
    cutblock <- cbind(c(1, cutpoints + 1), c(cutpoints, dim(geno)[2]))
    cutblock <- cutblock[-(dim(cutblock)[1]), , drop = FALSE]
    LDblocks <- matrix(NA, dim(SNPinfo)[1], 2)
    # partition each segment into LD blocks
    for (i in 1:dim(cutblock)[1]) {
        nowst <- cutblock[i, 1]
        nowed <- cutblock[i, 2]
        subgeno <- geno[, nowst:nowed]
        subSNPinfo <- SNPinfo[nowst:nowed, ]
        subbinvec <- CLQD(subgeno, subSNPinfo, CLQcut, clstgap, CLQmode = "Density", codechange = FALSE)
        # print('CLQ done!')
        bins <- c(1:max(subbinvec[which(!is.na(subbinvec))]))
        clstlist <- sapply(bins, function(x) which(subbinvec == x), simplify = FALSE)
        clstlist <- lapply(clstlist, sort)  ###
        clstlist <- clstlist[order(sapply(clstlist, min))]  ###
        nowLDblocks <- constructLDblock(clstlist, SNPinfo[nowst:nowed, ])
        # print('constructLDblock done!')
        nowLDblocks <- nowLDblocks + (cutblock[i, 1] - 1)
        nowLDblocks <- nowLDblocks[order(nowLDblocks[, 1]), , drop = FALSE]
        preleng1 <- length(which(!is.na(LDblocks[, 1])))
        LDblocks[(preleng1 + 1):(preleng1 + dim(nowLDblocks)[1]), ] <- nowLDblocks
        print(c(i, dim(cutblock)[1]))
    }
    doneLDblocks <- LDblocks[which(!is.na(LDblocks[, 1])), , drop = FALSE]
    if (length(atfcut) != 0) {
        cutblock <- sapply(atfcut, function(x) {
            st <- max(which(doneLDblocks[, 1] < x))
            stp <- doneLDblocks[, 1][st]
            ed <- min(which(doneLDblocks[, 2] > x))
            edp <- doneLDblocks[, 2][ed]
            c(stp, edp)
        })
        cutblock <- t(cutblock)
        removeSNPs <- unlist(apply(cutblock, 1, function(x) min(x):max(x)))
        removeBlock <- which(sapply(LDblocks, function(x) length(intersect(x, removeSNPs)) > 0) == TRUE)
        LDblocks <- LDblocks[-removeBlock, ]
        for (k in 1:dim(cutblock)[1]) {
            nowst <- cutblock[k, 1]
            nowed <- cutblock[k, 2]
            subgeno <- geno[, nowst:nowed]
            subSNPinfo <- SNPinfo[nowst:nowed, ]
            subbinvec <- CLQD(subgeno, subSNPinfo, CLQcut, clstgap, CLQmode = "Density", codechange = FALSE)
            # print('CLQ done!')
            bins <- c(1:max(subbinvec[which(!is.na(subbinvec))]))
            clstlist <- sapply(bins, function(x) which(subbinvec == x), simplify = FALSE)
            clstlist <- lapply(clstlist, sort)  ###
            clstlist <- clstlist[order(sapply(clstlist, min))]  ###
            nowLDblocks <- constructLDblock(clstlist, SNPinfo[nowst:nowed, ])
            # print('constructLDblock done!')
            nowLDblocks <- nowLDblocks + (cutblock[k, 1] - 1)
            nowLDblocks <- nowLDblocks[order(nowLDblocks[, 1]), , drop = FALSE]
            preleng1 <- length(which(!is.na(LDblocks[, 1])))
            LDblocks[(preleng1 + 1):(preleng1 + dim(nowLDblocks)[1]), ] <- nowLDblocks
            print(c(k, dim(cutblock)[1]))
        }
    }
    LDblocks <- LDblocks[which(!is.na(LDblocks[, 1])), , drop = FALSE]
    LDblocks <- LDblocks[order(LDblocks[, 1]), , drop = FALSE]
    start.bp <- SNPinfo[, 2][LDblocks[, 1]]
    end.bp <- SNPinfo[, 2][LDblocks[, 2]]
    start.rsID <- as.character(SNPinfo[, 1][LDblocks[, 1]])
    end.rsID <- as.character(SNPinfo[, 1][LDblocks[, 2]])
    LDblocks <- data.frame(LDblocks, start.rsID, end.rsID, start.bp, end.bp)
    colnames(LDblocks) <- c("start", "end", "start.rsID", "end.rsID", "start.bp", "end.bp")
    return(LDblocks)
}
