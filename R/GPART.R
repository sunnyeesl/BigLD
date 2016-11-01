

#################################################################################################################
# GPART < input >
#' @title Partitioning genodata based on the result obtained by using Big-LD and gene region information.
#' @name GPART
#' @description
#' \code{GPART} partition the given genodata using the result obtained by Big-LD and gene region information.
#' The algorithm partition the whole sequence into sub sequences of which size do not exceed the given threshold.
#'
#' @param geno A data frame or matrix of additive genotype data, each column is additive genotype of each SNP.
#' @param SNPinfo A data frame or matrix of SNPs information.  1st column is rsID and 2nd column is bp position.
#' @param BigLDresult A data frame obtained by \code{Big_LD} function.
#'  If \code{NULL}(default), the \code{GPART} function first excute \code{Big_LD} function to obtain LD blocks estimation result.
#'
#' @param chrN A integer value to specify chromosome number of the given data.
#' @param allgenelist A data frame or matrix of Gene info data.
#' (1st col : Genename, 2nd col : chromosome, 3rd col : start bp,  4th col : end bp)
#' @param minsize A integer value specifying the lower bound of number of SNPs in a partition.
#' @param maxsize A integer value specifying the upper bound of number of SNPs in a partition.
# < output >
#' @return
#'
#' @author Sun Ah Kim <sunny03@snu.ac.kr>
#' @seealso \code{\link{Big_LD}}
#'
#' @examples
#' data(geno)
#' data(SNPinfo)
#' data(allgenelist)
#' GPART(geno, SNPinfo, 22, allgenelist)
#'

#' @importFrom plyr alply

#' @export
GPART <- function(geno, SNPinfo, chrN, allgenelist, BigLDresult = NULL, minsize = 4, maxsize = 50) {
    #######################################################################################################
    # sub-Functions 1. Big_LD 2. LDblockSplit 3. Merge_overlap_Gene 4.LDblock_Gene_Merge 5. split_Big_LD
    # 6. Merge_small_region 7.Naming_Region split LD
    #######################################################################################################
    # blocks with Big-LD result
    #             library(plyr)
    LDblockSplit = function(geno, LDblocks, maxsize) {
        LDblockSizes = apply(LDblocks, 1, diff) + 1
        LargeBlocksN = which(LDblockSizes > maxsize)
        LargeBlocks = LDblocks[LargeBlocksN, , drop = FALSE]
        Newblocks = NULL
        if (length(LargeBlocksN) != 0) {
            while (dim(LargeBlocks)[1] > 0) {
                nowblocks = LargeBlocks[1, ]
                nowgeno = geno[, nowblocks[1]:nowblocks[2]]
                nowr2 = cor(nowgeno)^2
                btwr2 = sapply(1:(dim(nowr2)[1] - 1), function(x) mean(nowr2[1:x, (x + 1):(dim(nowr2)[1] - 1)]))
                if (dim(nowr2)[1] <= ceiling(maxsize * 0.8 * 2)) {
                  subbtwr2 = btwr2[(ceiling(dim(nowr2)[1]/2) - (maxsize * 0.2)):(ceiling(dim(nowr2)[1]/2) + (maxsize * 0.2))]
                  addN = which(subbtwr2 == min(subbtwr2)) + (ceiling(dim(nowr2)[1]/2) - (maxsize * 0.2)) - 1
                  Newblocks = rbind(Newblocks, c(nowblocks[1], (nowblocks[1] + addN)))
                  Newblocks = rbind(Newblocks, c(nowblocks[1] + addN + 1, nowblocks[2]))
                  LargeBlocks <- LargeBlocks[-1, , drop = FALSE]
                } else {
                  subbtwr2 = btwr2[(maxsize * 0.4):maxsize]
                  addN = which(subbtwr2 == min(subbtwr2)) + (maxsize * 0.4) - 1
                  Newblocks = rbind(Newblocks, c(nowblocks[1], (nowblocks[1] + addN - 1)))
                  LargeBlocks[1, 1] <- (nowblocks[1] + addN)
                  if (diff(LargeBlocks[1, ] < maxsize)) {
                    Newblocks = rbind(Newblocks, LargeBlocks[1, ])
                  }
                }
            }  # end while
            FinalLDblocks = rbind(LDblocks[-LargeBlocksN, ], Newblocks)
            FinalLDblocks = FinalLDblocks[order(FinalLDblocks[, 1]), ]
            return(FinalLDblocks)
        } else {
            return(LDblocks)
        }
    }
    # merge overlapped gene regions
    Merge_overlap_Gene = function(Geneblocks) {
        overlaplist = NULL
        for (i in 1:(dim(Geneblocks)[1] - 1)) {
            for (j in (i + 1):dim(Geneblocks)[1]) {
                gene1 = Geneblocks[i, 1]:Geneblocks[i, 2]
                gene2 = Geneblocks[j, 1]:Geneblocks[j, 2]
                if (length(intersect(gene1, gene2)) > 0) {
                  overlaplist = rbind(overlaplist, c(i, j))
                }
            }
        }
        if(is.null(overlaplist)){
          return(Geneblocks)
        }else{
          overlaplist = as.list(data.frame(t(overlaplist)))
          newlist = NULL
          while (length(overlaplist) > 0) {
            now = overlaplist[[1]]
            N = sapply(overlaplist, function(x) length(intersect(x, now)) != 0)
            if (sum(N) == 1) {
              newlist <- c(newlist, list(now))
              overlaplist = overlaplist[-which(N == TRUE)]
            } else {
              newlist <- c(newlist, list(unique(unlist(overlaplist[N]))))
              overlaplist = overlaplist[-which(N == TRUE)]
            }
          }

          overlapN = unlist(newlist)
          RemainGeneblocks = Geneblocks[-overlapN, ]
          AddGeneblocks = NULL
          for (i in 1:length(newlist)) {
            nowoverlap = newlist[[i]]
            overlapgenes = Geneblocks[nowoverlap, ]
            Totalrange = range(as.vector(overlapgenes))
            newGeneblocks = matrix(Totalrange, ncol = 2)
            rownames(newGeneblocks) = paste(rownames(overlapgenes), collapse = "/")
            AddGeneblocks = rbind(AddGeneblocks, newGeneblocks)
          }
          Geneblocks = rbind(RemainGeneblocks, AddGeneblocks)
          Geneblocks = Geneblocks[order(Geneblocks[, 1]), ]
          return(Geneblocks)
        }
    }
    # LDblock construction and split Large regions
    LDblock_Gene_Merge = function(LDblocks.T, Geneblocks) {
        LDblocks.T = data.frame(LDblocks.T, "No", 0)
        Geneblocks = data.frame(Geneblocks, rownames(Geneblocks), 1)
        colnames(LDblocks.T) = c("st", "ed", "gname", "gnameN")
        colnames(Geneblocks) = c("st", "ed", "gname", "gnameN")
        remainedLDPart = NULL
        nowGeneN = 1
        while (dim(LDblocks.T)[1] > 0) {
            nowst = LDblocks.T[1, 1]
            nowed = LDblocks.T[1, 2]
            nowLD = LDblocks.T[1, ]
            # nowGst = as.numeric(nowGene[1]) nowGed = as.numeric(nowGene[2])
            nowGene = Geneblocks[nowGeneN, , drop = FALSE]
            minGeneSNPid = Geneblocks[1, 1]
            maxGeneSNPid = max(Geneblocks[, 2])
            intersectTrue = (length(intersect(c(nowst:nowed), c(nowGene[1, 1]:nowGene[1, 2]))) > 0)
            if (intersectTrue == FALSE) {
                if (nowed < nowGene[1, 1]) {
                  remainedLDPart = rbind(remainedLDPart, nowLD)
                  LDblocks.T = LDblocks.T[-1, , drop = FALSE]
                } else if (nowst > maxGeneSNPid) {
                  remainedLDPart = rbind(remainedLDPart, LDblocks.T)
                  LDblocks.T = LDblocks.T[-(dim(LDblocks.T)[1]), , drop = FALSE]
                } else if (nowst > nowGene[1, 2]) {
                  nowGeneN = nowGeneN + 1
                }
            } else {
                geneset = c(nowGene[1, 1]:nowGene[1, 2])
                LDset = c(nowLD[1, 1]:nowLD[1, 2])
                # print(paste('small!!',length(unique(c(geneset,LDset)))))
                if (nowGene[1, 1] < nowLD[1, 1]) {
                  g1 = as.character(nowGene[1, 3])
                  g2 = as.character(nowLD[1, 3])
                } else {
                  g2 = as.character(nowGene[1, 3])
                  g1 = as.character(nowLD[1, 3])
                }

                if (g1 == "No" & g2 == "No") {
                  gname = "No"
                } else if (g1 == "No" & g2 != "No") {
                  gname <- g2
                } else if (g1 != "No" & g2 == "No") {
                  gname <- g1
                } else if (g1 != "No" & g2 != "No") {
                  g1 = unlist(strsplit(g1, split = "-"))
                  g2 = unlist(strsplit(g2, split = "-"))
                  gname = paste(unique(c(g1, g2)), collapse = "-")
                }
                gnameN = nowLD[1, 4] + nowGene[1, 4]
                st = min(c(nowGene[1, 1], nowGene[1, 2], nowLD[1, 1], nowLD[1, 2]))
                ed = max(c(nowGene[1, 1], nowGene[1, 2], nowLD[1, 1], nowLD[1, 2]))
                if (nowed > nowGene[2]) {
                  levels(LDblocks.T$gname) <- c(levels(LDblocks.T$gname), as.character(gname))
                  levels(LDblocks.T$gnameN) <- c(levels(LDblocks.T$gname), gnameN)
                  LDblocks.T[1, ] <- data.frame(st, ed, as.character(as.factor(gname)), gnameN)
                  Geneblocks <- Geneblocks[-nowGeneN, ]
                  # print(LDblocks.T[1,])
                } else {
                  levels(Geneblocks$gname) <- c(levels(Geneblocks$gname), gname)
                  levels(Geneblocks$gnameN) <- c(levels(Geneblocks$gname), gnameN)
                  Geneblocks[nowGeneN, ] <- data.frame(st, ed, gname, gnameN)
                  LDblocks.T = LDblocks.T[-1, , drop = FALSE]
                  # print(Geneblocks[nowGeneN,])
                }
            }
            if (nowGeneN > dim(Geneblocks)[1]) {
                remainedLDPart = rbind(remainedLDPart, LDblocks.T)
                break
            }
        }
        return(list(remainedLDPart, Geneblocks))
    }

    # split big Gene region
    split_Big_LD = function(GeneLDblocks, LDblocks.T, Geneblocks, maxsize) {
        BigN = which(GeneLDblocks[, 5] > maxsize)
        BigGeneblocks = GeneLDblocks[BigN, ]
        GeneLDblocks = GeneLDblocks[-BigN, ]
        newblocks = NULL
        for (i in 1:dim(BigGeneblocks)[1]) {
            nowBlock = BigGeneblocks[i, ]
            nowSNPs = nowBlock[1, 1]:nowBlock[1, 2]
            intersectLD = apply(LDblocks.T, 1, function(x) length(intersect(c(x[1]:x[2]), nowSNPs)) > 0)
            intersectLD = LDblocks.T[intersectLD, ]
            intersectLD[1, 1] <- min(nowSNPs)
            intersectLD[dim(intersectLD)[1], 2] <- max(nowSNPs)
            NintersectLD = NULL
            while (dim(intersectLD)[1] > 0) {
                st = intersectLD[1, 1]
                edposi = max(which(intersectLD[, 2] < (st + maxsize)))
                ed = intersectLD[edposi, 2]
                NintersectLD = rbind(NintersectLD, c(st, ed))
                intersectLD = intersectLD[-(1:edposi), , drop = FALSE]
            }
            intersectLD = NintersectLD
            Genenames = alply(intersectLD, 1, function(x) {
                GeneR = apply(Geneblocks, 1, function(y) length(intersect(x[1]:x[2], y[1]:y[2])) > 0)
                list(rownames(Geneblocks[GeneR, , drop = FALSE]))
            })
            gname = sapply(Genenames, function(x) paste(x[[1]], collapse = "/"))
            gnameN = sapply(Genenames, length)
            Addblocks = data.frame(intersectLD, gname, gnameN)
            blockL = apply(intersectLD, 1, diff) + 1
            Addblocks = cbind(Addblocks, blockL)
            newblocks = rbind(newblocks, Addblocks)
        }
        GeneLDblocks = rbind(GeneLDblocks, newblocks)
        GeneLDblocks = GeneLDblocks[order(GeneLDblocks[, 1]), ]
        return(GeneLDblocks)

    }
    # small region merging
    Merge_small_region = function(GeneLDblocks, maxsize, minsize) {
        smallbin = NULL
        completeBin = NULL
        while (dim(GeneLDblocks)[1] > 0) {
            # print(dim(GeneLDblocks))
            nowbin = GeneLDblocks[1, ]
            if (is.null(smallbin) & nowbin[1, 5] >= minsize) {
                completeBin = rbind(completeBin, nowbin)
                GeneLDblocks = GeneLDblocks[-1, ]
            } else if (is.null(smallbin) & nowbin[1, 5] < minsize) {
                smallbin <- nowbin
                GeneLDblocks = GeneLDblocks[-1, ]
            } else if (smallbin[1, 5] + nowbin[1, 5] < minsize) {
                st = min(smallbin[1, 1], smallbin[1, 2], nowbin[1, 1], nowbin[1, 2])
                ed = max(smallbin[1, 1], smallbin[1, 2], nowbin[1, 1], nowbin[1, 2])
                g1 = as.character(smallbin[1, 3])
                g2 = as.character(nowbin[1, 3])
                if (g1 == "No" & g2 == "No") {
                  gname = "No"
                } else if (g1 == "No" & g2 != "No") {
                  gname <- g2
                } else if (g1 != "No" & g2 == "No") {
                  gname <- g1
                } else if (g1 != "No" & g2 != "No") {
                  g1 = unlist(strsplit(g1, split = "-"))
                  g2 = unlist(strsplit(g2, split = "-"))
                  gname = paste(unique(c(g1, g2)), collapse = "-")
                }
                gnameN = smallbin[1, 4] + nowbin[1, 4]
                blockL = smallbin[1, 5] + nowbin[1, 5]
                smallbin = data.frame(st, ed, gname, gnameN, blockL)
                GeneLDblocks = GeneLDblocks[-1, ]
            } else if (smallbin[1, 5] + nowbin[1, 5] >= minsize & smallbin[1, 5] + nowbin[1, 5] <= maxsize) {
                st = min(smallbin[1, 1], smallbin[1, 2], nowbin[1, 1], nowbin[1, 2])
                ed = max(smallbin[1, 1], smallbin[1, 2], nowbin[1, 1], nowbin[1, 2])
                g1 = as.character(smallbin[1, 3])
                g2 = as.character(nowbin[1, 3])
                if (g1 == "No" & g2 == "No") {
                  gname = "No"
                } else if (g1 == "No" & g2 != "No") {
                  gname <- g2
                } else if (g1 != "No" & g2 == "No") {
                  gname <- g1
                } else if (g1 != "No" & g2 != "No") {
                  g1 = unlist(strsplit(g1, split = "-"))
                  g2 = unlist(strsplit(g2, split = "-"))
                  gname = paste(unique(c(g1, g2)), collapse = "-")
                }
                gnameN = smallbin[1, 4] + nowbin[1, 4]
                blockL = smallbin[1, 5] + nowbin[1, 5]
                completeBin = rbind(completeBin, data.frame(st, ed, gname, gnameN, blockL))
                smallbin = NULL
                GeneLDblocks = GeneLDblocks[-1, ]
            } else if (smallbin[1, 5] + nowbin[1, 5] > maxsize) {
                lastbin = completeBin[dim(completeBin)[1], ]
                if (lastbin[1, 5] + smallbin[1, 5] <= maxsize) {
                  st = min(smallbin[1, 1], smallbin[1, 2], lastbin[1, 1], lastbin[1, 2])
                  ed = max(smallbin[1, 1], smallbin[1, 2], lastbin[1, 1], lastbin[1, 2])
                  g2 = as.character(smallbin[1, 3])
                  g1 = as.character(lastbin[1, 3])
                  if (g1 == "No" & g2 == "No") {
                    gname = "No"
                  } else if (g1 == "No" & g2 != "No") {
                    gname <- g2
                  } else if (g1 != "No" & g2 == "No") {
                    gname <- g1
                  } else if (g1 != "No" & g2 != "No") {
                    g1 = unlist(strsplit(g1, split = "-"))
                    g2 = unlist(strsplit(g2, split = "-"))
                    gname = paste(unique(c(g1, g2)), collapse = "-")
                  }
                  gnameN = smallbin[1, 4] + lastbin[1, 4]
                  blockL = smallbin[1, 5] + lastbin[1, 5]
                  completeBin = completeBin[-dim(completeBin)[1], ]
                  completeBin = rbind(completeBin, data.frame(st, ed, gname, gnameN, blockL))
                  smallbin = NULL
                } else {
                  print(c("Large-small-Large"))
                  print(rbind(lastbin, smallbin, nowbin))
                  completeBin = rbind(completeBin, smallbin, nowbin)
                  smallbin = NULL
                  GeneLDblocks = GeneLDblocks[-1, ]
                }
            }
        }
        return(completeBin)
    }
    # naming each region
    Naming_Region = function(GeneLDblocks) {
        FinalGeneLDblocks = NULL
        PartN = 1
        Pgene = NULL
        gene = NULL
        Ngene = NULL
        for (i in 1:dim(GeneLDblocks)[1]) {
            if (i%%10 == 0)
                print(i)
            nowBlock = GeneLDblocks[i, ]
            if (is.null(Pgene) & nowBlock$gname == "No") {
                # before First Gene region
                if (is.null(Ngene)) {
                  Ngene = GeneLDblocks$gname[min(which(GeneLDblocks$gname != "No"))]
                  Ngene = unlist(strsplit(as.character(Ngene), "-"))[1]
                  Ngene = unlist(strsplit(as.character(Ngene), "/"))[1]
                }
                Rname = paste("before-", Ngene, sep = "")
                nowFinal = data.frame(nowBlock, Rname)
                FinalGeneLDblocks = rbind(FinalGeneLDblocks, nowFinal)
            } else if (nowBlock$gname != "No") {
                gene = nowBlock$gname
                Rname = gene
                nowFinal = data.frame(nowBlock, Rname)
                FinalGeneLDblocks = rbind(FinalGeneLDblocks, nowFinal)
                Pgene = tail(unlist(strsplit(as.character(gene), "-")), n = 1)
                Pgene = tail(unlist(strsplit(as.character(Pgene), "/")), n = 1)
            } else if (nowBlock$gname == "No") {
                geneNames = GeneLDblocks$gname[(i:dim(GeneLDblocks)[1])]
                if (all(geneNames == "No")) {
                  # after Last Gene region
                  LastPart = GeneLDblocks[i:(dim(GeneLDblocks)[1]), ]
                  Rname = paste("after-", Pgene, sep = "")
                  nowFinal = data.frame(LastPart, Rname)
                  FinalGeneLDblocks = rbind(FinalGeneLDblocks, nowFinal)
                  break
                } else {
                  Ngene = geneNames[min(which(geneNames != "No"))]
                  Ngene = unlist(strsplit(as.character(Ngene), "-"))[1]
                  Ngene = unlist(strsplit(as.character(Ngene), "/"))[1]
                  Rname = paste("inter-", Pgene, "-", Ngene, sep = "")
                  nowFinal = data.frame(nowBlock, Rname)
                  FinalGeneLDblocks = rbind(FinalGeneLDblocks, nowFinal)
                }
            }
        }
        # part numbering
        FinalGeneLDblocks = data.frame(FinalGeneLDblocks, 0)
        FinalGeneLDblocks[1, 7] <- 1
        for (i in 2:dim(FinalGeneLDblocks)[1]) {
            ifelse(FinalGeneLDblocks[i - 1, 6] == FinalGeneLDblocks[i, 6],
                   FinalGeneLDblocks[i, 7] <- FinalGeneLDblocks[(i - 1), 7] + 1, FinalGeneLDblocks[i, 7] <- 1)
            # if(i%%10==0) print(i)
        }


        Finalnames = apply(FinalGeneLDblocks, 1, function(x) {
            # paste(as.character(x[6]), '-part', as.character(x[7]), sep='')
            paste(c(x[6], "-part", as.numeric(x[7])), collapse = "")
        })
        FinalGeneLDblocks = cbind(FinalGeneLDblocks[, 1:5], Finalnames)
        return(FinalGeneLDblocks)
    }

    #######################################################################################################
    # Main part
    if(is.null(BigLDresult)){
      print("Start to execute Big_LD!")
      BigLDblocks = Big_LD(geno, SNPinfo)
      print("Big-LD, done!")
    }else{
      BigLDblocks = BigLDresult
    }
    LDblocks = cbind(as.integer(as.character(BigLDblocks[, 1])), as.integer(as.character(BigLDblocks[, 2])))
    # Split Large LDblocks
    FinalLDblocks = LDblockSplit(geno, LDblocks, maxsize)
    # gene-base block partitioning
    genelist <- allgenelist[which(allgenelist[,2] == chrN), ]
    GeneRegionSNPs1 = NULL
    for (i in 1:dim(genelist)[1]) {
        test = (which(SNPinfo[, 2] >= genelist[i, 3] & SNPinfo[, 2] <= genelist[i, 4]))
        ifelse(length(test) > 0, test <- range(test), test <- c(0, 0))
        GeneRegionSNPs1 = rbind(GeneRegionSNPs1, test)
        # if(i%%10) print(i)
    }

    rownames(GeneRegionSNPs1) <- genelist[, 1]
    SNPexist = apply(GeneRegionSNPs1, 1, function(x) (x[1] != 0 & x[2] != 0))
    SNPexist = unlist(SNPexist)
    Geneblocks = GeneRegionSNPs1[SNPexist, ]
    Geneblocks.M = Merge_overlap_Gene(Geneblocks)
    LDblocks.sgt = setdiff(1:dim(SNPinfo)[1], unlist(apply(FinalLDblocks, 1, function(x) min(x):max(x))))
    LDblocks.sgt = cbind(LDblocks.sgt, LDblocks.sgt)
    LDblocks.T = rbind(FinalLDblocks, LDblocks.sgt)
    LDblocks.T = LDblocks.T[order(LDblocks.T[, 1]), ]
    colnames(LDblocks.T) = c("st", "ed")
    GeneLDblocks = LDblock_Gene_Merge(LDblocks.T, Geneblocks.M)
    GeneLDblocks = rbind(GeneLDblocks[[1]], GeneLDblocks[[2]])
    GeneLDblocks = GeneLDblocks[order(GeneLDblocks[, 1]), ]
    blockL = GeneLDblocks[, 2] - GeneLDblocks[, 1] + 1
    GeneLDblocks = cbind(GeneLDblocks, blockL)
    # split big block
    GeneLDblocks = split_Big_LD(GeneLDblocks, LDblocks.T, Geneblocks, maxsize)
    GeneLDblocks = Merge_small_region(GeneLDblocks, maxsize, minsize)
    GeneLDblocks = Naming_Region(GeneLDblocks)
    st.rsID = SNPinfo[GeneLDblocks[, 1], 1]
    ed.rsID = SNPinfo[GeneLDblocks[, 2], 1]
    st.bp = SNPinfo[GeneLDblocks[, 1], 2]
    ed.bp = SNPinfo[GeneLDblocks[, 2], 2]
    Finalresult = data.frame(GeneLDblocks$st, GeneLDblocks$ed, st.rsID, ed.rsID, st.bp, ed.bp, GeneLDblocks$blockL, GeneLDblocks$Finalnames)
    colnames(Finalresult) = c("st", "ed", "st.rsID", "ed.rsID", "st.bp", "ed.bp", "blocksize", "Name")
    return(Finalresult)
}
