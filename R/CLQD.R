#################################################################################################################
# CLQD <input>
#' @title partitioning into cliques
#' @name CLQD
#'
#' @description \code{CLQD} partitioning the given data into subgroups that contain SNPs which are highly correlated.
#'
#' @param subgeno  A data frame or matrix of additive genotype data, each column is additive genotype of each SNP.
#' @param subSNPinfo  A data frame or matrix of SNPs information. 1st column is rsID and 2nd column is bp position.
#' @param CLQcut A numeric value of threshold for the correlation value |r|, between 0 to 1.
#' @param clstgap An integer value to specifing the threshold of physical distance (bp) between two consecutive SNPs
#' which do not belong to the same clique, i.e., if a physical distance between two consecutive SNPs in a clique
#' greater than \code{clstgap}, then the algorithm split the cliques satisfying each
#' clique do not contain such consecutive SNPs
#' @param CLQmode A character string to specify the way to give priority among detected cliques.
#' if \code{CLQmode = "Density"} then the algorithm gives priority to the clique of largest value of \eqn{(Number of SNPs)/(range of clique)},
#' else if \code{CLQmode = "Maximal"}, then the algorithm gives priority to the largest clique.
#' @param codechange If \code{TRUE}, choose the cliques after code change procedure.
#' @param checkLargest If \code{checkLargest = TRUE}, the algorithm use heuristic procedure to reduce runtime
# <output>
#' @return A vector of cluster numbers of all SNPs (\code{NA} represents singleton cluster)
#' @examples
#'
#' data(geno)
#' data(SNPinfo)
#' CLQD(geno,SNPinfo,CLQcut = 0.5, clstgap= 40000, CLQmode = 'Maximal', codechange = FALSE)
#' CLQD(geno,SNPinfo,CLQcut = 0.5, clstgap= 40000, CLQmode = 'Density', codechange = FALSE)
#'
#' @author Sun-Ah Kim <sunny03@snu.ac.kr>, Yun Joo Yoo <yyoo@snu.ac.kr>
#'
#'

#' @import igraph
#'
# [1] 12 9 13 9 13 6 7 NA 7 7 6 4 8 4 NA 6 13 9 7 7 6 8 12 9 12 13 8 8 12 13
# [31] 7 9 13 8 12 13 7 7 12 NA 10 10 NA 13 10 10 NA 10 18 NA 18 18 NA 18 18 18 18 NA 15 15
# [61] NA 15 16 15 15 15 16 16 16 16 19 16 16 19 NA 17 17 NA NA NA NA 5 5 NA 2 2 1 1 NA 11
# [91] 11 NA NA NA 3 3 NA NA 14 14 subfunctions
# < built-in > 1.CliqueDecision, 2.ChooseMaximal, 3.CodeChangeV, 4.new.split.cliques
#' @export
#' 
CLQD <- function(subgeno, subSNPinfo, CLQcut = 0.5, clstgap = 40000, CLQmode = c("Density", "Maximal"), 
                  codechange = FALSE,  checkLargest = TRUE) {
  # packages
  # library(igraph)
  #######################################################################################################
  # subfunctions : 1.CliqueDecision, 2.ChooseMaximal, 3.CodeChangeV, 4.new.split.cliques 1
  CliqueDecision = function(x, CLQmode) {
    if (length(CLQmode) == 2) {
      CLQmode <<- "Maximal"
      return(length(x))
    } else if (CLQmode == "Maximal") {
      return(length(x))
    } else if (CLQmode == "Density") {
      return(length(x)/(diff(range(x))/1000))
    }
  }
  # 2
  ChooseMaximal = function(vt, cut, OCM) {
    
    codeW <- CodeChangeV(vt, OCM)[[1]]  #use CodeChangeV function
    codeW[codeW < cut] <- 0
    subg <- graph.adjacency(codeW, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NULL)
    lgstcliq <- largest.cliques(subg)
    if (length(lgstcliq) == 1) {
      FC <- unlist(largest.cliques(subg))
    } else {
      sumvt <- sapply(lgstcliq, function(x) {
        sum(codeW[x, x])
      })
      cliqno <- which(sumvt == max(sumvt))[1]
      FC <- lgstcliq[[cliqno]]
    }
    return(vt[FC])
  }
  # 3
  CodeChangeV = function(vt, OCM) {
    rin <- OCM[vt, vt]
    rin <- as.matrix(rin)
    nr <- dim(rin)[1]
    code = rep(1, nr)
    change = 1
    iter = 0
    while (change == 1) {
      change = 0
      if (iter > 2^nr) {
        print("maximum iteration")
        break
      }
      # count number of negative r for each SNP & find the SNP with max neg r
      cneg = apply(rin, 1, function(x) {
        sum(x < 0)
      })
      maxneg = which.max(cneg)
      
      if (cneg[maxneg] > (nr - 1)/2) {
        # change signs of r for that row and column
        rin[maxneg, ] = -1 * rin[maxneg, ]
        rin[, maxneg] = -1 * rin[, maxneg]
        code[maxneg] = 1 - code[maxneg]
        change = 1
      }
      iter = iter + 1
      
    }
    return(list(rin, code))
  }
  # 4
  new.split.cliques <- function(cliques.bp, gapdist) {
    nowlist <- lapply(cliques.bp, sort)
    fixlist <- NULL
    repeat {
      need.split = which(sapply(nowlist, function(x) max(diff(x)) > gapdist) == TRUE)
      need.fix <- which(sapply(nowlist, function(x) max(diff(x)) > gapdist) == FALSE)
      addlist <- nowlist[need.fix]
      fixlist <- c(fixlist, addlist)
      if (length(need.split) == 0) {
        break
      }
      nowlist <- nowlist[need.split]
      nowlength <- length(nowlist)
      newlist <- as.list(rep(NA, nowlength))
      for (i in 1:nowlength) {
        gap = diff(nowlist[[i]])
        frontpart <- nowlist[[i]][1:min(which(gap > gapdist))]
        restpart <- nowlist[[i]][-(1:min(which(gap > gapdist)))]
        nowlist[[i]] <- frontpart
        newlist[[i]] <- restpart
      }
      addlist <- nowlist[sapply(nowlist, function(x) length(x) > 1)]
      fixlist <- c(fixlist, addlist)
      nowlist <- newlist[sapply(newlist, function(x) length(x) > 1)]
    }
    return(fixlist)
  }
  ########################################################################################################
  if (length(CLQmode) == 2) {
    print(" You do not choose CLQ mode! Defalt mode is 'Density'.")
    CLQmode <- "Density"
  }
  # Main Function
  SNPbps = subSNPinfo[, 2]
  OCM <- cor(subgeno)
  diag(OCM) <- 0
  OCM[abs(OCM) < CLQcut] <- 0
  r2Mat <- OCM^2
  r2Mat[r2Mat < CLQcut^2] <- 0
  r2Mat[r2Mat >= CLQcut^2] <- 1
  # Nr2Mat <- r2Mat
  # Nr2Mat[abs(Nr2Mat)<CLQcut] <- 0
  # Nr2Mat[Nr2Mat>0]<-1
  binvector = rep(NA, dim(r2Mat)[2])
  binnum = 1
  re.SNPbps <- SNPbps
  if(all(OCM==0)) return(1:length(binvector))
  ##
  # take Toooo Big block First!
  # g <- graph_from_adjacency_matrix(r2Mat, mode = "undirected", weighted = TRUE, diag = NULL, add.colnames = NA)
  subregionLeng = dim(OCM)[1]
  if(subregionLeng<500) checkLargest = FALSE
  while(checkLargest == TRUE){
    if(checkLargest == TRUE){
      g <- graph_from_adjacency_matrix(r2Mat, mode = "undirected", weighted = TRUE, diag = NULL, add.colnames = NA)
      compo = components(g)
      componum = which(compo$csize==max(compo$csize))
      compov = which(compo$membership==componum)
      compadjM = OCM[compov, compov]
      cg = graph_from_adjacency_matrix(compadjM, mode = "undirected", weighted = TRUE, diag = NULL, add.colnames = NA)
      if((median(coreness(cg))>80 & max(coreness(cg))>100)| (quantile(coreness(cg), 0.75)>100 & max(coreness(cg))>100)){
        degrees = apply(r2Mat, 1, sum)
        maxdegv = which(degrees >=(quantile(degrees, 0.7)))
        # if(length(maxdegv)>=1){
        maxdegvs = maxdegv
        edgeDens = NULL
        for(maxdegv in maxdegvs){
          Bignbds = which(r2Mat[maxdegv,, drop = FALSE]>0, arr.ind = TRUE)
          Bignbds.c = unique(Bignbds[,2])
          newr2Mat = r2Mat[Bignbds.c,Bignbds.c]
          EdgeDen = sum(newr2Mat)/((dim(newr2Mat)[1])*(dim(newr2Mat)[1]-1))
          edgeDens = c(edgeDens, EdgeDen)
        }
        maxdegvs = maxdegvs[order(edgeDens, decreasing = TRUE)]
        edgeDens = edgeDens[order(edgeDens, decreasing = TRUE)]
        degv = maxdegvs[1]
        edgeD = edgeDens[1]
        Bignbds = which(r2Mat[degv,, drop = FALSE]>0, arr.ind = TRUE)
        Bignbds.c = unique(Bignbds[,2])
        # maxiC = maximal.cliques(g, min = dim(OCM)[1]*0.9)
        # largestOneRange = range(Bignbds.c)
        # largestSNPn = diff(largestOneRange)
        # largestCsize = length(Bignbds.c)
        nowSNPsbp = re.SNPbps[Bignbds.c]
        nowSNPsbploca = match(nowSNPsbp, SNPbps)
        binvector[nowSNPsbploca] <- binnum
        binnum = binnum + 1
        r2Mat <- r2Mat[-Bignbds.c, -Bignbds.c]
        OCM <- OCM[-Bignbds.c, -Bignbds.c]
        re.SNPbps <- re.SNPbps[-Bignbds.c]
        # print("case2")
        checkLargest = TRUE
       
      }else{
        checkLargest = FALSE
      }
    }
  }
  print("end pre-steps")
  g <- graph_from_adjacency_matrix(r2Mat, mode = "undirected", weighted = TRUE, diag = NULL, add.colnames = NA)
  max.cliques <- max_cliques(g, min = 2)
  bp.cliques <- lapply(max.cliques, function(x) re.SNPbps[x])
  split.bp.cliques <- new.split.cliques(bp.cliques, clstgap)
  repeat {
    if (all(is.na(binvector) == FALSE)) {
      break
    }
    if(length(split.bp.cliques)==0) break
    density.v <- sapply(split.bp.cliques, function(x) CliqueDecision(x, CLQmode), simplify = TRUE)
    max.d <- which(density.v == max(density.v))
    max.cluster <- split.bp.cliques[max.d]
    if (length(max.cluster) > 1) {
      # if there are two bins of same density, then we choose the bigger one.
      max.cluster <- max.cluster[order(sapply(max.cluster, length), decreasing = TRUE)]
    }
    max.cluster <- max.cluster[[1]]
    max.cluster.od <- match(max.cluster, re.SNPbps)
    if (codechange == TRUE) {
      max.cluster.od <- ChooseMaximal(max.cluster.od, CLQcut, OCM)
      max.cluster <- re.SNPbps[max.cluster.od]
    }
    ## excluding all SNPs in max.cluster from re.SNPbps
    split.bp.cliques <- lapply(split.bp.cliques, function(x) setdiff(x, max.cluster))
    split.bp.cliques <- split.bp.cliques[which(sapply(split.bp.cliques, length) > 1)]
    binvector[match(max.cluster, SNPbps)] <- binnum
    binnum = binnum + 1
    # r2Mat <- r2Mat[-max.cluster.od, -max.cluster.od]
    # OCM <- OCM[-max.cluster.od, -max.cluster.od]
    # re.SNPbps <- setdiff(re.SNPbps, max.cluster)
    if (length(re.SNPbps) < 2) {
      break
    } 
    if(length(split.bp.cliques)==0) break
    # print(sum(is.na(binvector)))
  }  ##end repeat
  
  if (all(is.na(binvector) == TRUE)) {
    binvector <- c(1:length(binvector))
  }
  return(binvector)
}