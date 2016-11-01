##################################################################################################################
#                LDheatmap
##################################################################################################################
# < input >
#' @title Draw LDblock results on LD heatmap.
#' @name LDblockHeatmap
#' @description
#' \code{LDblockHeatmap} shows the LDblock regions obtained by Big-LD algorithm.
#'
#' @param geno A data frame or matrix of additive genotype data, each column is additive genotype of each SNP.
#' @param SNPinfo A data frame or matrix of SNPs information.  1st column is rsID and 2nd column is bp position.
#' @param chrN A integer value to specify chromosome number of the given data.
#' @param showSNPs A data frame which is part of \code{SNPinfo} that you want to show in the result LDblock heatmap.
#' The default is \code{NULL}
#' @param LDblockResult  A data frame obtained by \code{Big_LD} function.
#'  If \code{NULL}(default), the \code{GPART} function first excute \code{Big_LD} function to obtain LD blocks estimation result.
#'
#' @param tick A character string to specify how to show first SNPs and last SNPs of LD blocks,
#' in \code{"bp"} or in \code{"rsID"}.
#' @param st.bp A integer value to specify starting bp position of the region to draw.
#' @param ed.bp A integer value to specify end bp position of the region to draw.
#' @param showLDsize A integer value to specify the size (number of SNPs) of LDblocks to show ticks.
#' if \code{showLDsize = \eqn{k}}, then blocks of size equal or greater than k is shown with boundaries and
#' the "rsID" or "bp" of the first and last SNPs, however,
#' the blocks whose size is less than \eqn{k} is shown with only boundaries.
#' @param savefile logical. If \code{TRUE}, save tif file into work directory and
#' the file is named after the chromosome and the physical range to draw. The default is \code{FALSE}
# < output >
#' @return A grid graphical object of LD block heatmap.
#' The LD block heatmap will be presented on the screen after execution of the function.
#' #'
#' @author Sun-Ah Kim <sunny03@snu.ac.kr>
#' @seealso \code{\link{Big_LD}}
#'
#' @examples
#'
#' data(geno)
#' data(SNPinfo)
#' LDblockHeatmap(geno, SNPinfo,chrN = 22, showSNPs = NULL)
#' LDblockHeatmap(geno, SNPinfo, 22, showSNPs = SNPinfo[c(100, 200), ], showLDsize = 10, savefile = TRUE)
#' @export
#' @import grid
# TotalMap : object of heatmap.
# If savefile TRUE, "LDblock_heatmap_chr[chrN]-[start bp]-[end bp].tif" file will be made.
# The function will draw a plot on your screen.
##################################################################################################################
# sub-Functions
# 1. Big-LD    2. CLQD
#  < built - in >
# 2. makeRect (built-in)    3. LDheatmap.Legend.add (built-in)
##################################################################################################################


LDblockHeatmap = function(geno, SNPinfo, chrN, showSNPs = NULL, LDblockResult=NULL, tick = c("bp", "rsID"), st.bp=0 , ed.bp = Inf,
                          showLDsize = 3, savefile = FALSE){
  # packagese
  # library(grid)
  ########################################################################################################
  # sub-Functions
  # 1. Big-LD    2. makeRect (built-in)    3. LDheatmap.Legend.add (built-in)
  ########################################################################################################
  # make rectangles
  makeRect<-function(nrow,ncol,cols,name,byrow=TRUE){
    xx<-(1:ncol)/ncol
    yy<-(1:nrow)/nrow
    right<-rep(xx,nrow)
    top<-rep(yy,each=ncol)
    rectGrob(x=right,y=top,width=1/ncol,height=1/nrow,just=c("right","top"),gp=gpar(col=NA,fill=cols),name=name)
  }
  # make Legend
  LDheatmap.Legend.add <- function(color, vp=VPheatmap){
    ImageRect <- makeRect(2, length(color), col = c(rep(NA, length(color)), color[length(color):1]), "colorKey")
    keyVP <- viewport(x = 0.5, y = 0.2, height = 0.4, width =0.5, just = c("centre", "bottom"), name = "keyVP")
    ttt <- expression(paste("r"^{2}, " Color Key"))
    title <- textGrob(ttt, x = 0.5, y = 1.25, name = "title", gp = gpar(cex = 0.6))
    labels <- textGrob(paste(0.2 * 0:5), x = 0.2 * 0:5, y = 0.25, gp = gpar(cex = 0.6), name = "labels")
    ticks <- segmentsGrob(x0 = c(0:5) * 0.2, y0 = rep(0.4,6), x1 = c(0:5) * 0.2, y1 = rep(0.5, 6), name = "ticks")
    box <- linesGrob(x = c(0, 0, 1, 1, 0), y = c(0.5, 1, 1, 0.5, 0.5), name = "box")
    key <- gTree(children = gList(ImageRect, title, labels, ticks, box), name = "Key", vp = keyVP)
    return(key)
  }
  ########################################################################################################
  if(length(tick) >1 ){
    tick = "rsID"
  }
  subSNPinfo = SNPinfo[which(SNPinfo[,2]>=st.bp & SNPinfo[,2]<=ed.bp),]
  if(dim(subSNPinfo)[1]>1000){
    print("There are Too many SNPs! We will draw only first 1000 SNPs")
    subSNPinfo<- subSNPinfo[1:1000,]
  }else if(dim(subSNPinfo)[1]<10){
    stop("Too short Region!")
  }
  if(!all(colnames(geno)==SNPinfo[,1])){
    stop("column names of geno data do not agree with 1st column of SNPinfo ")
  }

  chosencol = as.vector(match(as.character(subSNPinfo[,1]), colnames(geno)))
  subgeno = geno[,chosencol]
  if(is.null(LDblockResult)){
    subLDblockRes = Big_LD(subgeno, subSNPinfo)
    print("Big_LD, done!")
  }else{
    subLDblockRes = LDblockResult[which(LDblockResult$start.bp >= min(subSNPinfo[,2]) & LDblockResult$end.bp <= max(subSNPinfo[,2])),]
  }

  s = subLDblockRes$start
  e = subLDblockRes$end
  s.rsID = as.character(subLDblockRes$start.rsID)
  e.rsID = as.character(subLDblockRes$end.rsID)
  s.bp = subLDblockRes$start.bp
  e.bp = subLDblockRes$end.bp

  if(s.bp[1]<subSNPinfo[1,2]) {
    s.bp[1]<-subSNPinfo[1,2]
    s.rsID[1]<-as.character(subSNPinfo[1,1])
  }
  if(max(e.bp)>max(subSNPinfo[,2])){
    e.bp[length(e.bp)]<-max(subSNPinfo[,2])
    e.rsID[length(e.bp)]<-as.character(subSNPinfo[dim(subSNPinfo)[1],1])
  }

  subLDblockRes =  data.frame(s, e, s.rsID, e.rsID, s.bp, e.bp)
  colnames(subLDblockRes)=c("start", "end", "start.rsID", "end.rsID", "start.bp", "end.bp")

  BlockstP = sapply(subLDblockRes$start.bp, function(x) which(min(abs(x-subSNPinfo[,2])) ==abs(x-subSNPinfo[,2]))[1])
  BlockedP = sapply(subLDblockRes$end.bp, function(x) which(min(abs(x-subSNPinfo[,2])) ==abs(x-subSNPinfo[,2]))[1])
  if(!is.null(showSNPs)){
    showSNPsP = sapply(showSNPs[,2], function(x) which(abs(subSNPinfo[,2] - x) == min(abs(subSNPinfo[,2] - x)))[1])
  } else {
    showSNPsP = NULL
  }


  if(tick == "rsID"){
    tickname.st = as.character(subLDblockRes$start.rsID)
    tickname.ed = as.character(subLDblockRes$end.rsID)
  } else if(tick =="bp"){
    tickname.st = (subLDblockRes$start.bp)
    tickname.ed = (subLDblockRes$end.bp)
  }
  ticksizeNsatisfy = ((BlockedP-BlockstP+1)<showLDsize)
  BlockstP.tick = BlockstP[!ticksizeNsatisfy]
  BlockedP.tick = BlockedP[!ticksizeNsatisfy]
  tickname.st = tickname.st[!ticksizeNsatisfy]
  tickname.ed = tickname.ed[!ticksizeNsatisfy]
  VPheatmap<-viewport(x=0.05,y=0.1,width=unit(0.85,"snpc"),height=unit(0.85,"snpc"),just=c("left","bottom"), name="VPheatmap")
  VPtick<-viewport(x=0.055,y=0.095,width=unit(0.85,"snpc"),height=unit(0.85,"snpc"),just=c("left","bottom"), name="VPtick")
  VPLegend<-viewport(x=0.65,y=0,width=unit(0.35,"snpc"),height=unit(0.15,"snpc"),just=c("left","bottom"), name="VPLegend")
  ###make rectangle
  color=heat.colors(50)
  M = cor(subgeno)^2
  M[lower.tri(M,diag=FALSE)]<-NA
  mybreak<-0:length(color)/length(color)
  colcut<-as.character(cut(1-M,mybreak,labels=as.character(color), include.lowest=TRUE, include.highest=FALSE))
  ldheatmap<-makeRect(dim(M)[1],dim(M)[2],colcut,"ldheatmap")
  WhiteTri = polygonGrob(x=c(0, 1, 1), y=c(0, 0, 1),
                         id=NULL, id.lengths=NULL, name=NULL, gp=gpar(col = "white", fill = "white"))
  x1s = BlockstP*(1/dim(M)[1])
  x2s = BlockstP*(1/dim(M)[1])
  x3s = BlockedP*(1/dim(M)[1])
  y1s = BlockstP*(1/dim(M)[1])
  y2s = BlockedP*(1/dim(M)[1])
  y3s = BlockedP*(1/dim(M)[1])
  trascolor =   adjustcolor( "red", alpha.f = 0)
  BlockBoundaries = polygonGrob(x = c(x1s, x2s, x3s), y = c(y1s, y2s, y3s), id = rep(1:length(x1s),3), gp = gpar(fill = trascolor))
  blockname1 = tickname.st
  blockname2 = tickname.ed
  tickposi1 =  BlockstP.tick*(1/dim(M)[1])+0.005
  tickposi2 = BlockedP.tick*(1/dim(M)[1])-0.005
  tickname = c(tickname.st, tickname.ed)
  tickposi = c(tickposi1, tickposi2)
  print(length(tickname))
  LDticks = textGrob(tickname, x = tickposi, y = tickposi, just = c("left", "centre"),rot = -45, gp = gpar(cex = 0.5), vp = VPtick)
  LDblockHeatmap = gTree(children=gList(ldheatmap, WhiteTri, BlockBoundaries),name="LDHeatmap", vp = VPheatmap)
  if(is.null(showSNPsP)){
    showSNPposi = NULL
    LDblockMap = gTree(children=gList(LDblockHeatmap, LDticks),name="LDHeatmap")
  }else{
    showSNPposi = showSNPsP*(1/dim(M)[1])
    SNPticks = textGrob(as.character(showSNPs[,1]), x = showSNPposi+0.075, y = showSNPposi-0.075,
                        just = c("left", "centre"), gp = gpar(cex = 0.5, col = "red"), vp = VPtick)
    SNPseg = segmentsGrob(x0 =showSNPposi , y0=showSNPposi, x1= showSNPposi+0.07, y1=showSNPposi-0.07, gp = gpar(col = "red", alpha = 0.5),vp = VPtick )
    LDblockMap = gTree(children=gList(LDblockHeatmap,  SNPseg, SNPticks ,LDticks),name="LDHeatmap")
  }

  Legend = gTree(children = gList(LDheatmap.Legend.add(color=heat.colors(50))), vp = VPLegend)

  regionname = paste("chr",chrN,": ",format(min(subSNPinfo[,2]), big.mark=",",scientific=FALSE) ,
                     "bp ~ ", format(max(subSNPinfo[,2]), big.mark=",",scientific=FALSE),"bp", sep="" )

  Figregion = textGrob(regionname, x = 0.05, y = 0.98, just = c("left", "centre"),gp = gpar(cex = 1))
  TotalMap = gTree(children = gList(LDblockMap, Figregion, Legend))
  # grid.draw(LDblockMap)

  if(savefile ==FALSE){
    grid.draw(TotalMap)
    return(TotalMap)
  } else {
    tiff(paste("LDblock_heatmap_chr",chrN,"-", min(subSNPinfo[,2]), "-", max(subSNPinfo[,2]), ".tif", sep=""),
         res = 300, height=210, width=210, units = 'mm')
    grid.draw(TotalMap)
    dev.off()
    grid.draw(TotalMap)
    return(TotalMap)
  }
}
