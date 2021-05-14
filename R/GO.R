if(!require(ggplot2))install.packages("ggplot2")
if(!require(tidyverse))install.packages("tidyverse")
if(!require(ggprism))remotes::install_github("csdaw/ggprism")
library(ggplot2)
library(tidyverse)
library(ggprism) #可以完善ggplot2的图使之达到发表级别
#' Return the GO results
#'
#' @param genelist a list like this"
#' MELO3C000003
#' MELO3C000548
#' MELO3C001803
#' @param genetype a string like "WTvsT1" or "Mut12vsWT" for the output file prefix
#' @param GOfile a filename "GO.info",default:"Go.info".GO.info must be only have 4 column,not have header. separator must is tab.
#' "GO.info" like this:
#' #GO:0016020	membrane	MELO3C000003	cellular_component
#' #GO:0016021	integral	MELO3C000003	cellular_component
#'
#' @return GO enrich results
#' @export the XXX_GO.csv , XXX.GO.tiff,XXX.GO.pdf
#'
#' @examples
#' GO_enrich("gene.xls","WTvsT1","GO.info")
#' GO_enrich("gene.xls","WTvsT1")
#' go_res <- GO_enrich("gene.xls","WTvsT1")
#' if gO_res$Count >1,mean GO enrich is successful.
#' if go_res$Count >1,mean no gene can be enrich.
#'
GO_enrich <- function(genelist,genetype,GOfile="GO.info"){
  golist <- read.delim(file=GOfile,sep = "\t",header = FALSE)
  gene <- genelist

  spo2gene <- data.frame(golist$V1,golist$V3) #GO,Geneid
  spo2name <- data.frame(golist$V1,golist$V2) #GO,GO描述
  gene1 <- t(gene)  #获取需要分析的基因列表
  #GO富集分析
  spo_GO <- enricher(gene1,TERM2GENE = spo2gene,TERM2NAME = spo2name)
  res <- length(spo_GO$Count)
  print(res)
  #判断是否有富集结果
  if (res > 1){
    p1 <- dotplot(spo_GO,showCategory=length(spo_GO$Count))+   #showCategory设置出图，显示的行数
      scale_y_discrete(labels=function(x) str_wrap(x, width=45)) #设置term的description的宽度，自动换行
    p1 #GO出图
    ggsave(p1,file=paste(genetype,".GO.tiff",sep=""))
    ggsave(p1,file=paste(genetype,".GO.pdf",sep=""))
    write.csv(spo_GO,file=paste(genetype,"_GO.csv",sep = ""))
  }
  return(spo_GO)
}




