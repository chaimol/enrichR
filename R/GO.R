#if(!require(ggplot2))install.packages("ggplot2")
#if(!require(tidyverse))install.packages("tidyverse")
#if(!require(ggprism))remotes::install_github("csdaw/ggprism")
#require(ggplot2)
#require(tidyverse)
#require(ggprism) #可以完善ggplot2的图使之达到发表级别
#' GOenrich: for GO enrich analysis
#' GOenrich ,KEGG_enrich,TF_enrich, have the some type input and output
#' @import clusterProfiler
#' @importFrom ggplot2 ggsave
#' @param genelist a gene list like this"
#' MELO3C000003
#' MELO3C000548
#' MELO3C001803
#' @param genetype a string like "WTvsT1" or "Mut12vsWT" for the output image file prefix
#' @param GOfile a filename "GO.info",default:"Go.info".GO.info must be only have 4 column,not have header. separator must is tab.
#' "GO.info" Demo:
#' #GO:0016020	membrane	MELO3C000003	cellular_component
#' #GO:0016021	integral	MELO3C000003	cellular_component
#'
#' @return return GO enrich results
#' @export GOenrich
#'
#' @examples
#' genelist <- read.delim("data/ST1.txt")
#' GOenrich(genelist,"WTvsT1","data/GO.info")
#' GOenrich(genelist,"WTvsT1")
#' go_res <- GOenrich("gene.xls","WTvsT1")
#' if gO_res$Count >1,mean GO enrich is successful.
#' if go_res$Count >1,mean no gene can be enrich.
#'
GOenrich <- function(genelist,genetype,GOfile="GO.info"){
  golist <- read.delim(file=GOfile,sep = "\t",header = FALSE)
  ###删除GO的描述字符中，字符长度超过60的，逗号分割的第一个逗号后的字符
  term <- c("term")
  for (i in 1:length(golist[,2])){
    if(nchar(golist[,2][i])>60){
      term <- append(term,golist[,2][i])
      #print(golist[,2][i])
      golist[,2][i] <- strsplit(golist[,2][i],split= ",")[[1]][1]
    }
  }
  unique(term)

  gene <- genelist
  spo2gene <- data.frame(golist[,1],golist[,3]) #GO,Geneid
  spo2name <- data.frame(golist[,1],golist[,2]) #GO,GO描述

  gene1 <- t(gene)  #获取需要分析的基因列表
  #GO富集分析
  spo_GO <- enricher(gene1,TERM2GENE = spo2gene,TERM2NAME = spo2name)
  res <- length(spo_GO$Count)
  print(res)
  #判断是否有富集结果
  if (res > 0){
    p1 <- dotplot(spo_GO,showCategory=length(spo_GO$Count))+   #showCategory设置出图，显示的行数
      scale_y_discrete(labels=function(x) str_wrap(x, width=45)) #设置term的description的宽度，自动换行
    p1 #GO出图
    ggsave(p1,file=paste(genetype,".GO.tiff",sep=""))
    ggsave(p1,file=paste(genetype,".GO.pdf",sep=""))
    p2 <- cnetplot(spo_GO, circular = TRUE, colorEdge = TRUE)
    ggsave(p2,file=paste(genetype,".GO.circle.tiff",sep=""))
    ggsave(p2,file=paste(genetype,".GO.circle.pdf",sep=""))
    write.csv(spo_GO,file=paste(genetype,"_GO.csv",sep = ""),row.names = FALSE)
    work_path <- getwd()
    message("Output file is in ",work_path,"/",paste0(genetype,"_GO.csv"))
    message("Output file is in ",work_path,"/",paste0(genetype,".GO.circle.pdf"))
    message("Output file is in ",work_path,"/",paste0(genetype,".GO.pdf"))
  }
  return(spo_GO)
}


