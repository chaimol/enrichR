### a function for barplot GO,KEGG results.
rm(list = ls())
sedwd("e:/pakwork/cancer/enrichR/R/")
if(!require(ggplot2))install.packages("ggplot2")
if(!require(tidyverse))install.packages("tidyverse")
#if(!require(ggprism))remotes::install_github("csdaw/ggprism")
library(ggplot2)
library(tidyverse)
#library(ggprism) #可以完善ggplot2的图使之达到发表级别
# if (!require(venn)) install.packages("venn")
# library(venn)
if (!require(devtools)) install.packages("devtools")
if(!require(ggvenn))devtools::install_github("yanlinlin82/ggvenn")
library(ggvenn)
if(!require(UpSetR))install.packages("UpSetR")
library("UpSetR")
require(ggplot2); require(plyr); require(gridExtra); require(grid);

#######输入基因的列，获取0/1或FALSE/TRUE矩阵的数据框########
#' Title get_matrix:return a 0/1 matrix dataframe.
#'
#' @param data_venn dataframe.input the genelist dataframe
#' @param getlogic TRUE or FALSE.Default:FALSE.if set true, the output will be logic dataframe.
#' @return dataframe. Default:contain 0 or 1 dataframe.
#' @export
#'
#' @examples
#'load("../data/veen.Rdata")
#' head(data_venn)
#'   PH6WC_XY335   PH4CVvsXY335   PH4CVvsPH6WC         RIL335       F2_3_335       F2_3_958         RIL958
#' 1 Zm00001d033896 Zm00001d044327 Zm00001d033896 Zm00001d053675 Zm00001d031899 Zm00001d034543 Zm00001d018779
#' 2 Zm00001d032359 Zm00001d012313 Zm00001d004348 Zm00001d043045 Zm00001d013795 Zm00001d045451 Zm00001d002000
#' 3 Zm00001d008977 Zm00001d035030 Zm00001d009840 Zm00001d047789 Zm00001d021310 Zm00001d021784 Zm00001d033091
#' 4 Zm00001d038546 Zm00001d013456 Zm00001d012313 Zm00001d018386 Zm00001d047789 Zm00001d047789 Zm00001d043299
#' get_matrix(data_venn)
get_matrix <- function(data_venn,getlogic = FALSE){
  nsets <- length(data_venn) #获取数据的列数
  dim(data_venn)
  #Reduce是批量执行函数.Reduce返回的是1个向量
  alldata <- Reduce(union,data_venn %>% select(1:nsets)) #所有列求并集并去重
  alldata <- sort(alldata) #排序
  alldata[which(alldata == "")] <- NA # 修改“”的为NA
  alldata <- na.omit(alldata) #过滤掉NA的行
  #alldata 是所有列的基因的合并，去重后的向量
  #使用str_sort函数对数据框的指定列进行逐列分别排序。Map是批量执行函数，返回的是list.
  data_venn_sort <- data.frame(Map(str_sort,data_venn %>% select(1:nsets)))

  #输入基因列（genelist)，返回所有基因中是否有输入的基因的1/0向量，
  ##返回值长度是所有基因的长度。1表示有，0表示没有。
  if(getlogic){
    gene2num <- function(genelist){
      return((is.element(alldata,genelist))) #如果需要返回逻辑矩阵，则此处就不转为数字矩阵
    }
  }else{
    gene2num <- function(genelist){
      return(as.numeric((is.element(alldata,genelist)))) #默认：返回的是0/1矩阵
    }
  }

  #获取最终的0/1矩阵。
  data_venn_logic <-
    data.frame(Map(gene2num,data_venn_sort %>% select(1:nsets)))  #获取从gene到0/1矩阵
  rownames(data_venn_logic) <- alldata
  #view(data_venn_logic)
  return(data_venn_logic)  #返回0/1矩阵
}



##################获取2-4元的venn图#########
#' Title:get_venn: input a dataframe ,output venn diagram.Only draw venn graphs of 2-4 dimensions
#'
#' @param data2
#' @param fill_color Default:color , can be color 1 to color4, or a color vector of length 4
#' @param percentage logic:TRUE or FALSE, control show the percentage .
#' @return
#' @export
#'
#' @examples
#' load("../data/veen.Rdata)
#' head(data2)
#'  PH6WC_XY335 PH4CVvsXY335 PH4CVvsPH6WC RIL335 F2_3_335
#' ZeamMp007           0            0            0      1        0
#' ZeamMp010           0            0            0      0        0
#' ZeamMp011           0            0            0      0        0
#'
#'
get_venn <- function(data2,fill_color = color,percentage = TRUE){
  #设置删除geneID的列，同时转为逻辑矩阵
  dat3 <- data.frame(Map(as.logical,data2))
  #select返回的值是list,rowname需要的是向量.所以用unlist
  rownames(dat3) <- rownames(data2)
  #view(dat3)
  head(dat3)
  color <- c("red","green","orange","pink","purple")
  if(!require(RColorBrewer))install.packages("RColorBrewer")
  library(RColorBrewer)
  color1 <- brewer.pal(4,"Set1")
  color2 <- brewer.pal(4,"Set2")
  color3 <- brewer.pal(4,"Dark2")
  color4 <- brewer.pal(4,"Accent")

  if(length(dat3)>4){
    dat3 <- select(dat3,1:4)
    warning("the input dataframe columns are more than 4! Only show the first 4 columns")
  }
  img <- ggvenn(dat3,set_name_size = 4,text_size = 4,fill_color = fill_color,fill_alpha=0.6,show_percentage = percentage)
  return(img)
}


#########直接从gene到输出venn图和upset图####
#' Title
#'
#' @param data_gene dataframe. contain the columns of each group.
#' @param outname string.Default:"demo",for output file name prefix.
#'
#' @return
#' @export
#'
#' @examples
#' data_venn
#' head(data_venn)
#'   PH6WC_XY335   PH4CVvsXY335   PH4CVvsPH6WC         RIL335       F2_3_335       F2_3_958         RIL958
#' 1 Zm00001d033896 Zm00001d044327 Zm00001d033896 Zm00001d053675 Zm00001d031899 Zm00001d034543 Zm00001d018779
#' 2 Zm00001d032359 Zm00001d012313 Zm00001d004348 Zm00001d043045 Zm00001d013795 Zm00001d045451 Zm00001d002000
#' 3 Zm00001d008977 Zm00001d035030 Zm00001d009840 Zm00001d047789 Zm00001d021310 Zm00001d021784 Zm00001d033091
#' 4 Zm00001d038546 Zm00001d013456 Zm00001d012313 Zm00001d018386 Zm00001d047789 Zm00001d047789 Zm00001d043299
#' gene2venn(data_venn)
gene2venn <- function(data_gene,outname="demo"){

  #获取矩阵
  data2 <- get_matrix(data_gene)
  csvfile <- paste0(outname,"_matrix.csv")
  write.csv(data2,file=csvfile,row.names = TRUE)
  wdpath <- getwd()
  message("output gene matrix in :",paste(wdpath,csvfile,sep = "/"))

  #获取venn图
  vennfile <- paste0(outname,"_venn.pdf")
  p1 <- get_venn(data2)
  #ggsave(filename = paste0(vennfile,".tiff"),p1)
  ggsave(filename = vennfile,p1)
  message("output venn image in :",paste(wdpath,vennfile,sep = "/"))

  #获取upset图
  data2$geneID <- rownames(data2)
  upsetfile <- paste0(outname,"_upset.pdf")
  pdf(paste0(upsetfile))
  p2 <- upset(data2,decreasing = c(FALSE,TRUE))
  print(p2) #使用print绘制，否则是空白
  dev.off()
  message("output upset image in :",paste(wdpath,upsetfile,sep = "/"))

  #获取交集基因
  intersectgene <-
    data2 %>%
    select(-(geneID)) %>%  #删除geneID列
    filter_all(all_vars(. == 1)) %>%  #获取每一列都是1的列
    rownames #获取行名，即为最终的交集基因
  return(intersectgene) #返回的是所有列的基因的交集
  #return(data2) #返回0/1矩阵
}

#############demo1 示例############
#############demo1 示例############
demo1 <- function(){
  load("../data/veen.Rdata")
  intersect_gene <- gene2venn(data_venn %>% select(1:4)) #获取所有列的交集
  write.table(intersect_gene,file = "intersect_gene.csv",row.names = FALSE)
}

demo1()


###############demo2示例################
demo2 <- function(){
  if(!require("ggpolypath"))install.packages("ggpolypath")
  require(ggplot2)
  require(ggpolypath)
  require(venn)
  require(dplyr)

  load("../data/veen.Rdata") #读取本地的示例数据
  data1 <- data_venn %>% select(1:6)
  listdata1 <- get_matrix(data1) %>% select(2:4)
  venn(listdata1, ilabels = TRUE, zcolor = "style",ilcs = 0.8,
       ggplot = TRUE,borders = FALSE,box = FALSE,
       par = FALSE,plotsize = 30)
  ggsave("demo2.tiff")
}
demo2()



