if(!require(ggplot2))BiocManager::install("ggplot2")
if(!require(devtools))install.packages("devtools")
library(devtools)
if(!require(ggord))install_github('fawda123/ggord')


##dat的输入矩阵是，行是基因名，列是材料/样本。group_list是材料/样本的分组信息的向量list
#' Title
#'
#' @param dat 
#' @param group_list 
#'
#' @return
#' @export
#'
#' @examples
get_pca <- function(dat,group_list){
  dat <- dat[which(rowSums(dat) > 0),]
  dat <- t(dat) #转置
  pca1=prcomp(dat,scale.=T,center=TRUE,retx=T)#prcomp是R自带的pca分析函数
  ggord(pca1, grp_in = group_list, arrow=0, vec_ext =0,txt=NULL)+theme_classic()
}

#dat格式如下
#RC1_readcount RC2_readcount RC3_readcount RT1_1_readcount RT1_2_readcount RT1_3_readcount RT2_1_readcount
#MELO3C000001             0             0             0               0               0               0               0
#MELO3C000002             0             0             0               0               0               0               0
#MELO3C000003          2310          2674          2525            2059            1992            2136            1791

#group_list格式如下
# "RC"  "RC"  "RC"  "RT1" "RT1" "RT1" "RT2"
demo <- function(){
  all_dat <- read.delim("e:/bioinformation_center/tiangua/updatedata/all_readcount.xls",row.names = 1,header = T)
  all_dat <- all_dat[1:500,]
  group_level <- substr(colnames(all_dat),1,3)
  group_level[1:3] <- rep("RC",3)
  group_level[16:18] <- rep("SC",3)
  source("E:/pakwork/cancer/pca.R")
  get_pca(all_dat,group_level)
}
