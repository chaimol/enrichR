#require(ggord)
#require(ggplot2)
## dat的输入矩阵是，行是基因名，列是材料/样本。group_list是材料/样本的分组信息的向量list
#' get_pca
#' @import ggord
#' @import ggplot2
#' @param dat a dataframe; row is samples/materials,column is genename
#' @param group_list a list: contain the samples/materials repeat.
#' group_list should like this:`"RC"  "RC"  "RC"  "RT1" "RT1" "RT1" "RT2" "RT2" "RT2" "RT3" "RT3" "RT3"`
#' @return a ggplot2 image.
#' @export getPCA
#'
#' @examples
#' load("PCA.Rdata")
#' group_level <- substr(colnames(read_counts), 1, 3)
#‘ group_level[1:3] <- rep("RC", 3)
#’ group_level[16:18] <- rep("SC", 3)
#' getPCA(read_counts,group_level)
getPCA <- function(dat, group_list) {
  dat <- dat[which(rowSums(dat) > 0), ]
  dat <- t(dat) # 转置
  pca1 <- prcomp(dat, scale. = T, center = TRUE, retx = T) # prcomp是R自带的pca分析函数
  p1 <- ggord(pca1, grp_in = group_list, arrow = 0, vec_ext = 0, txt = NULL) + theme_classic()
  return(p1)
}

# dat格式如下
# RC1_readcount RC2_readcount RC3_readcount RT1_1_readcount RT1_2_readcount RT1_3_readcount RT2_1_readcount
# MELO3C000001             0             0             0               0               0               0               0
# MELO3C000002             0             0             0               0               0               0               0
# MELO3C000003          2310          2674          2525            2059            1992            2136            1791

# group_list格式如下
# "RC"  "RC"  "RC"  "RT1" "RT1" "RT1" "RT2"
demoPCA <- function() {
  all_dat <- read.delim("e:/bioinformation_center/tiangua/updatedata/all_readcount.xls", row.names = 1, header = T)
  all_dat <- all_dat[1:500, ]
  group_level <- substr(colnames(all_dat), 1, 3)
  group_level[1:3] <- rep("RC", 3)
  group_level[16:18] <- rep("SC", 3)
  source("E:/pakwork/cancer/pca.R")
  p2 <- getPCA(all_dat, group_level)
  ggsave("PCA.pdf",p2)
}
