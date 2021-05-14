### a function for barplot GO,KEGG results.

if(!require(ggplot2))install.packages("ggplot2")
if(!require(tidyverse))install.packages("tidyverse")
if(!require(ggprism))remotes::install_github("csdaw/ggprism")
library(ggplot2)
library(tidyverse)
library(ggprism) #可以完善ggplot2的图使之达到发表级别


#view(data_GO_all)
#head(data_GO_all)
# ID                    Description count type
# 1 GO:0046906           tetrapyrrole binding    60   MF
# 2 GO:0015979                 photosynthesis    31   BP
# 3 GO:0009535 chloroplast thylakoid membrane    26   BP
# 4 GO:0055035     plastid thylakoid membrane    26   CC
# 5 GO:0044436                 thylakoid part    31   CC
# 6 GO:0042651             thylakoid membrane    28   MF

##可视化GO和KEGG的富集结果.输入的矩阵如上所示，要求列名必须完全一致。
#' Title barplot: a function for GO,KEGG enrich results visual in barplot
#'
#' @description  a function for GO,KEGG or TF enrich results visual in barplot
#' @param data_GO_all a dataframe like this
#' head(data_GO_all)
# ID                    Description count type
# 1 GO:0046906           tetrapyrrole binding    60   MF
# 2 GO:0015979                 photosynthesis    31   BP
# 3 GO:0009535 chloroplast thylakoid membrane    26   BP
# 4 GO:0055035     plastid thylakoid membrane    26   CC
# 5 GO:0044436                 thylakoid part    31   CC
# 6 GO:0042651             thylakoid membrane    28   MF
#' @param fill_level vector. for show the terms' type.Default:c("MF","BP","CC","KEGG")
#' if you also enrich TF or other.you can use c("MF","BP","CC","KEGG","TF").
#' this vector's string must be the same with the input dataframe type row.
#' @param term_num numeric. for show terms count number。Default:8
#' @param term_width numeric. for show terms' description width.Default:40
#'
#' @return a ggplot2 results
#' @export
#'
#' @examples
#' dat1 <- read.csv("./data/enrich_result.csv")
#' barplot(dat1,term_num = 4)
#' barplot(dat1)
#' barplot(dat1,term_num = 5,term_width = 45)
#'
#'
#'
barplot <- function(data_GO_all,fill_level=c("MF","BP","CC","KEGG"),term_num=8,term_width=40){
  data1 <- data_GO_all %>% arrange(type,desc(count)) %>% group_by(type) %>%
    slice(1:term_num) #筛选出前n项用于展示
  #分组排序
  data_GO_all <- group_by(data1,type) %>%
    arrange(desc(type),count)
  #根据分组排序的结果进行factor。用于固定顺序，绘图时按照这个顺序绘图
  data_GO_all$type <- factor(data_GO_all$type,levels = fill_level)
  data_GO_all$Description <- factor(data_GO_all$Description,levels = data_GO_all$Description)

  p1 <- ggplot(data_GO_all) +
    geom_bar(stat="identity",aes(x=count,y=Description,fill=type))+
    #facet_wrap(~group)+ #设置分面
    labs(y='Term',x='Count',fill='Type')+ #设置标题
    #coord_flip() + 颠倒xy轴
    scale_y_discrete(labels=function(x) str_wrap(x, width=term_width))+ #设置term的description的宽度，自动换行
    theme_prism(base_size=10) #使用prism标准化图，达到出版级别
  p1
  return(p1)
}


