#get_volcano
#setwd("E:/pakwork/cancer/")
if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
if (!requireNamespace("ggplot2", quietly = TRUE))install.packages("ggplot2")
if(!require(ggprism))remotes::install_github("csdaw/ggprism")
library(ggprism) #可以完善ggplot2的图使之达到发表级别

##输入的resdata 第一列必须是foldchange,第二列必须是pvalue
#定义生成火山图的函数，参数：read_count矩阵，类型1，类型2，差异倍数(默认是1)，x轴范围（默认是-10,10），y轴范围（默认是0，100），
#' Title get a volcano image.
#'
#' @param resdata a dataframe , fist col is Foldchange,second col is pvalue.
#' @param foldchange Default:1
#' @param xrange Default:c(-10:10)
#' @param yrange Default:c(0,100)
#' @param legends Ture|FALSE Default:FALSE
#'
#' @return volcano
#' @export the img of volcano
#'
#' @examples
#' get_volcano(testdata,foldchange=1,xrange=c(-10,10),yrange=c(0,120),legends=TRUE)
#' get_volcano(testdata,foldchange=1.5)
#' get_volcano(testdata)
get_volcano <- function(resdata,foldchange=1,xrange=c(-10,10),yrange=c(0,100),legends=FALSE){
  #删除含有Na的行
  resdata <- na.omit(resdata)
  colnames(resdata)[1:2] <- c("log2FoldChange","padj") #指定第1,2列的列名分别是差异倍数和padj
  threshold <- as.factor(ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange) >= foldchange ,ifelse(resdata$log2FoldChange >= foldchange ,'Up','Down'),'Not'))
  deg_img <- ggplot(resdata,aes(x=resdata$log2FoldChange,y=-log10(resdata$padj),colour=threshold)) +
    xlab("log2(Fold Change)")+ylab("-log10(padj)") +
    geom_point(size = 0.5,alpha=1) + ylim(yrange[1],yrange[2]) + xlim(xrange[1],xrange[2]) +
    scale_color_manual(values=c("green","grey", "red"))+
    #scale_shape_prism() +
    theme_prism(base_size=10)
    #theme(legend.position = "none")
  bold_deg_img <- deg_img+scale_y_continuous(guide = "prism_offset_minor",limits = yrange)+ scale_x_continuous(guide = "prism_offset_minor",limits = xrange)
  ##添加阈值线
  line_valcao <- bold_deg_img+
    #geom_vline(xintercept=c(-foldchange,foldchange),colour="#990000",linetype="dashed")+
    #geom_hline(yintercept=1.3,colour="#0000FF",linetype="dashed")
    geom_hline(yintercept = 1.3,linetype="dotted")+geom_vline(xintercept=c(-foldchange,foldchange),linetype="dotted")
  ##如果lengends是true,则有图例，默认是无图例
  if (!legends){
    line_valcao <- line_valcao+theme(legend.position = "none")
  }

  line_valcao
  # ggsave(paste(type1,"_",type2,".deg.pdf",sep=""),line_valcao)
  # ggsave(paste(type1,"_",type2,".deg.tiff",sep=""),line_valcao)
  return(line_valcao)
}

devtools::document()


file.edit("DESCRIPTION")
