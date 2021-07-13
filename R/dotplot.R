#气泡图可视化富集结果
#' dotplot: Combine multiple groups of GO/KEGG/TF analysis to generate a bubble chart
#' @import ggplot2
#' @param data_GO_all  a dataframe,contain columns name of "type","Description","adjust","Count"
#' @param xlabels  a vector for x-axis labels
#' @param output  output file prefix
#' @param ranges  Default: c(4,7); a range,like c(3,5), for control the size of the bubble
#'
#' @return a ggplot object
#' @export dotplot
#'
#' @examples
#' load("data/data_GO_all.Rdata")
#' dotplot(data_GO_all,xlabels=xlabels,output="All_GO")
#' dotplot(data_GO_all,xlabels=xlabels,output="All_GO",range=c(3,5))
#' dotplot(data_GO_all,xlabels=xlabels,output="All_GO",range=c(4,7))
dotplot <- function(data_GO_all,xlabels=data_GO_all$type,output,ranges=c(4,7)){
  # filter Description string length over 60 will be split by ",",and extract the first sentence
  term <- c("term")
  for (i in 1:length(data_GO_all$Description)){
    if(nchar(data_GO_all$Description[i])>60){
      term <- append(term,data_GO_all$Description[i])
      #print(data_GO_all$Description[i])
      data_GO_all$Description[i] <- strsplit(data_GO_all$Description[i],split= ",")[[1]][1]
    }
  }
  #unique(term)
  #draw dot plot
  p1 <- ggplot(data_GO_all) +
  geom_point(aes(x=type,y=Description,fill = adjust,size=Count),alpha=0.85,pch=21)+
  #facet_wrap(~ group)+
  #fill对应点的填充色，colour对应点的边框色
  #scale_fill_gradient(low='Yellow',high='Red') + #设定颜色的变化范围
  scale_fill_gradient(low='SpringGreen',high='DeepPink') + #设定颜色的变化范围
  #scale_size_area() + #设定点的大小比例和图例上显示的间隔
  labs(y='GO term',x='Type',fill='P.adjust',size='Count number')+
  scale_size_continuous(range=ranges)+ #设置气泡的收缩比例，差值越大，最大气泡越大
  scale_x_discrete(labels= xlabels)+
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),axis.title.y = element_text(size = 7))
#+coord_flip()
  print(p1)
  ggsave(paste0(output,".tiff"),width = 170,height = 150,units = c("mm"),dpi=300)
  ggsave(paste0(output,".pdf"),width = 170,height = 230,units = c("mm"),dpi=300)
  write.csv(data_GO_all,file=paste0(output,".csv"),row.names = FALSE)
  work_path <- getwd()
  message("Output file is in ",work_path,"/",paste0(output,".csv"))
  message("Output file is in ",work_path,"/",paste0(output,".tiff"))
  message("Output file is in ",work_path,"/",paste0(output,".pdf"))
  return(p1)
}

