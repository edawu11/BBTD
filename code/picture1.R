getwd()
library(tidyverse)
library(igraph)
library(ggraph)
library(RColorBrewer)
library(ggtext)
library(lemon)
library(ggpubr)
library(caret)
library(cowplot)
allbosscell = c("AB","C","D","E","MS")
#####our method#####
tempdata1 = data.frame(TPR=rep(0,5),
                       PPR=rep(0,5),
                       subtree=rep("AB",5),
                       method=rep("BBTD",5))
datanum = 50
for(k in 1:5){
  bosscell = allbosscell[k]
  dir_simu = paste("./simu_data_new/",bosscell,"/",sep="")
  dir_result = paste("./simu_result_new/1127/",bosscell,"/",sep="")
  allresult = matrix(0,datanum,2)
  for(p in 1:datanum){
    realaijmat = read.csv(file = paste(dir_simu,"realaijmat_",p,"_",bosscell,"_.csv",sep=""))
    realaijmat = abs(as.matrix(realaijmat))
    preaijmat = read.csv(file = paste(dir_result,bosscell,"_",p,"_preaij.csv",sep=""))
    preaijmat = abs(as.matrix(preaijmat))
    conmat = table(Prediction = factor(as.numeric(preaijmat),levels = c(1,0),labels = c(1,0)),
                   Reference = factor(as.numeric(realaijmat),levels = c(1,0),labels = c(1,0)))
    
    TPR = conmat[1,1]/sum(conmat[,1])
    PPR = conmat[1,1]/sum(conmat[1,])
    allresult[p,1] = TPR
    allresult[p,2] = PPR
  }
  tempdata1[k,1:2] = colMeans(allresult)
  tempdata1[k,3] = bosscell
}

tempdata2 = data.frame(TPR=rep(0,5),
                       PPR=rep(0,5),
                       subtree=rep("AB",5),
                       method=rep("BIBN",5))
for(k in 1:5){
  bosscell = allbosscell[k]
  dir_comp = paste("./simu_result_new/1127/compare_result/",bosscell,"/",sep="")
  thisresult = read.csv(file = paste(dir_comp,bosscell,"_bestfit.csv",sep=""))
  tempdata2[k,1:2] = thisresult[1,]
  tempdata2[k,3] = bosscell
}


tempdata3 = data.frame(TPR=rep(0,5),
                       PPR=rep(0,5),
                       subtree=rep("AB",5),
                       method=rep("BFE",5))
for(k in 1:5){
  bosscell = allbosscell[k]
  dir_comp = paste("./simu_result_new/1127/compare_result/",bosscell,"/",sep="")
  thisresult = read.csv(file = paste(dir_comp,bosscell,"_reveal.csv",sep=""))
  tempdata3[k,1:2] = thisresult[1,]
  tempdata3[k,3] = bosscell
}

tempdata4 = data.frame(TPR=rep(0,5),
                       PPR=rep(0,5),
                       subtree=rep("AB",5),
                       method=rep("ATEN",5))
for(k in 1:5){
  bosscell = allbosscell[k]
  dir_comp = paste("./simu_result_new/1127/compare_result/",bosscell,"/",sep="")
  thisresult = read.csv(file = paste(dir_comp,bosscell,"_reveal.csv",sep=""))
  tempdata4[k,1:2] = thisresult[1,]
  tempdata4[k,3] = bosscell
}

testdata = rbind(tempdata1,tempdata2,tempdata3,tempdata4)



# testdata = data.frame(TPR = rep(c(0.9,0.8,0.3,0.2),5),
#                       PPR = rep(c(0.9,0.8,0.3,0.2),5),
#                       subtree = factor(rep(c("AB","C","D","E","MS"),each=4),
#                                        levels = c("AB","C","D","E","MS")),
#                       method = factor(rep(c("BIBNTD","BIBN","BFE","ATEN"),5),
#                                        levels = c("BIBNTD","BIBN","BFE","ATEN")))
####methods comparison####
testdata = read.csv(file="testdata.csv")
testdata =gather(testdata,key = "index",value = "value",-subtree,-Method) %>% add_column(temp = "aaa")
testdata$index = factor(testdata$index,levels = c("TPR","PPR"))
testdata$subtree = factor(testdata$subtree,levels = c("AB","C","D","E","MS"))
# testdata$method = factor(testdata$method,levels = c("BIBNTD","BIBN","BFE","ATEN"))
testdata$Method = factor(testdata$Method,levels = c("ATEN","BFE","BB","BBTD"))
ggplot(data=testdata,aes(x =value,y= temp,fill=Method))+
  geom_col(position = 'dodge',width = 0.6)+
  guides(fill = guide_legend(reverse = TRUE))+
  scale_fill_manual(values=c("#F8E71C", "#A0E8AF", "#FFC107","#4CAF50"))+
  facet_grid(rows = vars(subtree), cols = vars(index))+
  scale_x_continuous(breaks = seq(0,1,by=0.2))+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks=element_blank())
ggsave(file = paste("./real_result_nona/picture/compare.pdf",sep=""),
       width = 14,
       height = 8,
       units = "cm")

####transform data####
allrelationfun = function(bosscell,p){
  realpath = paste("./real_data/readydata/",bosscell,"/",sep="")
  resultpath = paste("./real_result_seed5/",bosscell,"/",sep="")
  newdata = read_csv(file=paste(realpath,"mer_",bosscell,"_",maxdelta_time,"_selectdata.csv",sep=""),show_col_types = FALSE)
  finalinfo = newdata[,1:6]
  finaldata = newdata[,-c(1:6)]
  genename = names(finaldata)
  connect = data.frame(from = "a",
                       to = "b",
                       preaij = 0,
                       predelta = 0)
  load(file=paste(realpath,bosscell,"_simu_info.RData",sep=""))
  Mc = simu_info$Mc
  eachdelta = 1/(round(Mc)-1)
  preaij = read.csv(file=paste(resultpath,bosscell,"_",p,"_preaij.csv",sep=""))
  predelta = read.csv(file=paste(resultpath,bosscell,"_",p,"_predelta.csv",sep=""))
  for(i in 1:length(genename)){
    nozero = which(preaij[i,]!=0)
    if(length(nozero)==0){
      next
    }
    else{
      for(j in nozero){
        if(preaij[i,j]==-0.1){
          preaij[i,j]=-1
        }
        tempconnect = data.frame(from = genename[j],
                                 to = genename[i],
                                 preaij = preaij[i,j],
                                 predelta = round(predelta[i,j]/eachdelta))
        connect = rbind(connect,tempconnect)
      }
    }
  }
  connect = connect[-1,]
  print(sort(unique(connect$predelta)))
  return(connect)
}

# geom_edge_loop(aes(edge_colour= factor(predelta,levels=c(0,1,2,3,4,5)),
#                      linetype = factor(preaij,levels=c(1,-1,-0.1)),
#                      direction =ifelse((from - 0.5) * 360 / length(mygraph)>180,
#                                        -90-(from - 0.5) * 360 / length(mygraph)+180,
#                                        90-(from - 0.5) * 360 / length(mygraph))),
#                  edge_alpha=0.8, edge_width=0.6,
#                  arrow = arrow(length = unit(4, "mm"),type="open",ends= "last"),
#                  start_cap = square(3, 'mm'),
#                  end_cap = circle(3, 'mm'))

####network results####
pic_fun = function(bosscell,connect,indexnum){
  maxdelta_time =5
  realpath = paste("./real_data/readydata/",bosscell,"/",sep="")
  resultpath = paste("./real_result_seed5/",bosscell,"/",sep="")
  picpath = paste("./real_result_nona/picture/",bosscell,"/",sep="")
  
  newdata = read_csv(file=paste(realpath,"mer_",bosscell,"_",maxdelta_time,"_selectdata.csv",sep=""),show_col_types = FALSE)
  finalinfo = newdata[,1:6]
  finaldata = newdata[,-c(1:6)]
  genename = names(finaldata)
  
  c( as.character(connect$from), as.character(connect$to)) %>%
    as_tibble() %>%
    group_by(value) %>%
    summarize(n=n()) -> vertices
  colnames(vertices) <- c("name", "n")
  
  connect <- connect %>%
    filter(from %in% vertices$name) %>%
    filter(to %in% vertices$name) %>%
    left_join(vertices,by=c('from'='name'))
  
  number_of_bar<-nrow(vertices)
  vertices$id = seq(1, nrow(vertices))
  angle= 360 * (vertices$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  vertices$hjust<-ifelse(angle>180, 1, 0)
  vertices$angle<-ifelse(angle>180, 90-angle+180, 90-angle)
  
  mygraph <- graph_from_data_frame( connect, vertices = vertices, directed = TRUE )
  # mycolor = brewer.pal(12, 'Paired')[c(2,4,6,8,10,12)]
  mycolor = brewer.pal(8, 'Dark2')[c(1:6)]
  g1 = ggraph(mygraph,layout = "linear",circular=T) +
    geom_edge_arc(aes(edge_colour=factor(predelta,levels=c(0,1,2,3,4,5)),
                      linetype = factor(preaij,levels=c(1,-1))), edge_alpha=0.8, 
                  edge_width=0.6,
                  arrow = arrow(length = unit(4, "mm"),type="open",ends= "last"),
                  start_cap = square(3, 'mm'),
                  end_cap = circle(3, 'mm'))+
    scale_edge_linetype_manual(values = c("solid","dashed"),labels = c("Positive","Negative"),drop=F)+
    scale_edge_color_manual(values=mycolor,
                            labels=c("0",expression(Delta*t),expression(2*Delta*t),
                                     expression(3*Delta*t),expression(4*Delta*t),
                                     expression(5*Delta*t)),drop=F) +
    geom_node_text(aes(x = x*1.16, y=y*1.16, label=name, angle=angle,hjust=hjust),size=5,fontface='bold')+
    geom_node_point(aes(shape = 1), fill = "#ffa400",shape=21,size=3,
                    color='#ffa400',alpha=0.6) +
    scale_size_continuous(range=c(0.5,10)) +
    scale_fill_manual(values=mycolor) +
    scale_color_manual(values=mycolor) +
    
    expand_limits(x = c(-1.6, 1.6), y = c(-1.6, 1.6))+
    coord_fixed()+
    theme_minimal() +
    labs(edge_colour="Time delay",
         edge_linetype="Regulation types"
    )+
    theme(
      legend.title = element_text(colour='black', size = 18,face="bold"),
      legend.text = element_text(margin = margin(l = 10),colour="black", size = 15,face="bold",hjust = 0),
      legend.direction = 'vertical',
      legend.box.background = element_rect(fill=NA,color = "black",linetype = 1),
      legend.key.width = unit(20,"pt"),
      legend.key.height = unit(20, "pt"),
      legend.margin = margin(t = 10, r = 15, b = 10, l = 15),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks =element_blank(),
      axis.text =element_blank(),
      axis.title = element_blank(),
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null"))+
    annotate("text", x = 0, y = 0,
             label = paste(bosscell),
             colour = "black",size=10,fontface = "bold")+
    annotate("richtext", x = -2, y = 1.2,
             label = indexnum,
             colour = "black",size=10,
             label.colour = NA,
             fill = NA,fontface = "bold")
  
  return(g1)
}

allbosscell = c("AB","C","D","E","MS")
allindex = c("a","b","c","d","e")
allseed = c(1,1,1,1,1)
alltreedata = data.frame(from="a",
                         to="b",
                         preaij=1,
                         predelta=1,
                         bosscell="AB")
maxdelta_time= 5
for(i in 1:5){
  bosscell = allbosscell[i]
  p = allseed[i]
  resultpath = paste("./real_result_self/",bosscell,"/",sep="")
  connect = allrelationfun(bosscell,p)
  write.csv(connect,file = paste(resultpath,bosscell,"_allrelation.csv",sep=""),row.names = F)
  connect= cbind(connect,data.frame(bosscell =bosscell))
  alltreedata = rbind(alltreedata,connect)
}

write.csv(alltreedata[-1,],file = "realrelationships.csv",row.names = F)
getwd()

for(i in 1:5){
  bosscell = allbosscell[i]
  thisindex = allindex[i]
  p = allseed[i]
  connect = allrelationfun(bosscell,p)
  print(bosscell)
  print(connect)
  assign(paste("g",i,sep=""),pic_fun(bosscell,connect,thisindex) +theme(legend.position = "none"))
}
legend <- g_legend(g1 + theme(legend.position='right'))
ggarrange(g1,g2,g3,g4,g5,legend)
# plot_grid(g1,g2,g3,g4,g5,legend,labels = c('a','b','c','d','e'))
ggsave(file = paste("./real_result_nona/picture/network.pdf",sep=""),
       width = 52,
       height = 30,
       units = "cm")


# convergence -------------------------------------------------------------
getwd()
allbosscell = c("AB","C","D","E","MS")
for(k in 1:5){
  bosscell = allbosscell[k]
  realpath = paste0("./real_result_seed5/",bosscell,"/")
  load(paste0(realpath,bosscell,"_1_allprob.RData"))
  ####transform data####
  allprobdf = data.frame(index= 1:2e4,
                         Chain1 = rep(0,2e4),
                         Chain2 = rep(0,2e4),
                         Chain3 = rep(0,2e4),
                         Chain4 = rep(0,2e4),
                         Chain5 = rep(0,2e4))
  for(i in 1:5){
    allprobdf[,i+1] = allprob[[i]]
  }
  p = ggplot()+geom_line(data = allprobdf,aes(x=index,y=Chain1,colour="Chain1"))+
    geom_line(data = allprobdf,aes(x=index,y=Chain2,colour="Chain2"))+
    geom_line(data = allprobdf,aes(x=index,y=Chain3,colour="Chain3"))+
    geom_line(data = allprobdf,aes(x=index,y=Chain4,colour="Chain4"))+
    geom_line(data = allprobdf,aes(x=index,y=Chain5,colour="Chain5"))+
    xlab("Iteration index")+ylab("Log-posterior")+
    theme_bw()+
    theme(legend.position = 'none',
          axis.text.x = element_text(size=8),
          axis.text.y = element_text(size=8))
  assign(paste0("p",k),p)
}
legend <- g_legend(p1 + theme(legend.position='right',
                              legend.title = element_blank(),
                              legend.text = element_text(size =12)))
plot_grid(p1,p2,p3,p4,p5,legend,labels = c('a','b','c','d','e'))
ggsave(file = paste("./real_result_nona/picture/conven.pdf",sep=""),
       width = 22,
       height = 12,
       units = "cm")

# missing -----------------------------------------------------------------
library(visdat)
library(cowplot)

napicFun = function(usedata,xname,yname){
  allnum = nrow(usedata)*(ncol(usedata)-1)
  nanum = sum(unlist(usedata[,-1])==0)
  naratio = round(nanum/allnum*100,2)
  nonaratio = 100-naratio
  usedata2 = usedata %>% gather(key = gene,value = ifna,-cell)
  usedata2$ifna = factor(usedata2$ifna)
  ggplot(data = usedata2, aes(x=gene, y=cell, fill=ifna)) + 
    geom_tile()+
    scale_fill_manual(values = c("#A7A9AC", "#FFE066"),
                      labels = c(paste0("NA"," ","(",naratio,"%)"),paste0("Not NA"," ","(",nonaratio,"%)")))+
    labs(x = xname,y=yname)+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size =15),
          axis.title.y = element_text(size =15),
          legend.title = element_blank(),
          legend.position = "bottom" ,
          legend.text= element_text(size = 12))
}

maxdelta_time = 5
allbosscell = c("AB","C","D","E","MS")
for(i in 1:5){
  bosscell = allbosscell[i]
  realpath = paste("E:/xmusta/R/GRN/code/usealldata_com/real_data/readydata/",bosscell,"/",sep="")
  newstatdata = read_csv(file=paste(realpath,"mer_",bosscell,"_",maxdelta_time,"_statdata.csv",sep=""),
                         show_col_types = F)
  assign(paste("g",i,sep=""), napicFun(newstatdata,"Gene","Cell"))
}
plot_grid(g1,g2,g3,g4,g5,labels =c("a","b","c","d","e"),ncol=3,nrow=2,label_size = 16)
ggsave(file = paste("./picture/allmissing.pdf",sep=""),
       width = 12,
       height = 8)

maxdelta_time = 5
allbosscell = c("AB","C","D","E","MS")
for(i in 1:5){
  bosscell = allbosscell[i]
  realpath = paste("E:/xmusta/R/GRN/code/usealldata_com/real_data/readydata/",bosscell,"/",sep="")
  newstatdata = read_csv(file=paste(realpath,"mer_",bosscell,"_",maxdelta_time,"_statdata.csv",sep=""),
                         show_col_types = F)
  selectdata = read_csv(file=paste(realpath,"mer_",bosscell,"_",maxdelta_time,"_selectdata.csv",sep=""),show_col_types = FALSE)
  finalinfo = selectdata[,1:6]
  finaldata = selectdata[,-c(1:6)]
  selectgene = names(finaldata)
  selectcell = unique(finalinfo$cell)
  newstatdata = newstatdata %>% filter(cell %in% selectcell) %>% select(cell,all_of(selectgene))
  assign(paste("g",i,sep=""), napicFun(newstatdata,"Candidate gene","Candidate cell"))
}
plot_grid(g1,g2,g3,g4,g5,labels =c("a","b","c","d","e"),ncol=3,nrow=2,label_size = 16)
ggsave(file = paste("./picture/allcandidatemissing.pdf",sep=""),
       width = 12,
       height = 8)

onepicFun = function(usedata,xname,yname){
  allnum = nrow(usedata)*(ncol(usedata)-1)
  nanum = sum(unlist(usedata[,-1])==0)
  naratio = round(nanum/allnum*100,2)
  nonaratio = 100-naratio
  usedata2 = usedata %>% gather(key = gene,value = ifna,-originrow)
  usedata2$ifna = factor(usedata2$ifna)
  ggplot(data = usedata2, aes(x=gene, y=originrow, fill=ifna)) + 
    geom_tile()+
    scale_fill_manual(values = c("#A7A9AC", "#FFE066"),
                      labels = c(paste0("1"," ","(",naratio,"%)"),paste0("0"," ","(",nonaratio,"%)")))+
    labs(x = xname,y=yname)+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size =15),
          axis.title.y = element_text(size =15),
          legend.title = element_blank(),
          legend.position = "bottom" ,
          legend.text= element_text(size = 12))
}
maxdelta_time = 5
allbosscell = c("AB","C","D","E","MS")
for(i in 1:5){
  bosscell = allbosscell[i]
  realpath = paste("E:/xmusta/R/GRN/code/usealldata_com/real_data/readydata/",bosscell,"/",sep="")
  selectdata = read_csv(file=paste(realpath,"mer_",bosscell,"_",maxdelta_time,"_selectdata.csv",sep=""),show_col_types = FALSE)
  finaldata = selectdata %>% select(-cell,-time,-lastrow,-treeindex,-ifchange)
  finaldata = map_dfc(finaldata,function(x) replace_na(x,0))
  assign(paste("g",i,sep=""), onepicFun(finaldata,"Candidate gene","Row"))
}
plot_grid(g1,g2,g3,g4,g5,labels =c("a","b","c","d","e"),ncol=3,nrow=2,label_size = 16)
ggsave(file = paste("./picture/allcandidatemissing.pdf",sep=""),
       width = 12,
       height = 8)

# enrichment --------------------------------------------------------------
# data = read_csv(file="./real_result/cluster/alldata_merge.csv")
# data = data  %>% group_by(bosscell)
# data$Count = factor(data$Count,levels = c(2,3,4))
# temp1 = matrix(as.numeric(str_split(data$GeneRatio,"//",simplify = T)),ncol = 2) 
# changevalue = temp1[,1]/temp1[,2]
# 
# temp2 = matrix(as.numeric(str_split(data$BgRatio,"//",simplify = T)),ncol = 2) 
# changevalue2= temp2[,1]/temp2[,2]
# data$GeneRatio = changevalue
# data$BgRatio = changevalue2
# data$`Fold enrichment` = data$GeneRatio/data$BgRatio
# 
# names(data)[8] = "P-value"
# ggplot(data,aes(x=`Fold enrichment`,y=Description,color=`P-value`,shape=Count))+
#   geom_point(size=4)+
#   scale_color_gradient(low = "#E6F598",high = "#369790",trans="reverse")+
#   theme_bw()+
#   theme(
#     axis.title.y  = element_blank(),
#     strip.background = element_rect(
#       color="black", fill="white",linewidth = 1),
#     strip.text.x = element_text(size = 15, colour = "black",face="bold"),
#     axis.text.x = element_text(size = 12, color = "black", face = "bold"),
#     axis.text.y = element_text(size = 14, color = "black", face = "bold"),
#     axis.title.x = element_text(size = 18, color = "black", face = "bold"))+
#   facet_grid( .~ bosscell)
# 
# ggsave(file = paste("./code/usealldata_com/real_result_nona/picture/goenrich.pdf",sep=""),
#        width = 40, height = 25, units = "cm")

clusterpath = "./real_result/cluster/"
mergedata = read.csv(file=paste(clusterpath,"alldata_merge.csv",sep=""))
newmergedata = mergedata %>% group_by(bosscell,regulator,regulated) %>% 
  summarise(GO = paste0(ID,collapse = ", "),.groups = "drop" ) %>% arrange(bosscell)
newmergedata$regulator = paste0("\\emph{",newmergedata$regulator,"}")
newmergedata$regulated = paste0("\\emph{",newmergedata$regulated,"}")
write_csv(newmergedata,file=paste(clusterpath,"alldata_merge_new.csv",sep=""))

mergedata = read.csv(file=paste(clusterpath,"alldata_merge.csv",sep=""))
newdesdata = mergedata %>% dplyr::select(ID,Description) %>% distinct() %>% arrange(ID)
newdesdata$Description = str_to_sentence(newdesdata$Description)
write_csv(newdesdata,file=paste(clusterpath,"alldata_des_new.csv",sep=""))
