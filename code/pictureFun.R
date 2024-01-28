# (1) allrelationFun ------------------------------------------------------
allrelationFun = function(subtree,p,realpath,realresultpath){
  newdata = read_csv(file=paste(realpath,"mer_",subtree,"_",maxdelta_time,"_selectdata.csv",sep=""),show_col_types = FALSE)
  finalinfo = newdata[,1:6]
  finaldata = newdata[,-c(1:6)]
  genename = names(finaldata)
  connect = data.frame(from = "a",
                       to = "b",
                       preaij = 0,
                       predelta = 0)
  load(file=paste(realpath,subtree,"_simu_info.RData",sep=""))
  Mc = simu_info$Mc
  eachdelta = 1/(round(Mc)-1)
  preaij = read.csv(file=paste(resultpath,subtree,"_",p,"_preaij.csv",sep=""))
  predelta = read.csv(file=paste(resultpath,subtree,"_",p,"_predelta.csv",sep=""))
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
  return(connect)
}

# (2) networkFun ----------------------------------------------------------
networkFun = function(subtree,connect,maxdelta_time,realpath){
  newdata = read_csv(file=paste(realpath,"mer_",subtree,"_",maxdelta_time,"_selectdata.csv",sep=""),show_col_types = FALSE)
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
  angle= 360 * (vertices$id-0.5) /number_of_bar
  vertices$hjust<-ifelse(angle>180, 1, 0)
  vertices$angle<-ifelse(angle>180, 90-angle+180, 90-angle)
  
  mygraph <- graph_from_data_frame( connect, vertices = vertices, directed = TRUE )
  mycolor = brewer.pal(8, 'Dark2')[c(1:6)]
  g = ggraph(mygraph,layout = "linear",circular=T) +
    geom_edge_arc(aes(edge_colour=factor(predelta,levels=c(0,1,2,3,4,5)),
                      linetype = factor(preaij,levels=c(1,-1))), edge_alpha=0.8, 
                  edge_width=0.6,
                  arrow = arrow(length = unit(1.5, "mm"),type="open",ends= "last"),
                  start_cap = square(1.5, 'mm'),
                  end_cap = circle(1.5, 'mm'))+
    scale_edge_linetype_manual(values = c("solid","dashed"),labels = c("Positive","Negative"),drop=F)+
    scale_edge_color_manual(values=mycolor,
                            labels=c("0",expression(Delta*t),expression(2*Delta*t),
                                     expression(3*Delta*t),expression(4*Delta*t),
                                     expression(5*Delta*t)),drop=F) +
    geom_node_text(aes(x = x*1.16, y=y*1.16, label=name, angle=angle,hjust=hjust),size=3,fontface='bold')+
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
      legend.title = element_text(colour='black', size = 14,face="bold"),
      legend.text = element_text(colour="black", size = 12,face="bold",hjust = 0),
      legend.direction = 'vertical',
      legend.box.background = element_rect(fill=NA,color = "black",linetype = 1),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks =element_blank(),
      axis.text =element_blank(),
      axis.title = element_blank(),
      plot.margin=unit(c(0,0,0,0), "null"),
      panel.spacing=unit(c(0,0,0,0), "null"))+
    annotate("text", x = 0, y = 0,
             label = paste(subtree),
             colour = "black",size=10,fontface = "bold")
  return(g)
}