# (1) expressData ---------------------------------------------------------
expressData = function(fileindex,oldfilename,hufilename,newfilename,
                        hufilepath,oldfilepath,expdatapath){
  hufile = read_excel(path=paste(hufilepath,hufilename[fileindex],sep=""),sheet = 3)
  choosecell = hufile$cell
  choosecell = str_replace_all(choosecell,"\"","")
  choosetime = hufile$time
  choosetime = as.numeric(str_replace_all(choosetime,"\"","")) 
  oldfile = read_csv(file=paste(oldfilepath,oldfilename[fileindex],sep=""),show_col_types = FALSE)
  oldfile = oldfile %>% select(cell,time,blot)
  oldcell = unique(oldfile$cell)
  choosedata = oldfile %>% filter(cell %in% str_subset(oldcell,paste("^",choosecell[1],sep="")))
  for(k in 2:length(choosecell)){
    tempdata = oldfile %>% filter(cell %in% str_subset(oldcell,paste("^",choosecell[k],sep="")))
    choosedata = choosedata %>% bind_rows(tempdata)
  }
  choosedata = choosedata %>% mutate(ifblot = 1)
  for(k in 1:length(choosecell)){
    subdata = choosedata %>% filter(cell == choosecell[k]) %>% 
      mutate(ifblot = ifelse(time>=choosetime[k],1,0))
    chooseindex = which(choosedata$cell==choosecell[k])
    choosedata[chooseindex,] = subdata
  } 
  write_csv(choosedata,file = paste(expdatapath,newfilename[fileindex],sep=""))
}

# (2) cellnumStat ---------------------------------------------------------
cellnumStat = function(subtree,usepath,savepath){
  
  if(subtree == "E"){
    usesubtree = paste("E","[^M]|E",sep="")
  }
  else{
    usesubtree = subtree
  }
  
  allgenefile = list.files(path = usepath,
                           pattern = ".csv")
  alldata = list()
  genename = str_sub(allgenefile[1],start = 1,end = nchar(allgenefile[1])-4)
  genedata = read_csv(file = paste(usepath,allgenefile[1],sep = ""),show_col_types = FALSE)
  allchoosedata = genedata %>% filter(cell %in% str_subset(genedata$cell,paste("^",usesubtree,sep=""))) %>% 
    group_by(cell) %>% summarise(cellname = sum(ifblot))
  names(allchoosedata) = c("cell",genename)
  
  for(i in allgenefile[-1]){
    genedata = read_csv(file = paste(usepath,i,sep = ""),show_col_types = FALSE)
    genename = str_sub(i,start = 1,end = nchar(i)-4)
    choosedata = genedata %>% filter(cell %in% str_subset(genedata$cell,paste("^",usesubtree,sep="")))  %>% 
      group_by(cell) %>% summarise(cellname = sum(ifblot))
    names(choosedata) = c("cell",genename)
    allchoosedata  = allchoosedata %>% full_join(choosedata,by="cell")
  }
  
  write_csv(allchoosedata,file = paste(savepath,"all_",subtree,"_cell.csv",sep=""))
  
  allchoosedata = allchoosedata %>% select_if(~!all(is.na(.)))
  write_csv(allchoosedata,file = paste(savepath,"all_",subtree,"_nona_cell.csv",sep=""))
}


# (3) geneCopyFun ---------------------------------------------------------
geneCopyFun = function(subtree,usedata,allfilename,savepath){
  copyname = str_sub(allfilename[[3]],1,nchar(allfilename[[3]])-4)
  genename = allfilename[[4]]
  genecopy_h = hash(copyname,genename)
  
  newdata = usedata %>% arrange(cell)
  allcellname = newdata$cell
  alltreenum = rep(1,nrow(newdata))
  nowmother = allcellname[1]
  for(i in 2:nrow(newdata)){
    if(str_sub(allcellname[i],1,nchar(nowmother))==nowmother){
      alltreenum[i] = alltreenum[i-1]
    }
    else{
      nowmother = allcellname[i]
      alltreenum[i] = alltreenum[i-1]+1
    }
  }
  
  tree_h = hash(allcellname,alltreenum)
  
  simu_info = list()
  simu_info[["genecopy"]] = genecopy_h
  simu_info[["copyname"]] = keys(genecopy_h)
  simu_info[["genename"]] = sort(as.character(unique(values(genecopy_h))))
  simu_info[["cellname"]] = allcellname
  simu_info[["Mc"]] = median(unlist(newdata[,-1]),na.rm = T)
  save(simu_info,file=paste(savepath,subtree,"_simu_info.RData",sep=""))
  
}


# (4) interData -----------------------------------------------------------
interData = function(cellname,copyname,Mc,expdatapath,subtree,
                         maxdelta_time,savepath){
  M = round(Mc)
  eachdelta = 1/(M-1)
  maxdelta = maxdelta_time*eachdelta
  
  allfun = list()
  copynum = length(copyname)
  length(allfun) = copynum
  names(allfun) = copyname

  allcelllist = list()
  length(allcelllist) = length(cellname)
  names(allcelllist) = cellname

  for(i in copyname){
    allfun[[i]] = allcelllist
  }


  for(i in copyname){
    genefile = read.csv(paste(expdatapath,i,".csv",sep = ""))
    for(j in cellname){
      choosedata = genefile %>% filter(cell == j)
      if(nrow(choosedata) == 0){
        next
      }
      if(nrow(choosedata) < 3 & nrow(choosedata)>0){
        lastcell = str_sub(j,1,nchar(j)-1)
        lastcelldata = genefile %>% filter(cell == lastcell)
        alldata = lastcelldata %>% add_row(choosedata) %>% mutate(gap = nchar(cell)) %>%
          mutate(scale_time = gap+(time-min(time))/(max(time)-min(time)+1)) %>%
          mutate(newblot = blot -lag(blot)) %>% na.omit()
        allfun[[i]][[j]] = splinefun(x = alldata$scale_time,y = alldata$newblot)
      }
      else{
        choosedata = choosedata %>% mutate(gap = nchar(cell)) %>%
          mutate(scale_time = gap+(time-min(time))/(max(time)-min(time)+1))  %>%
          mutate(newblot = blot-lag(blot)) %>% na.omit()
        allfun[[i]][[j]] = splinefun(x = choosedata$scale_time,y = choosedata$newblot)
      }
    }
  }
  
  save(allfun,file = paste(savepath,subtree,"_allfun.RData",sep=""))
  # load(file = paste(realpath,subtree,"_allfun.RData",sep=""))
  newdata = matrix(0,nrow = length(cellname)*(M-1),
                   ncol = 2)
  newdata = data.frame(newdata)
  names(newdata) = c("cell","time")
  newdata$cell = rep(cellname,each=M-1)
  newdata = newdata %>% mutate(gap = nchar(cell)) %>% group_by(cell) %>% 
    mutate(time = seq(min(gap),min(gap)+1,length.out=M)[-M]) %>% ungroup() %>% select(-gap) 
  
  infodata = newdata
  
  singlecellinter_fun = function(thiscellname,infodata,thiscopyname,allfun){
    thisfun = allfun[[thiscopyname]][[thiscellname]]
    thiscelldata = infodata %>% filter(cell == thiscellname)
    if(is.null(thisfun)){
      tempvalue = rep(NA,nrow(thiscelldata))
    }
    else{
      tempvalue = thisfun(thiscelldata[["time"]])
    }
    return(tempvalue)
  }
  
  singlecopyinter_fun = function(thiscopyname,cellname,infodata,allfun){
    thiscopyvalue = unlist(sapply(cellname,singlecellinter_fun,infodata,thiscopyname,allfun,simplify = F))
    return(thiscopyvalue)
  }
  
  adddata = sapply(copyname,singlecopyinter_fun,cellname,infodata,allfun)
  newdata = newdata %>% bind_cols(adddata)
  names(newdata) = c("cell","time",copyname)
  
  allcellname = newdata$cell
  alltreenum = rep(1,nrow(newdata))
  nowmother = allcellname[1]
  for(i in 2:nrow(newdata)){
    if(str_sub(allcellname[i],1,nchar(nowmother))==nowmother){
      alltreenum[i] = alltreenum[i-1]
    }
    else{
      nowmother = allcellname[i]
      alltreenum[i] = alltreenum[i-1]+1
    }
  }
  
  tree_h = hash(allcellname,alltreenum)
  
  newdata$treeindex = 0
  for(i in 1:nrow(newdata)){
    newdata$treeindex[i] = tree_h[[newdata$cell[i]]]
  }
  
  newdata = newdata %>% arrange(cell,time) %>%  
    mutate(originrow = 1:n(),lastrow = originrow -1) %>% 
    group_by(treeindex) %>% mutate(ifchange = as.numeric(lag(cell)!=cell))  %>%
    select(cell,time,originrow,lastrow,treeindex,ifchange,everything()) %>% ungroup() %>% 
    mutate(ifchange = replace_na(ifchange,9))
  
  changeindex = which(newdata$ifchange==1)
  c_cell = newdata$cell
  lastrow = newdata$lastrow
  for(i in changeindex){
    thiscell = c_cell[i]
    lastcell = substr(thiscell,1,nchar(thiscell)-1)
    choosedata = newdata %>% filter(cell==lastcell) %>% select(originrow)
    temprow = choosedata$originrow
    lastrow[i] = temprow[length(temprow)]
  }
  newdata$lastrow = lastrow
  
  write_csv(newdata,file=paste(savepath,subtree,"_",maxdelta_time,"_condata.csv",sep=""))
  
}


# (5) hashData ------------------------------------------------------------
hashData = function(allfun,realdata,copyname,cellname,maxdelta_time,Mc,
                        savepath,subtree){
  M = round(Mc)
  eachdelta = 1/(M-1)
  maxdelta = maxdelta_time*eachdelta
  
  allhash = list()
  copynum = length(copyname)
  length(allhash) = copynum
  names(allhash) = copyname
  
  for(i in copyname){
    for(j in cellname){
      allhash[[i]][[j]] = hash() 
    }
  }
  
  allseq = seq(0,maxdelta,by=eachdelta)
  
  for(i in copyname){
    thiscopydata = realdata[[i]]
    if(all(is.na(thiscopydata))){
      next
    }
    for(j in cellname){
      choosedata = realdata %>% filter(cell==j)
      if(any(!is.na(choosedata[[i]]))){
        oldtime = choosedata$time
        if(any(is.na(choosedata[[i]]))){
          nonaindex = which(!is.na(choosedata[[i]]))
          oldtime = oldtime[nonaindex]
        }
        ini_point = floor(min(oldtime))

        newx = oldtime[1]
        newtime = c(sort(newx),oldtime)
        
        if(any(newtime<ini_point)){
          chooindex = which(newtime<ini_point)
          newtime_1 = newtime[-chooindex]
          newtime_2 = newtime[chooindex]
          lastcell = str_sub(j,1,nchar(j)-1)
          
          fcg1 = allfun[[i]][[j]]
          allhash[[i]][[j]] = hash(keys = round(newtime_1,3),values =fcg1(newtime_1))
          
          fcg2 = allfun[[i]][[lastcell]]
          .set(allhash[[i]][[lastcell]],keys = round(newtime_2,3),values = fcg2(newtime_2))
        }
        else{
          fcg = allfun[[i]][[j]]
          allhash[[i]][[j]] = hash(keys = round(newtime,3),values = fcg(newtime))
        }
      }
    }
  }
  
  save(allhash,file = paste(savepath,subtree,"_",maxdelta_time,"_allhash.RData",sep=""))
}


# (6) computeM ------------------------------------------------------------
computeM = function(hashdata,subtree,maxdelta_time,savepath){
  
  firstlist_fun = function(subfirst){
    tempvalue = as.numeric(values(subfirst))
    return(tempvalue)
  }
  
  secondlist_fun = function(sublist){
    tempvalue = sapply(sublist,firstlist_fun)
    return(median(unlist(tempvalue)))
  }
  
  M = sapply(hashdata,secondlist_fun) 
  
  save(M,file=paste(savepath,subtree,"_M_",maxdelta_time,"_.RData",sep=""))
}


# (7) mergeCopyFun --------------------------------------------------------
mergeCopyFun = function(allhash,genecopy_h,M,realdata,GeneDataset,InfoDataset,subtree,deletecellrate,
                        deletegenerate,copyname,cellname,mergerate,savepath){
  for(i in 1:length(allhash)){
    M_i = M[i]
    if(is.na(M_i)){
      next
    }
    for(j in cellname){
      if(length(keys(allhash[[i]][[j]]))==0){
        next
      }
      thekey = keys(allhash[[i]][[j]])
      thevalue = as.numeric(values(allhash[[i]][[j]]))
      biovalue = ifelse(thevalue>M_i,1,0)
      allhash[[i]][[j]] = hash(keys= thekey,values=biovalue)
    }
  }
  
  save(allhash,file = paste(savepath,"mer_",subtree,"_",maxdelta_time,"_bioallhash.RData",sep=""))
  # load(file = paste(realpath,"mer_",subtree,"_",maxdelta_time,"_bioallhash.RData",sep=""))
  allcopygene = c()
  for(i in copyname){
    thiscopygene = genecopy_h[[i]]
    allcopygene = c(allcopygene,thiscopygene)
  }
  uniquegene = sort(unique(allcopygene))

  newallhash = list()

  for(i in uniquegene){
    colindex = which(allcopygene == i)
    copynum = length(colindex)
    nacopy = 0
    for(q in colindex){
      thiscopy = copyname[q]
      if(all(is.na(GeneDataset[[thiscopy]]))){
        nacopy = nacopy+1
        next
      }
    }
    if(nacopy==copynum){
      next
    }

    newallhash[[i]] = list()

    for(k in cellname){
      mergedata = tibble(time = 0)
      for(j in colindex){
        thiscopy = copyname[j]
        nacopy = 0
        thishash = allhash[[j]][[k]]
        thiskeys = keys(thishash)
        thisvalues = values(thishash)
        if(length(thiskeys)==0){
          next
        }
        tempdata = data.frame(time = as.numeric(thiskeys),thisvalues)
        names(tempdata) = c("time",copyname[j])
        mergedata = mergedata %>% full_join(tempdata,by="time")
      }
      if(nrow(mergedata)==1){next}
      mergedata = mergedata %>% filter(time>0) %>% arrange(time) %>%
        mutate(blot = rowMeans(select(.,contains(i)),na.rm = T)) %>%
        mutate(blot = ifelse(blot >=mergerate,1,0))
      newallhash[[i]][[k]] = hash(keys = mergedata$time,values = mergedata$blot)
    }
  }
  
  save(newallhash,file = paste(savepath,"mer_",subtree,"_",maxdelta_time,"_origenehash.RData",sep=""))
  # load(file = paste(savepath,"mer_",subtree,"_",maxdelta_time,"_origenehash.RData",sep=""))
  
  selectgene = names(newallhash)
  bioGeneDataset = matrix(0,nrow =nrow(GeneDataset),ncol=length(selectgene))
  bioGeneDataset = as_tibble(bioGeneDataset)
  names(bioGeneDataset) = selectgene
  eachgenevalue = function(eachgene,newallhash,InfoDataset,cellname){
    thisvale = rep(0,nrow(InfoDataset))
    for(i in cellname){
      chooseindex = which(InfoDataset$cell==i)
      choosetime = round(InfoDataset$time[chooseindex],3)

      if(is.null(newallhash[[eachgene]][[i]])){
        thisvale[chooseindex] = rep(NA,length(chooseindex))
      }
      else{
        thisvale[chooseindex] = as.numeric(values(newallhash[[eachgene]][[i]],keys=choosetime))
      }
    }
    return(thisvale)
  }

  bioGeneDataset[,] = sapply(selectgene,eachgenevalue,newallhash,InfoDataset,cellname)

  write_csv(bioGeneDataset,file=paste(savepath,"mer_",subtree,"_",maxdelta_time,"_taskdata.csv",sep=""))
  # bioGeneDataset = read_csv(file=paste(savepath,"mer_",subtree,"_",maxdelta_time,"_taskdata.csv",sep=""),show_col_types = F)
  # vis_expect(bioGeneDataset, ~.x == 1)
  # ggsave(file = paste("./picture/",subtree,"_conven.pdf",sep=""),
  #        width = 15,
  #        height = 6)
  
  newstatdata = InfoDataset[,1] %>% add_column(bioGeneDataset) %>%
    gather(key="genename",value = "num",-cell) %>%group_by(cell,genename) %>%
    summarise(count = sum(!is.na(num)),.groups = "drop") %>% spread(key = genename,value = count)
  write_csv(newstatdata,file=paste(savepath,"mer_",subtree,"_",maxdelta_time,"_statdata.csv",sep=""))
  # newstatdata = read_csv(file=paste(savepath,"mer_",subtree,"_",maxdelta_time,"_statdata.csv",sep=""),show_col_types = F)
  # vis_expect(newstatdata[,-1], ~.x == 0)
  # ggsave(file = paste("./picture/",subtree,"_statdata.pdf",sep=""),
  #        width = 15,
  #        height = 6)
  
  selectinfo = delectFun(newstatdata,deletecellrate,deletegenerate)
  selectcell = selectinfo$selectcell
  selectgene = selectinfo$selectgene

  selectindex = which(names(bioGeneDataset)%in%selectgene)
  newbioGeneDataset = bioGeneDataset[,selectindex]
  mergedata = InfoDataset %>% add_column(newbioGeneDataset) %>%
    filter(cell %in% selectcell)
  allcellname = mergedata$cell
  alltreenum = rep(1,nrow(mergedata))
  nowmother = allcellname[1]
  for(i in 2:nrow(mergedata)){
    if(str_sub(allcellname[i],1,nchar(nowmother))==nowmother){
      alltreenum[i] = alltreenum[i-1]
    }
    else{
      nowmother = allcellname[i]
      alltreenum[i] = alltreenum[i-1]+1
    }
  }

  tree_h = hash(allcellname,alltreenum)

  mergedata$treeindex = 0
  for(i in 1:nrow(mergedata)){
    mergedata$treeindex[i] = tree_h[[mergedata$cell[i]]]
  }

  newmergedata = mergedata %>% arrange(cell,time) %>%
    mutate(originrow = 1:n(),lastrow = originrow -1) %>%
    group_by(treeindex) %>% mutate(ifchange = as.numeric(lag(cell)!=cell))  %>%
    select(cell,time,originrow,lastrow,treeindex,ifchange,everything()) %>% ungroup() %>%
    mutate(ifchange = replace_na(ifchange,9))

  changeindex = which(newmergedata$ifchange==1)
  c_cell = newmergedata$cell
  lastrow = newmergedata$lastrow
  for(i in changeindex){
    thiscell = c_cell[i]
    lastcell = substr(thiscell,1,nchar(thiscell)-1)
    choosedata = newmergedata %>% filter(cell==lastcell) %>% select(originrow)
    temprow = choosedata$originrow
    lastrow[i] = temprow[length(temprow)]
  }
  newmergedata$lastrow = lastrow

  write_csv(newmergedata,file=paste(savepath,"mer_",subtree,"_",maxdelta_time,"_selectdata.csv",sep=""))

  allstartrow = which(newmergedata[["ifchange"]]==9)
  allstartrow = as.numeric(sapply(allstartrow,function(x) seq(from=x,length.out=5)))
  allchooserow = setdiff(1:nrow(newmergedata),allstartrow)

  save(allstartrow,file=paste(savepath,"mer_",subtree,"_",maxdelta_time,"_allstartrow.RData",sep=""))
  save(allchooserow,file=paste(savepath,"mer_",subtree,"_",maxdelta_time,"_allchooserow.RData",sep=""))
}


# (8) delectFun -----------------------------------------------------------
delectFun = function(newstatdata,deletecellrate,deletegenerate){
  genename = names(newstatdata)[-1]
  cellname = newstatdata[[1]]
  tempdata = as.matrix(newstatdata[,-1])
  rowstat = apply(tempdata,1,function(x) sum(x==0)/length(x))
  derowindex = which(rowstat>deletecellrate)
  if(length(derowindex)>0){
    selectcell = cellname[-derowindex]
  }
  else{
    selectcell = cellname
  }
  newstatdata = newstatdata %>% filter(cell%in%selectcell)
  tempdata = as.matrix(newstatdata[,-1])
  colstat = apply(tempdata,2,function(x) sum(x==0)/length(x))
  decolindex = which(colstat>deletegenerate)
  if(length(decolindex)>0){
    selectgene = genename[-decolindex]
  }
  else{
    selectgene = genename
  }

  selectdata = newstatdata %>% select(all_of(selectgene))

  # vis_expect(selectdata, ~.x == 0)
  #
  # ggsave(file = paste("./picture/",subtree,"_cellgenesum.pdf",sep=""),
  #        width = 15,
  #        height = 6)

  return(list(selectcell = selectcell,
              selectgene = selectgene))
}