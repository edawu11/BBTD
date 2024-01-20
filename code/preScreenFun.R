# (1) preScreen -----------------------------------------------------------
preScreen = function(aijindexlist,selectgene,GoalDataset,finalinfo,allchooserow,realaijmat,
                                       maxdelta_time,Mc,yuzhipoint,ifallpvalue,ifreal,simupath,p,subtree,thishash){
  
  eachdelta = 1/(round(Mc)-1)
  maxdelta = maxdelta_time*eachdelta
  genenum = length(selectgene)
  
  if(ifallpvalue==T){
    all_pvalue = aijPvalue(aijindexlist,selectgene,genenum,GoalDataset,finalinfo,allchooserow,
                                     maxdelta,eachdelta,thishash)
    save(all_pvalue,file=paste(simupath,"mer_",subtree,"_all_pvalue_",p,"_","_.RData",sep=""))
  }
  else{
    load(paste(simupath,"mer_",subtree,"_all_pvalue_",p,"_","_.RData",sep=""))
  }
  
  newaijindex = aijindexlist
  for(i in 1:length(all_pvalue)){
    thepvalue = all_pvalue[[i]]
    tempindex = which(thepvalue < yuzhipoint)
    newaijindex[[i]] = aijindexlist[[i]][tempindex]
  }
  
  if(ifreal == F){
    
    initialresult = list()
    initialresult[[1]] = length(unlist(newaijindex))
    initialmat = matrix(0,nrow = genenum,ncol =genenum)
    for(i in 1:length(newaijindex)){
      initialmat[i,newaijindex[[i]]]=1
    }
    colnowgene = factor(as.numeric(initialmat),levels = c(0,1),labels = c(0,1))
    conmat = table(Prediction = factor(as.numeric(initialmat),levels = c(0,1),labels = c(0,1)),
                   Reference = factor(as.numeric(as.matrix(abs(realaijmat[,1:genenum]))),levels = c(0,1),labels = c(0,1)))
    conmat[1,1] = conmat[1,1] - genenum
    
    if(sum(conmat[,2])==0){
      FPR=0
    }
    else{
      FPR = conmat[1,2]/sum(conmat[,2])
    }
    conf_matrix = confusionMatrix(conmat)
    initialresult[[2]] = conf_matrix
    initialresult[[3]] = genenum*(genenum-1)-initialresult[[1]]
    initialresult[[4]] = FPR
    initialresult[[5]] = conmat[,2]
    
    save(initialresult,file=paste(simupath,"mer_",subtree,"_initial_result_",p,"_.RData",sep=""))
  }
  
  save(newaijindex,file=paste(simupath,"mer_",subtree,"_newaijindex_",p,"_.RData",sep=""))
}

# (2) aijPvalue -----------------------------------------------------------
aijPvalue = function(aijindexlist,selectgene,genenum,GoalDataset,finalinfo,allchooserow,
                     maxdelta,eachdelta,thishash){
  oldcell = finalinfo$cell
  oldtime = finalinfo$time
  alldelta = seq(0,maxdelta,by=eachdelta)
  genenum = length(selectgene)
  subPvalue = function(numindex,aijindexlist,selectgene,genenum,allchooserow,oldcell,oldtime,alldelta,thishash){
    thisgenename = selectgene[numindex]
    
    colfuturegene = GoalDataset[[numindex]][allchooserow]
    colfuturegene = replace_na(colfuturegene,0)
    
    eachMin = function(eachindex,alldelta,thisgenename,eachdelta,colfuturegene,
                       finalinfo,allchooserow,oldcell,oldtime,thishash){
      origenename = selectgene[eachindex]
      origenenowvalue = sapply(alldelta,orivalueFun,origenename,thisgenename,eachdelta,
                               finalinfo,allchooserow,oldcell,oldtime,thishash)
      thisgenepvalue = apply(origenenowvalue,2,inGeneCol,colfuturegene)
      return(min(thisgenepvalue))
    }
    
    thisgeneallpvalue = sapply(aijindexlist[[numindex]],eachMin,alldelta,thisgenename,eachdelta,colfuturegene,
                               finalinfo,allchooserow,oldcell,oldtime,thishash)
    
    return(thisgeneallpvalue)
  }
  
  allpvalue = sapply(1:genenum,subPvalue,aijindexlist,selectgene,genenum,
                     allchooserow,oldcell,oldtime,alldelta,thishash,simplify = F)
  names(allpvalue) = selectgene
  return(allpvalue)
}

# (3) orivalueFun -------------------------------------------------------------
orivalueFun = function(choosedelta,origenename,futgenename,eachdelta,
                        finalinfo,allchooserow,oldcell,oldtime,thishash){
  
  alllastvalue = rep(0,length(allchooserow))
  lastrow = finalinfo$lastrow
  
  singleGeneOrivalue = function(k,allchooserow,choosedelta,origenename,lastrow){
    oldrow = allchooserow[k]
    thistime = oldtime[oldrow]
    thiscell = oldcell[oldrow]
    ini_point = floor(thistime)
    lasttime = round(thistime - choosedelta,3)
    if(lasttime<ini_point){
      lastcell = str_sub(thiscell,1,-2)
      thish = thishash[[origenename]][[lastcell]]
      if(is.null(thish)){
        lastvalue = 0
      }
      else if(length(intersect(keys(thish),as.character(lasttime)))==0){
        lastvalue = 0
      }
      else{
        lastvalue = values(thish,round(lasttime,3))
      }
      
    }
    else{
      thish = thishash[[origenename]][[thiscell]]
      if(is.null(thish)){
        lastvalue = 0
      }
      else if(length(intersect(keys(thish),as.character(lasttime)))==0){
        lastvalue = 0
      }
      else{
        lastvalue = values(thish,round(lasttime,3))
      }
    }
    return(as.numeric(lastvalue))
  }
  
  alllastvalue = sapply(1:length(allchooserow),singleGeneOrivalue,allchooserow,
                        choosedelta,origenename,lastrow)
}


# (4) inGeneCol -----------------------------------------------------------
inGeneCol = function(colnowgene,colfuturegene){
  colnowgene = factor(colnowgene,levels = c(0,1),labels = c(0,1))
  colfuturegene = factor(colfuturegene,levels = c(0,1),labels = c(0,1))
  ktab = table(colnowgene,colfuturegene)
  colnowgene = as.numeric(as.character(colnowgene))
  colfuturegene = as.numeric(as.character(colfuturegene))
  
  
  if(any(ktab==0)){
    if(ktab[1]==0){
      colnowgene = c(colnowgene,0)
      colfuturegene = c(colfuturegene,0)
    }
    if(ktab[2]==0){
      colnowgene = c(colnowgene,1)
      colfuturegene = c(colfuturegene,0)
    }
    if(ktab[3]==0){
      colnowgene = c(colnowgene,0)
      colfuturegene = c(colfuturegene,1)
    }
    if(ktab[4]==0){
      colnowgene = c(colnowgene,1)
      colfuturegene = c(colfuturegene,1)
    }
    
    result = fisher.test(colnowgene,colfuturegene)
    return(result$p.value)
    
  }
  else{
    result = fisher.test(colnowgene,colfuturegene)
    return(result$p.value)
  }
}