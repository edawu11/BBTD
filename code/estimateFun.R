# (1) catchTable ---------------------------------------------------------
catchTable = function(a,num){
  a = as.numeric(a)
  a = factor(a,levels = c(-1,0,1),labels = c(-1,0,1))
  temp = table(a)
  namenum = which(names(temp)==as.character(num))
  return(as.numeric(temp)[namenum])
}

# (2) singleMode ---------------------------------------------------------
singleMode = function(allchain,mcmcnum,genenum,eachdelta,dir_result){
  startindex = (nrow(allchain[[1]])/2+1)
  endindex = nrow(allchain[[1]])
  temp = allchain[[1]][startindex:endindex,1:(2*genenum^2)]
  if(mcmcnum>2){
    for(k in 2:mcmcnum){
      temp = rbind(temp,allchain[[k]][startindex:endindex,1:(2*genenum^2)])
    }
  }
  
  temp[,(genenum^2+1):(2*genenum^2)] = round(temp[,(genenum^2+1):(2*genenum^2)]/eachdelta,0)
  
  rowcbind_f = function(eachrow){
    rowvector = paste0(round(eachrow,2),collapse = "")
    return(rowvector)
  }
  
  rowstat = apply(temp,1,rowcbind_f)
  stattable = table(rowstat)
  maxindex = which.max(stattable)
  choosevalue = names(stattable[maxindex])
  
  choosevalue = as.numeric(str_split(choosevalue,pattern = "",simplify = F)[[1]]) 
  preaij = matrix(choosevalue[1:(genenum^2)],genenum,genenum)  
  predelta = matrix(choosevalue[(genenum^2+1):(2*genenum^2)]*eachdelta,genenum,genenum) 
  
  # return(list(preaij = preaij,predelta=predelta))
  
  write.csv(preaij,file=paste(dir_result,bosscell,"_",p,"_preaij.csv",sep=""),row.names = F)
  write.csv(predelta,file=paste(dir_result,bosscell,"_",p,"_predelta.csv",sep=""),row.names = F)
}

# (3) maxEstimate --------------------------------------------------------
maxEstimate = function(allchain,allprob,genenum){
  mcmcnum = length(allprob)
  eachlen = nrow(allchain[[1]])
  for(i in 1:mcmcnum){
    allchain[[i]] = allchain[[i]][(eachlen/2+1):eachlen,,drop=F]
    allprob[[i]] = allprob[[i]][(eachlen/2+1):eachlen]
  }
  
  listvalue = c()
  for(i in 1:mcmcnum){
    listvalue = c(listvalue,max(allprob[[i]]))
  }
  maxlist_index = which.max(listvalue)
  chooseprob = allprob[[maxlist_index]]
  maxprob_index = which.max(chooseprob)
  choosevalue = allchain[[maxlist_index]][maxprob_index,]
  
  chooseaijmat = matrix(choosevalue[1:(genenum^2)],genenum,genenum)
  choosedelta = matrix(choosevalue[((genenum^2)+1):(2*genenum^2)],genenum,genenum)
  choosebeta = choosevalue[2*genenum^2+1]
  choosealpha = choosevalue[2*genenum^2+2]
  chooselambda = choosevalue[2*genenum^2+3]
  chooseka = choosevalue[2*genenum^2+4]
  
  finalvalue = data.frame(beta = choosebeta,
                          alpha = choosealpha,
                          lambda = chooselambda,
                          ka = chooseka)
  result = list()
  result[[1]] = chooseaijmat
  result[[2]] = choosedelta
  result[[3]] = finalvalue
  return(result)
} 

# (4) evaAijDelta ----------------------------------------------------------------
evaAijDelta = function(a_gene,b_gene,a_delta,b_delta,genenum){
  pingjia=list()
  
  choosenum=0
  allnum_cujin = catchTable(a_gene,1)
  if(allnum_cujin==0){
    pingjia$TPR = 99
  }
  else{
    for(i in 1:genenum){
      for(j in 1:genenum){
        if((a_gene[i,j]==1)&(b_gene[i,j]==1)){
          choosenum = choosenum+1
        }
      }
    }
    pingjia$TPR = choosenum/allnum_cujin
  }
  
  choosenum=0
  allnum_cujin = catchTable(b_gene,1)
  if(allnum_cujin==0){
    pingjia$PPR = 99
  }
  else{
    for(i in 1:genenum){
      for(j in 1:genenum){
        if((a_gene[i,j]==1)&(b_gene[i,j]==1)){
          choosenum = choosenum+1
        }
      }
    }
    pingjia$PPR = choosenum/allnum_cujin
  }
  
  choosenum=0
  allnum_yizhi = catchTable(a_gene,-1)
  if(allnum_yizhi==0){
    pingjia$TNR = 99
  }
  else{
    for(i in 1:genenum){
      for(j in 1:genenum){
        if((a_gene[i,j]==-1)&(b_gene[i,j]==-1)){
          choosenum = choosenum+1
        }
      }
    }
    pingjia$TNR = choosenum/allnum_yizhi
  }
  
  choosenum=0
  allnum_yizhi = catchTable(b_gene,-1)
  if(allnum_yizhi==0){
    pingjia$PNR = 99
  }
  else{
    for(i in 1:genenum){
      for(j in 1:genenum){
        if((a_gene[i,j]==-1)&(b_gene[i,j]==-1)){
          choosenum = choosenum+1
        }
      }
    }
    pingjia$PNR = choosenum/allnum_yizhi
  }
  
  a_delta = round(a_delta,3)
  b_delta = round(b_delta,3)
  choosenum=0
  allnum_nozero = catchTable(b_gene,1)+catchTable(b_gene,-1)
  if(allnum_nozero==0){
    pingjia$ACC = 1
  }
  else{
    for(i in 1:genenum){
      for(j in setdiff(1:genenum,i)){
        if(a_delta[i,j]==b_delta[i,j] & b_gene[i,j]!=0 ){
          choosenum = choosenum+1
        }
      }
    }
    
    pingjia$ACC = choosenum/allnum_nozero
  }
  
  pingjia=as.data.frame(pingjia)
  return(pingjia)
}