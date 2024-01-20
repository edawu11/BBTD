# (1) mainMCMC ------------------------------------------------------------
mainMCMC = function(iters,genenum,genename,lastvaluelist,newallhash,
                           
                           GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,newaijindex,
                           
                           updateaijmat,updatebeta,updatealpha,updatelambda,updateka,updatedelta,
                           
                           betasd,alphasd,lambdasd,kasd,maxdelta,eachdelta){
  
  allpara = sparseMatrix(i=1,j=1,x=1,dims=c(iters,genenum^2*2+4))
  
  yesalpha = log(1/(1+exp(-updatealpha)))
  noalpha = log(1/(1+exp(updatealpha)))
  acc_beta = 0
  acc_alpha = 0
  acc_ka = 0
  acc_lambda = 0
  betaaccept = 0
  alphaaccept = 0
  lambdaaccept = 0
  
  oldcell = InfoDataset$cell
  oldtime = InfoDataset$time
  
  allprob = rep(0,iters)
  alllogprior = singleAijPrior(genenum,updatelambda,updateka)
  
  for(t_iter in 1:iters){
    

    updatetwo = updateTwoMat(updateaijmat,updatebeta,updatedelta,maxdelta,eachdelta,
                              
                              alllogprior,GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,
                              
                              genenum,genename,lastvaluelist,newallhash,oldcell,oldtime,
                              
                              yesalpha,noalpha,newaijindex)

    
    updateaijmat = updatetwo[[1]]
    updatedelta = updatetwo[[2]]
    
    betatemp = updateBeta(updatebeta,betasd,updateaijmat,updatedelta,
                           
                           GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,genename,genenum,lastvaluelist,
                           
                           oldcell,oldtime,newallhash,eachdelta,
                           
                           yesalpha,noalpha)
    
    updatebeta = betatemp[[1]]
    betaaccept = betatemp[[2]]
    acc_beta = acc_beta+betaaccept
    
    alphatemp = updateAlpha(updatealpha,alphasd,updateaijmat,updatedelta,updatebeta,
                             
                             GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,genename,genenum,lastvaluelist,
                             
                             oldcell,oldtime,newallhash,eachdelta)
    
    updatealpha = alphatemp[[1]]
    alphaaccept = alphatemp[[2]]
    acc_alpha = acc_alpha+alphaaccept
    yesalpha = log(1/(1+exp(-updatealpha)))
    noalpha = log(1/(1+exp(updatealpha)))

    lambdatemp = updateLambda(updatelambda,lambdasd,updateaijmat,genenum,alllogprior)
    updatelambda = lambdatemp[[1]]
    lambdaaccept = lambdatemp[[2]]
    acc_lambda = acc_lambda+lambdaaccept
    
    katemp = updateKa(updateka,kasd,updateaijmat,genenum)
    updateka = katemp[[1]]
    kaaccept = katemp[[2]]
    acc_ka = acc_ka+kaaccept
    
    alllogprior = singleAijPrior(genenum,updatelambda,updateka)
    
    allprob[t_iter] = postProb(updateaijmat,updatedelta,updatebeta,updatealpha,
                               updatelambda,updateka,GoalDataset,InfoDataset,
                               statdiff,allchooserow,allstartrow,
                               genename,genenum,lastvaluelist,oldcell,oldtime,maxdelta,eachdelta,alllogprior,newallhash,
                               yesalpha,noalpha)
    
    updateaijvec = as.numeric(updateaijmat)
    updatedeltavec = as.numeric(updatedelta)
    allpara[t_iter,1:(genenum^2)] = updateaijvec
    allpara[t_iter,(genenum^2+1):(2*genenum^2)] = updatedeltavec
    allpara[t_iter,2*genenum^2+1] = updatebeta
    allpara[t_iter,2*genenum^2+2] = updatealpha
    allpara[t_iter,2*genenum^2+3] = updatelambda
    allpara[t_iter,2*genenum^2+4] = updateka
  }
  result = list()
  allaccept = data.frame(acc_beta = acc_beta/iters,
                         acc_alpha = acc_alpha/iters,
                         acc_lambda = acc_lambda/iters,
                         acc_ka = acc_ka/iters)
  result[[1]] = allpara
  result[[2]] = allaccept
  result[[3]] = allprob
  return(result)
}

# (2) transExp ------------------------------------------------------------
transExp = function(x){
  if(any(x<(-700))){
    x_max = max(x)
    x_diff = (-50) - x_max
    y = x + x_diff
    a = exp(y)
    return(a)
  }
  else{
    a = exp(x)
    return(a)
  }
}

# (3) noZeroProb ---------------------------------------------------------
noZeroProb = function(lastvaluelist,singlegene,deltaindex,useaijvec,nozeroindex,
                         yesalpha,noalpha,updatebeta){
  valuemat = lastvaluelist[[singlegene]]
  thisnow = valuemat[,1]
  thislast = valuemat[,2]
  
  revaluemat = lastvaluelist[[nozeroindex[1]]]
  relast = revaluemat[,deltaindex[1],drop=F]
  
  if(length(nozeroindex)!=1){
    for(i in 2:length(nozeroindex)){
      revaluemat = lastvaluelist[[nozeroindex[i]]]
      relast = cbind(relast,revaluemat[,deltaindex[i],drop=F])
    }
  }

  H = as.numeric(useaijvec[nozeroindex]%*%t(relast))
  
  zeroindex = which(H==0)
  zerostat = abs(thisnow[zeroindex] - thislast[zeroindex])
  zerostat = table(factor(zerostat,levels=c(0,1)))
  lnzero =  zerostat[1]*yesalpha + zerostat[2]*noalpha

  oneindex = which(H!=0)
  Hbeta = updatebeta*H[oneindex]
  lnone = sum(2*thisnow[oneindex]*Hbeta - log(exp(2*Hbeta)+1))
  return(lnzero+lnone)
}


# (4) singleK ------------------------------------------------------------
singleK = function(k,oldtime,oldcell,nozeroindex,genename,choosedelta,
                       
                       newallhash,useaijvec,GoalDataset,singlegene,lastrow,
                       
                       yesalpha,noalpha,updatebeta){
  thistime = oldtime[k]
  thiscell = oldcell[k]
  ini_point = floor(thistime)
  lastvalue = rep(0,length(nozeroindex))
  for(j in 1:length(nozeroindex)){
    geneindex = nozeroindex[j]
    origene = genename[geneindex]
    lasttime = round(thistime - choosedelta[j],3)
    if(lasttime<ini_point){
      lastcell = str_sub(thiscell,1,nchar(thiscell)-1)
      thish = newallhash[[origene]][[lastcell]]
    }
    else{
      thish = newallhash[[origene]][[thiscell]]
    }

    if(is.null(thish)){
      lastvalue[j] = 0
    }
    else if(length(intersect(keys(thish),as.character(lasttime)))==0){
      lastvalue[j] = 0
    }
    else{
      lastvalue[j] = values(thish,lasttime)
    }
  }

  H_m=as.numeric(useaijvec[nozeroindex]%*%lastvalue)
  
  if(H_m==0){
    thislastrow = lastrow[k]
    lny = ifelse((GoalDataset[[singlegene]][k]==GoalDataset[[singlegene]][thislastrow]),yesalpha,noalpha)
  }
  else{
    lny = 2*updatebeta*GoalDataset[[singlegene]][k]*H_m-log(exp(2*updatebeta*H_m)+1)
  }
  return(lny)
}



# (5) singleLikelihood ---------------------------------------------------
singleLikelihood = function(singlegene,updateaijmat,updatedelta,updatebeta,
                            
                            GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,
                            
                            genename,genenum,lastvaluelist,
                            
                            oldcell,oldtime,eachdelta,
                            
                            newallhash,rowAijsum,
                            
                            yesalpha,noalpha){
  
  if(length(rowAijsum)==1){
    singlerowAijsum = rowAijsum
  }
  else{
    singlerowAijsum = rowAijsum[singlegene]
  }
  
  useaijvec = updateaijmat[singlegene,]
  if(singlerowAijsum==0){
    temp = statdiff[,singlegene]
    lny =  temp[1]*yesalpha + temp[2]*noalpha
  }
  
  else{
    nozeroindex = which(useaijvec!=0)
    choosedelta = updatedelta[singlegene,nozeroindex]
    deltaindex = round(choosedelta/eachdelta)+1
    
    lny = noZeroProb(lastvaluelist,singlegene,deltaindex,useaijvec,nozeroindex,
                        yesalpha,noalpha,updatebeta)
  }
  return(lny)
}




# (6) singleDelta --------------------------------------------------------
singleDelta = function(singleDelta,i_iter,j_iter,updateaijmat,updatedelta,updatebeta,
                          
                          GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,
                          
                          genename,genenum,lastvaluelist,oldcell,oldtime,eachdelta,
                          
                          newallhash,singlerowAijsum,yesalpha,noalpha){
  
  updatedelta[i_iter,j_iter] = singleDelta
  
  thisprob = singleLikelihood(i_iter,updateaijmat,updatedelta,updatebeta,
                              
                              GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,
                              
                              genename,genenum,lastvaluelist,
                              
                              oldcell,oldtime,eachdelta,
                              
                              newallhash,singlerowAijsum,
                              
                              yesalpha,noalpha)
  
  return(thisprob)
}


# (7) newTwoProb ---------------------------------------------------------
newTwoProb = function(newaij,alldelta,updateaijmat,updatedelta,updatebeta,
                          
                          i_iter,j_iter,
                          
                          GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,
                          
                          genename,genenum,lastvaluelist,
                          
                          oldcell,oldtime,eachdelta,
                          
                          alllogprior,newallhash,logzerodelta,logdelta,
                          
                          yesalpha,noalpha){

  if(newaij == 0){
    updateaijmat[i_iter,j_iter] = 0
    singlerowAijsum = sum(abs(updateaijmat[i_iter,]))
    allprob = singleLikelihood(i_iter,updateaijmat,updatedelta,updatebeta,
                               
                               GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,
                               
                               genename,genenum,lastvaluelist,
                               
                               oldcell,oldtime,eachdelta,
                               
                               newallhash,singlerowAijsum,
                               
                               yesalpha,noalpha)
    return(allprob +invAijPrior(updateaijmat,alllogprior,genenum)+ logzerodelta)
  }
  else{
    updateaijmat[i_iter,j_iter] = newaij
    singlerowAijsum = 1
    allprob = sapply(alldelta,singleDelta,i_iter,j_iter,updateaijmat,updatedelta,updatebeta,
                     
                     GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,
                     
                     genename,genenum,lastvaluelist,oldcell,oldtime,eachdelta,
                     
                     newallhash,singlerowAijsum,yesalpha,noalpha)
    return(allprob+invAijPrior(updateaijmat,alllogprior,genenum)+ logdelta)
  }
}


# (8) rowAij ------------------------------------------------------------------
rowAij = function(i_iter,newaijindex,allaij,alldelta,updateaijmat,updatedelta,updatebeta,
                     GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,
                     genename,genenum,lastvaluelist,
                     oldcell,oldtime,eachdelta,
                     alllogprior,newallhash,logzerodelta,logdelta,allnewtwo,
                     yesalpha,noalpha){
  for(j_iter in newaijindex[[i_iter]]){
    newaij_prob = sapply(allaij,newTwoProb,alldelta,updateaijmat,updatedelta,updatebeta,
                         
                         i_iter,j_iter,
                         
                         GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,
                         
                         genename,genenum,lastvaluelist,
                         
                         oldcell,oldtime,eachdelta,
                         
                         alllogprior,newallhash,logzerodelta,logdelta,
                         
                         yesalpha,noalpha)
    
    newaij_prob = as.numeric(newaij_prob)
    newaij_prob = transExp(newaij_prob)
    choosenum = sample(1:nrow(allnewtwo),1,F,newaij_prob)
    
    updateaijmat[i_iter,j_iter] = allnewtwo[choosenum,1]
    
    updatedelta[i_iter,j_iter] = allnewtwo[choosenum,2]
  }
  templist = list(rowAij = updateaijmat[i_iter,],
                  rowdelta = updatedelta[i_iter,])
  return(templist)
}


# (9) updateTwoMat -------------------------------------------------------
updateTwoMat =function(updateaijmat,updatebeta,updatedelta,maxdelta,eachdelta,
                        
                        alllogprior,GoalDataset,InfoDataset,statdiff,
                        
                        allchooserow,allstartrow,genenum,genename,lastvaluelist,
                        
                        newallhash,oldcell,oldtime,
                        
                        yesalpha,noalpha,newaijindex){
  
  alltwonum = round(maxdelta/eachdelta)+1
  
  allaij = c(0,1,-1)
  alldelta = seq(0,maxdelta,by=eachdelta)
  
  allnewtwo = matrix(c(rep(allaij,each=6),rep(alldelta,3)),ncol=2)
  logzerodelta = c(log(0.999),rep(log(0.0002),5))
  logdelta = c(rep(log(1/6),6))
  templist = lapply(1:genenum,rowAij,newaijindex,allaij,alldelta,updateaijmat,updatedelta,updatebeta,
         GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,
         genename,genenum,lastvaluelist,
         oldcell,oldtime,eachdelta,
         alllogprior,newallhash,logzerodelta,logdelta,allnewtwo,
         yesalpha,noalpha)
  
  updateaijmat = t(sapply(1:genenum,function(x) templist[[x]][[1]]))
  updatedelta = t(sapply(1:genenum,function(x) templist[[x]][[2]]))
  updatetwo = list()
  updatetwo[[1]] = updateaijmat
  updatetwo[[2]] = updatedelta
  return(updatetwo)
}


# (10) piBeta -------------------------------------------------------------
piBeta = function(beta,updateaijmat,updatedelta,
                     
                     GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,genename,
                     
                     oldcell,oldtime,newallhash,eachdelta,
                     
                     genenum,lastvaluelist,yesalpha,noalpha){
  
  rowAijsum = apply(abs(updateaijmat),1,sum)
  
  allprob = sum(sapply(1:genenum,singleLikelihood,updateaijmat,updatedelta,beta,
                       
                       GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,genename,genenum,lastvaluelist,
                       
                       oldcell,oldtime,eachdelta,
                       
                       newallhash,rowAijsum,
                       
                       yesalpha,noalpha))
  
  I = dgamma(beta,shape = 100,rate = 100,log = T)
  
  return(allprob+I)
}

# (11) updateBeta ---------------------------------------------------------
updateBeta = function(updatebeta,betasd,updateaijmat,updatedelta,
                       
                       GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,genename,genenum,lastvaluelist,
                       
                       oldcell,oldtime,newallhash,eachdelta,
                       
                       yesalpha,noalpha){
  
  newbeta = rtruncnorm(1,a=0,b=Inf,mean = updatebeta,sd=betasd)
  
  ratio = piBeta(newbeta,updateaijmat,updatedelta,
                    
                    GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,genename,
                    
                    oldcell,oldtime,newallhash,eachdelta,
                    
                    genenum,lastvaluelist, yesalpha,noalpha)-
    piBeta(updatebeta,updateaijmat,updatedelta,
              
              GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,genename,
              
              oldcell,oldtime,newallhash,eachdelta,
              
              genenum,lastvaluelist,yesalpha,noalpha)+
    
    log(dtruncnorm(updatebeta,a=0,b=Inf,mean = newbeta,sd=betasd))-
    
    log(dtruncnorm(newbeta,a=0,b=Inf,mean = updatebeta,sd=betasd))
  
  accept = c(0,ratio)
  result = list()
  if(log(runif(1)) <= min(accept)){
    result[[1]] = newbeta
    result[[2]] = 1
    return(result)
  }
  else{
    result[[1]] = updatebeta
    result[[2]] = 0
    return(result)
  }
}


# (12) piAlpha ------------------------------------------------------------
piAlpha = function(alpha,updateaijmat,updatedelta,updatebeta,
                      
                      GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,genename,
                      
                      oldcell,oldtime,newallhash,eachdelta,
                      
                      genenum,lastvaluelist){
  
  yesalpha = log(1/(1+exp(-alpha)))
  noalpha = log(1/(1+exp(alpha)))
  rowAijsum = apply(abs(updateaijmat),1,sum)
  
  allprob = sum(sapply(1:genenum,singleLikelihood,updateaijmat,updatedelta,updatebeta,
                       
                       GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,genename,genenum,lastvaluelist,
                       
                       oldcell,oldtime,eachdelta,
                       
                       newallhash,rowAijsum,
                       
                       yesalpha,noalpha))
  
  I = dgamma(alpha,shape = 1,rate = 10,log = T)
  
  return(allprob+I)
}


# (13) updateAlpha --------------------------------------------------------
updateAlpha = function(updatealpha,alphasd,updateaijmat,updatedelta,updatebeta,
                        
                        GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,genename,genenum,lastvaluelist,
                        
                        oldcell,oldtime,newallhash,eachdelta){
  
  newalpha = rtruncnorm(1,a=0,b=Inf,mean = updatealpha,sd=alphasd)
  ratio = piAlpha(newalpha,updateaijmat,updatedelta,updatebeta,
                     
                     GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,genename,
                     
                     oldcell,oldtime,newallhash,eachdelta,
                     
                     genenum,lastvaluelist)-
    piAlpha(updatealpha,updateaijmat,updatedelta,updatebeta,
               
               GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,genename,
               
               oldcell,oldtime,newallhash,eachdelta,
               
               genenum,lastvaluelist)+
    log(dtruncnorm(updatealpha,a=0,b=Inf,mean = newalpha,sd=alphasd))-
    log(dtruncnorm(newalpha,a=0,b=Inf,mean = updatealpha,sd=alphasd))
  
  accept = c(0,ratio)
  result = list()
  if(log(runif(1)) <= min(accept)){
    result[[1]] = newalpha
    result[[2]] = 1
    return(result)
  }
  else{
    result[[1]] = updatealpha
    result[[2]] = 0
    return(result)
  }
}


# (14) piLambda -----------------------------------------------------------
piLambda = function(lambda,updateaijmat,genenum,alllogprior){
  allprob = invAijPrior(updateaijmat,alllogprior,genenum)
  I = dgamma(lambda,shape = 490,rate= 70)
  return(allprob+log(I))
}


# (15) updateLambda -------------------------------------------------------
updateLambda = function(updatelambda,lambdasd,updateaijmat,genenum,alllogprior){
  newlambda = rtruncnorm(1,a=0,b=Inf,mean = updatelambda,sd=lambdasd)
  ratio = piLambda(newlambda,updateaijmat,genenum,alllogprior)-
    piLambda(updatelambda,updateaijmat,genenum,alllogprior)+
    log(dtruncnorm(updatelambda,a=0,b=Inf,mean = newlambda,sd=lambdasd))-
    log(dtruncnorm(newlambda,a=0,b=Inf,mean = updatelambda,sd=lambdasd))
  accept = c(0,ratio)
  result = list()
  if(log(runif(1)) <= min(accept)){
    result[[1]] = newlambda
    result[[2]] = 1
    return(result)
  }
  else{
    result[[1]] = updatelambda
    result[[2]] = 0
    return(result)
  }
}


# (16) piKa ---------------------------------------------------------------
piKa = function(updateka,updateaijmat,genenum){
  aijrowsum = apply(abs(updateaijmat),1,sum)
  aijrowsum[aijrowsum>0] = 1
  onesum = sum(aijrowsum)
  zerosum = genenum - onesum
  allprob  = log((updateka^onesum)*((1-updateka)^zerosum))
  I = dbeta(updateka,shape1 = 24,shape2 = 6,log = T)
  return(allprob+I)
}


# (17) updateKa -----------------------------------------------------------
updateKa = function(updateka,kasd,updateaijmat,genenum){
  newka = rtruncnorm(1,a=0,b=1,mean = updateka,sd=kasd)
  ratio = piKa(newka,updateaijmat,genenum)-
    piKa(updateka,updateaijmat,genenum)+
    log(dtruncnorm(updateka,a=0,b=1,mean = newka,sd=kasd))-
    log(dtruncnorm(newka,a=0,b=1,mean = updateka,sd=kasd))
  accept = c(0,ratio)
  result = list()
  if(log(runif(1)) <= min(accept)){
    result[[1]] = newka
    result[[2]] = 1
    return(result)
  }
  else{
    result[[1]] = updateka
    result[[2]] = 0
    return(result)
  }
}

# (18) postProb -----------------------------------------------------------
postProb = function(updateaijmat,updatedelta,updatebeta,updatealpha,
                        updatelambda,updateka,GoalDataset,InfoDataset,
                        statdiff,allchooserow,allstartrow,
                        genename,genenum,lastvaluelist,oldcell,oldtime,
                        maxdelta,eachdelta,alllogprior,newallhash,
                        yesalpha,noalpha){
  
  allrowAijsum = apply(abs(updateaijmat),1,sum)
  
  likelihood_prob = sum(sapply(1:genenum,singleLikelihood,updateaijmat,updatedelta,updatebeta,
                               
                               GoalDataset,InfoDataset,statdiff,allchooserow,allstartrow,genename,genenum,lastvaluelist,
                               
                               oldcell,oldtime,eachdelta,
                               
                               newallhash,allrowAijsum,
                               
                               yesalpha,noalpha))
  
  aijmat_prob = invAijPrior(updateaijmat,alllogprior,genenum)

  delta_prob = InvDeltaPrior(updatedelta,updateaijmat,genenum)
  
  beta_prob = dgamma(updatebeta,shape =100,rate = 100,log = T)
  
  alpha_prob = dgamma(updatealpha,shape =1,rate = 10,log = T)
  
  lambda_prob = dgamma(updatelambda,shape = 490,rate= 70,log = T)
  
  ka_prob = dbeta(updateka,shape1 = 24,shape2= 6,log = T)
  
  allprob = likelihood_prob+aijmat_prob+delta_prob+beta_prob+alpha_prob+lambda_prob+ka_prob
  return(allprob)
}

# (19) realMCMC -----------------------------------------------------------
realMCMC = function(subtree,p=1,mcmcnum=5,num_maxiters=20000,iters=500,
                    betasd=0.1,alphasd=0.1,lambdasd=0.1,kasd=0.01,maxdelta_time,realpath,realresultpath){
  load(file=paste(realpath,subtree,"_simu_info.RData",sep=""))
  Mc = simu_info$Mc
  eachdelta = 1/(round(Mc)-1)
  maxdelta = maxdelta_time*eachdelta
  load(paste(realpath,"mer_",subtree,"_",maxdelta_time,"_allchooserow.RData",sep=""))
  load(paste(realpath,"mer_",subtree,"_",maxdelta_time,"_allstartrow.RData",sep=""))
  load(paste(realpath,"mer_",subtree,"_",maxdelta_time,"_origenehash.RData",sep=""))
  newdata = read_csv(file=paste(realpath,"mer_",subtree,"_",maxdelta_time,"_selectdata.csv",sep=""),show_col_types = FALSE)
  finalinfo = newdata[,1:6]
  finaldata = newdata[,-c(1:6)]
  oldtime = finalinfo$time
  oldcell = finalinfo$cell
  selectgene = names(finaldata)
  genenum = length(selectgene)
  aijindexlist = lapply(1:genenum,function(x,genenum) list=setdiff(1:genenum,x),genenum)
  names(aijindexlist) = selectgene
  load(file=paste(realpath,"mer_",subtree,"_newaijindex_",p,"_.RData",sep=""))
  load(paste(realpath,"mer_",subtree,"_",p,"_diff",".RData",sep=""))
  load(paste(realpath,subtree,"_",p,"_lastvaluelist.RData",sep=""))
  finaldata = map_dfc(finaldata,function(x) replace_na(x,0))
  alliters = 0
  firstsinglemcmc = function(i,p,genenum,lastvaluelist,iters,selectgene,newallhash,aijindexlist,
                             
                             finaldata,finalinfo,statdiff,allchooserow,allstartrow,newaijindex,
                             
                             betasd,alphasd,lambdasd,kasd,maxdelta,eachdelta){
    set.seed(p+i)
    updatebeta = rgamma(1,shape = 100,rate =100)
    updatealpha = rgamma(1,shape = 1,rate  = 10)
    updatelambda = rgamma(1,shape =490,rate= 70)
    updateka = rbeta(1,shape1 =24,shape2= 6)
    
    aijmat = geneAij(genenum,updatelambda,updateka,aijindexlist)
    updateaijmat = aijmat[,-(genenum+1)]
    for(j in 1:nrow(updateaijmat)){
      zeroindex = setdiff(aijindexlist[[j]],newaijindex[[j]])
      updateaijmat[j,zeroindex] = 0
    }
    updatedelta = geneDelta(updateaijmat,maxdelta,eachdelta,genenum)
    
    chainresult = mainMCMC(iters,genenum,selectgene,lastvaluelist,newallhash,
                           
                           finaldata,finalinfo,statdiff,allchooserow,allstartrow,newaijindex,
                           
                           updateaijmat,updatebeta,updatealpha,updatelambda,updateka,updatedelta,
                           
                           betasd,alphasd,lambdasd,kasd,maxdelta,eachdelta)
    
    chain = chainresult[[1]]
    prob = chainresult[[3]]
    updateaijmat = matrix(chain[iters,1:(genenum^2)],genenum,genenum)
    updatedelta = matrix(chain[iters,(genenum^2+1):(2*genenum^2)],genenum,genenum)
    updatebeta = chain[iters,2*genenum^2+1]
    updatealpha = chain[iters,2*genenum^2+2]
    updatelambda = chain[iters,2*genenum^2+3]
    updateka = chain[iters,2*genenum^2+4]
    
    tmpresult = list(chain = chain,
                     prob = prob,
                     updateaijmat = updateaijmat,
                     updatedelta = updatedelta,
                     updatebeta = updatebeta,
                     updatealpha = updatealpha,
                     updatelambda = updatelambda,
                     updateka = updateka)
    
    return(tmpresult)
  }
  
  allchain = list()
  allprob = list()
  allpara = list()
  firstlist = sfLapply(1:mcmcnum,firstsinglemcmc,p,genenum,lastvaluelist,iters,selectgene,newallhash,aijindexlist,
                       
                       finaldata,finalinfo,statdiff,allchooserow,allstartrow,newaijindex,
                       
                       betasd,alphasd,lambdasd,kasd,maxdelta,eachdelta)
  for(i in 1:mcmcnum){
    thislist = firstlist[[i]]
    chain = thislist[[1]]
    thisprob = thislist[[2]]
    allpara[[i]] = list(updateaijmat = thislist[[3]],
                        updatedelta = thislist[[4]],
                        updatebeta = thislist[[5]],
                        updatealpha = thislist[[6]],
                        updatelambda = thislist[[7]],
                        updateka = thislist[[8]])
    allchain[[i]] = sparseMatrix(i=1,j=1,x=1,dims=c(iters,2*genenum^2+4))
    allchain[[i]][,1:(genenum^2)] = chain[,1:(genenum^2)]
    allchain[[i]][,(genenum^2+1):(2*genenum^2)] = chain[,(genenum^2+1):(2*genenum^2)]
    allchain[[i]][,2*genenum^2+1] = chain[,2*genenum^2+1]
    allchain[[i]][,2*genenum^2+2] = chain[,2*genenum^2+2]
    allchain[[i]][,2*genenum^2+3] = chain[,2*genenum^2+3]
    allchain[[i]][,2*genenum^2+4] = chain[,2*genenum^2+4]
    allprob[[i]] = thisprob
  }
  
  save(allchain,file = paste(realresultpath,subtree,"_",p,"_result",".RData",sep=""))
  save(allprob,file= paste(realresultpath,subtree,"_",p,"_allprob",".RData",sep=""))
  choosepara = maxEstimate(allchain,allprob,genenum)
  preaij = choosepara[[1]]
  predelta = choosepara[[2]]
  preparas = choosepara[[3]]
  write.csv(preaij,file=paste(realresultpath,subtree,"_",p,"_preaij.csv",sep=""),row.names = F)
  write.csv(predelta,file=paste(realresultpath,subtree,"_",p,"_predelta.csv",sep=""),row.names = F)
  write.csv(preparas,file=paste(realresultpath,subtree,"_",p,"_preparas.csv",sep=""),row.names = F)
  alliters = alliters+ iters
  subchain = list()
  eachsinglemcmc = function(i,p,allpara,genenum,lastvaluelist,iters,selectgene,newallhash,
                            
                            finaldata,finalinfo,statdiff,allchooserow,allstartrow,newaijindex,
                            
                            betasd,alphasd,lambdasd,kasd,maxdelta,eachdelta){
    set.seed(p+i)
    updatebeta = allpara[[i]][["updatebeta"]]
    updatealpha = allpara[[i]][["updatealpha"]]
    updatelambda = allpara[[i]][["updatelambda"]]
    updateka = allpara[[i]][["updateka"]]
    updateaijmat = allpara[[i]][["updateaijmat"]]
    updatedelta = allpara[[i]][["updatedelta"]]
    chainresult = mainMCMC(iters,genenum,selectgene,lastvaluelist,newallhash,
                           
                           finaldata,finalinfo,statdiff,allchooserow,allstartrow,newaijindex,
                           
                           updateaijmat,updatebeta,updatealpha,updatelambda,updateka,updatedelta,
                           
                           betasd,alphasd,lambdasd,kasd,maxdelta,eachdelta)
    
    chain = chainresult[[1]]
    prob = chainresult[[3]]
    updateaijmat = matrix(chain[iters,1:(genenum^2)],genenum,genenum)
    updatedelta = matrix(chain[iters,(genenum^2+1):(2*genenum^2)],genenum,genenum)
    updatebeta = chain[iters,2*genenum^2+1]
    updatealpha = chain[iters,2*genenum^2+2]
    updatelambda = chain[iters,2*genenum^2+3]
    updateka = chain[iters,2*genenum^2+4]
    tmpresult = list(chain = chain,
                     prob = prob,
                     updateaijmat = updateaijmat,
                     updatedelta = updatedelta,
                     updatebeta = updatebeta,
                     updatealpha = updatealpha,
                     updatelambda = updatelambda,
                     updateka = updateka)
    
    return(tmpresult)
  }
  while(alliters<num_maxiters){
    alliters = alliters+ iters
    eachlist = sfLapply(1:mcmcnum,eachsinglemcmc,p,allpara,genenum,lastvaluelist,iters,selectgene,newallhash,
                        
                        finaldata,finalinfo,statdiff,allchooserow,allstartrow,newaijindex,
                        
                        betasd,alphasd,lambdasd,kasd,maxdelta,eachdelta)
    for(i in 1:mcmcnum){
      thislist = eachlist[[i]]
      chain = thislist[[1]]
      thisprob = thislist[[2]]
      allpara[[i]] = list(updateaijmat = thislist[[3]],
                          updatedelta = thislist[[4]],
                          updatebeta = thislist[[5]],
                          updatealpha = thislist[[6]],
                          updatelambda = thislist[[7]],
                          updateka = thislist[[8]])
      subchain[[i]] = sparseMatrix(i=1,j=1,x=1,dims=c(iters,2*genenum^2+4))
      subchain[[i]][,1:(genenum^2)] = chain[,1:(genenum^2)]
      subchain[[i]][,(genenum^2+1):(2*genenum^2)] = chain[,(genenum^2+1):(2*genenum^2)]
      subchain[[i]][,2*genenum^2+1] = chain[,2*genenum^2+1]
      subchain[[i]][,2*genenum^2+2] = chain[,2*genenum^2+2]
      subchain[[i]][,2*genenum^2+3] = chain[,2*genenum^2+3]
      subchain[[i]][,2*genenum^2+4] = chain[,2*genenum^2+4]
      allchain[[i]] = rbind(allchain[[i]],subchain[[i]])
      allprob[[i]] = c(allprob[[i]],thisprob)
    }
    save(allchain,file = paste(realresultpath,subtree,"_",p,"_result",".RData",sep=""))
    save(allprob,file= paste(realresultpath,subtree,"_",p,"_allprob",".RData",sep=""))
    choosepara = maxEstimate(allchain,allprob,genenum)
    preaij = choosepara[[1]]
    predelta = choosepara[[2]]
    preparas = choosepara[[3]]
    write.csv(preaij,file=paste(realresultpath,subtree,"_",p,"_preaij.csv",sep=""),row.names = F)
    write.csv(predelta,file=paste(realresultpath,subtree,"_",p,"_predelta.csv",sep=""),row.names = F)
    write.csv(preparas,file=paste(realresultpath,subtree,"_",p,"_preparas.csv",sep=""),row.names = F)
  }
}

