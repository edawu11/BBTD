# (1) singleAijPrior --------------------------------------------------------
singleAijPrior = function(genenum,lambda,ka){
  m = genenum-1
  z = 1:m
  cc = sum(exp(-lambda*z))
  row_logprior = function(x,m,lambda,ka,cc){
    if(x == 0 ){
      log_p = log(1-ka)
    }
    else{
      log_p = log(ka) +(-lambda*x)-log(cc)-log(choose(m,x))-x*log(2)
    }
    return(log_p)
  }
  allprior = sapply(0:m,row_logprior,m,lambda,ka,cc)
  priorhash = hash(0:m,allprior)
  return(priorhash)
}


# (2) invAijPrior ---------------------------------------------------------
invAijPrior = function(aijmat,alllogprior,genenum){
  rowmatrix_fun = function(rowmatrix,alllogprior,genenum){
    nozeronum = as.character(sum(abs(rowmatrix)))
    return(alllogprior[[nozeronum]])
  }
  temp = apply(aijmat,1,rowmatrix_fun,alllogprior,genenum)
  return(sum(temp))
}

# (3) invDeltaPrior -------------------------------------------------------
invDeltaPrior = function(updatedelta,updateaijmat,genenum){
  allaij = as.numeric(updateaijmat)
  alldelta = as.numeric(updatedelta)
  
  nozeroindex = which(allaij!=0)
  nozeroprob = length(nozeroindex)*log(1/6)
  
  allzeroindex = which(allaij==0)
  selfindex = seq(1,genenum^2,by=genenum+1)
  yeszeroindex = setdiff(allzeroindex,selfindex)
  zeronum = length(which(alldelta[yeszeroindex]==0.0))
  if(zeronum == length(yeszeroindex)){
    zeroprob = zeronum*log(0.999)
  }
  else{
    temp1 = zeronum*log(0.999)
    temp2 = (length(yeszeroindex)-zeronum)*log(0.0002)
    zeroprob = temp1 + temp2
  }
  return(nozeroprob+zeroprob)
}
# (4) geneAijRow -------------------------------------------------------
geneAijRow = function(aijindex,lambda,ka,genenum){
  bengenrow = rep(0,genenum)
  ifzero = rbinom(1,1,prob = ka)
  if(ifzero==0){
    return(bengenrow)
  }
  else{
    m = length(aijindex)
    z = 1:m
    sumexp = sum(exp(-lambda*z))
    p = exp(-lambda*z)/sumexp
    choosenum = sample(z,1,replace = F,prob = p)
    choosepos = sample(aijindex,choosenum,replace = F)
    bengenrow[choosepos] = sample(c(-1,1),choosenum,replace = T,prob = c(0.5,0.5))
    return(bengenrow)
  }
}


# (5) geneAij ----------------------------------------------------------
geneAij = function(genenum,lambda,ka,aijindexlist){
  bengenlist = lapply(aijindexlist,geneAijRow,lambda,ka,genenum)
  bengen = unlist(bengenlist)
  bengen = matrix(bengen,genenum,genenum,byrow = T)
  bengen = cbind(bengen,1:genenum)
  return(bengen)
}

# (6) subDeltaPrior -------------------------------------------------------
subDeltaPrior = function(subaij,maxdelta,eachdelta){
  if(subaij==0){
    delta = sample(seq(0,maxdelta,by=eachdelta),1,replace = F,prob = c(0.999,rep(0.0002,5)))
  }
  else{
    delta = sample(seq(0,maxdelta,by=eachdelta),1,replace = F,prob = rep(1/6,6))
  }
  return(delta)
}

# (7) geneDelta -----------------------------------------------------------
geneDelta = function(updateaijmat,maxdelta,eachdelta,genenum){
  aijnum = as.numeric(updateaijmat)
  
  alldelta = sapply(aijnum,subDeltaPrior,maxdelta,eachdelta)
  
  return(matrix(alldelta,genenum,genenum))
}
