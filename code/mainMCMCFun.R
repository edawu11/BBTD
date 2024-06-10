# (1) mainMCMC ------------------------------------------------------------
# Define the function mainMCMC to perform MCMC simulation for sampling parameters
# Inputs:
#   iters: Integer, number of iterations for the MCMC simulation
#   genenum: Integer, number of genes
#   genename: Vector, names of the genes
#   lastvaluelist: List, contains last observed values for each gene
#   newallhash: List, contains hashed tables of prior interactions
#   GoalDataset: DataFrame, target data for the MCMC optimization
#   InfoDataset: DataFrame, contains timing and cell information for each observation
#   statdiff: Numeric, a statistic differential used in parameter updates
#   allchooserow: Vector, indices of rows chosen for analysis
#   allstartrow: Vector, indices of the start rows for each segment of data
#   newaijindex: List, indices of potential interactions for each gene
#   updateaijmat: Matrix, current estimates of interaction strengths
#   updatebeta: Numeric, continuous parameter
#   updatedelta: Matrix, current estimates of delta parameters
#   betasd, alphasd, lambdasd, kasd: Numerics, standard deviations for parameter updates
#   maxdelta, eachdelta: Numerics, parameters to define the range and steps of delta for the simulations
# Output:
#   List containing the parameter estimates, acceptance rates, and probability calculations across all iterations
mainMCMC = function(iters, genenum, genename, lastvaluelist, newallhash,
                    GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow, newaijindex,
                    updateaijmat, updatebeta, updatealpha, updatelambda, updateka, updatedelta,
                    betasd, alphasd, lambdasd, kasd, maxdelta, eachdelta){
  
  # Initialize a sparse matrix to store all parameters and their updates for each iteration
  allpara = sparseMatrix(i = 1, j = 1, x = 1, dims = c(iters, genenum^2 * 2 + 4))
  
  # Compute log-odds for alpha
  yesalpha = log(1 / (1 + exp(-updatealpha)))
  noalpha = log(1 / (1 + exp(updatealpha)))
  # Initialize counters for parameter acceptance
  acc_beta = 0
  acc_alpha = 0
  acc_ka = 0
  acc_lambda = 0
  betaaccept = 0
  alphaaccept = 0
  lambdaaccept = 0
  
  # Extract cell and timing information from InfoDataset
  oldcell = InfoDataset$cell
  oldtime = InfoDataset$time
  
  # Prepare to collect probabilities and log-priors
  allprob = rep(0, iters)
  alllogprior = singleAijPrior(genenum, updatelambda, updateka)
  
  # MCMC iteration loop
  for(t_iter in 1:iters){
    # Update aij and delta matrices
    updatetwo = updateTwoMat(updateaijmat, updatebeta, updatedelta, maxdelta, eachdelta,
                             alllogprior, GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow,
                             genenum, genename, lastvaluelist, newallhash, oldcell, oldtime,
                             yesalpha, noalpha, newaijindex)

    updateaijmat = updatetwo[[1]]
    updatedelta = updatetwo[[2]]
    
    # Update beta parameter
    betatemp = updateBeta(updatebeta, betasd, updateaijmat, updatedelta,
                          GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow, genename, genenum, lastvaluelist,
                          oldcell, oldtime, newallhash, eachdelta,
                          yesalpha, noalpha)

    updatebeta = betatemp[[1]]
    betaaccept = betatemp[[2]]
    acc_beta = acc_beta + betaaccept

    # Update alpha parameter
    alphatemp = updateAlpha(updatealpha, alphasd, updateaijmat, updatedelta, updatebeta,
                            GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow, genename, genenum, lastvaluelist,
                            oldcell, oldtime, newallhash, eachdelta)

    updatealpha = alphatemp[[1]]
    alphaaccept = alphatemp[[2]]
    acc_alpha = acc_alpha + alphaaccept
    yesalpha = log(1 / (1 + exp(-updatealpha)))
    noalpha = log(1 / (1 + exp(updatealpha)))

    # Update lambda parameter
    lambdatemp = updateLambda(updatelambda, lambdasd, updateaijmat, alllogprior)
    updatelambda = lambdatemp[[1]]
    lambdaaccept = lambdatemp[[2]]
    acc_lambda = acc_lambda + lambdaaccept
    
    # Update kappa parameter
    katemp = updateKa(updateka, kasd, updateaijmat, genenum)
    updateka = katemp[[1]]
    kaaccept = katemp[[2]]
    acc_ka = acc_ka + kaaccept

    # Recalculate log-priors
    alllogprior = singleAijPrior(genenum, updatelambda, updateka)
    
    # Compute post-probabilities
    allprob[t_iter] = postProb(updateaijmat, updatedelta, updatebeta, updatealpha,
                               updatelambda, updateka, GoalDataset, InfoDataset,
                               statdiff, allchooserow, allstartrow,
                               genename, genenum, lastvaluelist, oldcell, oldtime, maxdelta, eachdelta, alllogprior, newallhash,
                               yesalpha, noalpha)
    
    # Store updated parameters in the sparse matrix
    updateaijvec = as.numeric(updateaijmat)
    updatedeltavec = as.numeric(updatedelta)
    allpara[t_iter, 1:(genenum^2)] = updateaijvec
    allpara[t_iter, (genenum^2 + 1):(2 * genenum^2)] = updatedeltavec
    allpara[t_iter, 2 * genenum^2 + 1] = updatebeta
    allpara[t_iter, 2 * genenum^2 + 2] = updatealpha
    allpara[t_iter, 2 * genenum^2 + 3] = updatelambda
    allpara[t_iter, 2 * genenum^2 + 4] = updateka
  }
  
  # Compile the results into a list for return
  result = list()
  allaccept = data.frame(acc_beta = acc_beta / iters,
                         acc_alpha = acc_alpha / iters,
                         acc_lambda = acc_lambda / iters,
                         acc_ka = acc_ka / iters)
  result[[1]] = allpara
  result[[2]] = allaccept
  result[[3]] = allprob
  
  return(result)
}

# (2) transExp ------------------------------------------------------------
# Define the function transExp to compute the exponential of a given numeric input
# while handling underflows typically occurring with large negative values.
# Inputs:
#   x: Numeric vector, contains values for which exponentials are to be calculated
# Output:
#   Numeric vector, contains the computed exponential values for each element in x
transExp = function(x){
  # Check if all values in x are less than -700, which can lead to underflow
  if(all(x < (-700))){
    # Adjust x by shifting it closer to zero to avoid underflow when computing exponentials
    x_max = max(x)  # Find the maximum value in x (least negative)
    x_diff = (-50) - x_max  # Calculate the adjustment needed to shift x values to -50
    y = x + x_diff  # Apply the adjustment
    a = exp(y)  # Compute the exponential of the adjusted values
    return(a)  # Return the adjusted exponential values
  }
  else{
    # Compute the exponential of x directly when there are no underflow risks
    a = exp(x)
    return(a)
  }
}

# (3) noZeroProb ---------------------------------------------------------
# Define the function noZeroProb to calculate the log probability of gene expression outcomes
# considering influences from other genes in a network.
# Inputs:
#   lastvaluelist: List of matrices, each containing last observed values for each gene
#   singlegene: String, the target gene for which probability is being calculated
#   deltaindex: Numeric vector, indices for delta time adjustments in the model
#   useaijvec: Numeric vector, interaction strengths between the target gene and others
#   nozeroindex: Numeric vector, indices of genes that have non-zero interactions with the target gene
#   yesalpha, noalpha: Numeric, log odds of observing a gene expression change given the state is unchanged
#   updatebeta: Numeric, regression coefficient reflecting the effect size of gene interactions
# Output:
#   Numeric, the calculated log probability of observing the gene expression data given the model parameters
noZeroProb = function(lastvaluelist, singlegene, deltaindex, useaijvec, nozeroindex,
                       yesalpha, noalpha, updatebeta){
  # Retrieve value matrices for the single gene of interest
  valuemat = lastvaluelist[[singlegene]]
  thisnow = valuemat[,1]  # Current gene expression values
  thislast = valuemat[,2]  # Previous gene expression values
  
  # Retrieve and combine value matrices for other influential genes based on nozeroindex
  revaluemat = lastvaluelist[[nozeroindex[1]]]
  relast = revaluemat[,deltaindex[1], drop=F]
  
  # If there are multiple influential genes, combine their value matrices
  if(length(nozeroindex) != 1){
    for(i in 2:length(nozeroindex)){
      revaluemat = lastvaluelist[[nozeroindex[i]]]
      relast = cbind(relast, revaluemat[,deltaindex[i], drop=F])
    }
  }

  # Calculate the influence score H by multiplying interaction strengths with gene values
  H = as.numeric(useaijvec[nozeroindex] %*% t(relast))
  
  # Calculate contributions to the probability from zeros in H
  zeroindex = which(H == 0)  # Indices where influence is zero
  zerostat = abs(thisnow[zeroindex] - thislast[zeroindex])  # Absolute differences where influence is zero
  zerostat = table(factor(zerostat, levels=c(0,1)))  # Frequency table of differences
  lnzero = zerostat[1] * yesalpha + zerostat[2] * noalpha  # Log probability for zero influence cases
  
  # Calculate contributions to the probability from non-zeros in H
  oneindex = which(H != 0)  # Indices where influence is non-zero
  Hbeta = updatebeta * H[oneindex]  # Scaled influence by beta for non-zero indices
  lnone = sum(2 * thisnow[oneindex] * Hbeta - log(exp(2 * Hbeta) + 1))  # Log probability for non-zero influence cases

  # Return the total log probability combining both zero and non-zero contributions
  return(lnzero + lnone)
}

# (4) singleLikelihood ---------------------------------------------------
# Define the function singleLikelihood to compute the log-likelihood for the expression of a single gene,
# taking into account its interactions with other genes and the impact of various model parameters.
# Inputs:
#   singlegene: Integer, index of the gene for which likelihood is being calculated
#   updateaijmat: Matrix, current estimates of aij
#   updatedelta: Matrix, current estimates of time delay
#   updatebeta: Numeric, continuous parameter
#   GoalDataset: DataFrame, target data for the MCMC optimization
#   InfoDataset: DataFrame, contains timing and cell information for each observation
#   statdiff: Matrix, statistics differential used in likelihood calculations
#   allchooserow: Vector, indices of rows chosen for analysis
#   allstartrow: Vector, indices of the start rows for each segment of data
#   genename: Vector, names of the genes
#   genenum: Integer, number of genes
#   lastvaluelist: List, last observed values for each gene
#   oldcell: Vector, cell identifiers corresponding to each row in InfoDataset
#   oldtime: Vector, timing information corresponding to each row in InfoDataset
#   eachdelta: Numeric, step size for delta calculations
#   newallhash: List, contains hashed tables of prior interactions
#   rowaijsum: Numeric or vector, sum of interaction strengths for each gene (or the single gene)
#   yesalpha, noalpha: Numeric, log odds parameters for gene expression state changes
# Output:
#   Numeric, the computed log-likelihood for the specified single gene
singleLikelihood = function(singlegene, updateaijmat, updatedelta, updatebeta,
                            GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow,
                            genename, genenum, lastvaluelist,
                            oldcell, oldtime, eachdelta,
                            newallhash, rowaijsum,
                            yesalpha, noalpha){

  # Handle input whether rowaijsum is a single value or a vector
  singlerowaijsum = if(length(rowaijsum) == 1) rowaijsum else rowaijsum[singlegene]
  
  # Extract the interaction vector for the single gene
  useaijvec = updateaijmat[singlegene,]
  
  # Calculate log-likelihood depending on whether the gene has any non-zero interactions
  if(singlerowaijsum == 0){
    # If no interactions, compute log-likelihood based on the statistical difference directly
    temp = statdiff[,singlegene]
    lny = temp[1] * yesalpha + temp[2] * noalpha
  }
  else{
    # Identify indices of non-zero interactions
    nozeroindex = which(useaijvec != 0)
    # Corresponding delta values for these interactions
    choosedelta = updatedelta[singlegene, nozeroindex]
    # Index these deltas to use in calculations
    deltaindex = round(choosedelta / eachdelta) + 1
    
    # Compute log-likelihood using the noZeroProb function
    lny = noZeroProb(lastvaluelist, singlegene, deltaindex, useaijvec, nozeroindex,
                     yesalpha, noalpha, updatebeta)
  }
  
  # Return the computed log-likelihood
  return(lny)
}

# (5) singleDelta --------------------------------------------------------
# Define the function singleDelta to compute the likelihood by updating the delta parameter
# for a specific pair of genes within a gene interaction network.
# Inputs:
#   singledelta: Numeric, the new delta value to be tested for the gene pair (i_iter, j_iter)
#   i_iter: Integer, index of the regulated gene in the gene pair
#   j_iter: Integer, index of the regulator gene in the gene pair
#   updateaijmat: Matrix, current estimates of aij
#   updatedelta: Matrix, current estimates of time delay
#   updatebeta: Numeric, regression coefficient reflecting the effect size of gene interactions
#   GoalDataset: DataFrame, target data for the MCMC optimization
#   InfoDataset: DataFrame, contains timing and cell information for each observation
#   statdiff: Matrix, statistics differential used in likelihood calculations
#   allchooserow: Vector, indices of rows chosen for analysis
#   allstartrow: Vector, indices of the start rows for each segment of data
#   genename: Vector, names of the genes
#   genenum: Integer, number of genes
#   lastvaluelist: List, last observed values for each gene
#   oldcell: Vector, cell identifiers corresponding to each row in InfoDataset
#   oldtime: Vector, timing information corresponding to each row in InfoDataset
#   eachdelta: Numeric, the increment between each possible delta value within the range from 0 to maxdelta
#   newallhash: List, contains hashed tables of prior interactions
#   singlerowaijsum: Numeric, sum of interaction strengths for the given row (gene)
#   yesalpha, noalpha: Numeric, log odds parameters for gene expression state changes
# Output:
#   Numeric, the computed likelihood for the gene pair under the new delta parameter
singleDelta = function(singledelta, i_iter, j_iter, updateaijmat, updatedelta, updatebeta,
                       GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow,
                       genename, genenum, lastvaluelist, oldcell, oldtime, eachdelta,
                       newallhash, singlerowaijsum, yesalpha, noalpha){
  
  # Set the new delta value for the specified gene pair in the updatedelta matrix
  updatedelta[i_iter, j_iter] = singledelta
  
  # Calculate the likelihood of the gene expression data under the new delta setting using the singleLikelihood function
  thisprob = singleLikelihood(i_iter, updateaijmat, updatedelta, updatebeta,
                              GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow,
                              genename, genenum, lastvaluelist,
                              oldcell, oldtime, eachdelta,
                              newallhash, singlerowaijsum,
                              yesalpha, noalpha)
  
  # Return the computed likelihood
  return(thisprob)
}

# (6) newTwoProb ---------------------------------------------------------
# Define the function newTwoProb to compute the total log-probability for a gene interaction model
# after potentially modifying an interaction strength and assessing the resulting model fit.
# Inputs:
#   newaij: Numeric, proposed new interaction strength for a gene pair (i_iter, j_iter)
#   alldelta: Numeric vector, range of delta values for the interaction strength adjustments
#   updateaijmat: Matrix, current estimates of aij
#   updatedelta: Matrix, current estimates of time delay
#   updatebeta: Numeric, regression coefficient reflecting the effect size of gene interactions
#   i_iter: Integer, index of the regulated gene in the gene pair
#   j_iter: Integer, index of the regulator gene in the gene pair
#   GoalDataset: DataFrame, target data for the MCMC optimization
#   InfoDataset: DataFrame, contains timing and cell information for each observation
#   statdiff: Matrix, statistics differential used in likelihood calculations
#   allchooserow: Vector, indices of rows chosen for analysis
#   allstartrow: Vector, indices of the start rows for each segment of data
#   genename: Vector, names of the genes
#   genenum: Integer, number of genes
#   lastvaluelist: List, last observed values for each gene
#   oldcell: Vector, cell identifiers corresponding to each row in InfoDataset
#   oldtime: Vector, timing information corresponding to each row in InfoDataset
#   eachdelta: Numeric, the increment between each possible delta value within the range from 0 to maxdelta
#   alllogprior: Vector, log-prior probabilities associated with each gene's interaction parameters
#   newallhash: List, contains hashed tables of prior interactions
#   logzerodelta: Numeric, log of the prior probability when the interaction strength is zero
#   logdelta: Numeric, log of the prior probability when there is non-zero interaction strength
#   yesalpha, noalpha: Numeric, log odds parameters for gene expression state changes
# Output:
#   Numeric, total log-probability of observing the gene expression data under the updated model parameters
newTwoProb = function(newaij, alldelta, updateaijmat, updatedelta, updatebeta,
                      i_iter, j_iter,
                      GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow,
                      genename, genenum, lastvaluelist,
                      oldcell, oldtime, eachdelta,
                      alllogprior, newallhash, logzerodelta, logdelta,
                      yesalpha, noalpha){

  if(newaij == 0){
    # If the new aij is zero, set it and calculate the probability of the model
    updateaijmat[i_iter, j_iter] = 0
    singlerowaijsum = sum(abs(updateaijmat[i_iter,]))
    allprob = singleLikelihood(i_iter, updateaijmat, updatedelta, updatebeta,
                               GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow,
                               genename, genenum, lastvaluelist,
                               oldcell, oldtime, eachdelta,
                               newallhash, singlerowaijsum,
                               yesalpha, noalpha)
    # Return the combined probability including the inverse prior and log probability of zero delta
    return(allprob + invAijPrior(updateaijmat, alllogprior) + logzerodelta)
  }
  else{
    # If there is a non-zero new aij, update the model and calculate probabilities for each delta
    updateaijmat[i_iter, j_iter] = newaij
    singlerowaijsum = 1
    allprob = sapply(alldelta, singleDelta, i_iter, j_iter, updateaijmat, updatedelta, updatebeta,
                     GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow,
                     genename, genenum, lastvaluelist, oldcell, oldtime, eachdelta,
                     newallhash, singlerowaijsum, yesalpha, noalpha)
    # Return the combined probability including the inverse prior and log probability of non-zero delta
    return(allprob + invAijPrior(updateaijmat, alllogprior) + logdelta)
  }
}

# (7) rowAij ------------------------------------------------------------------
# Define the function rowAij to update interaction strengths and delta values for a specific gene (row) based on 
# their probability contributions to model fit.
# Inputs:
#   i_iter: Integer, index of the gene for which interactions are being updated
#   newaijindex: List, indices of genes that have interactions with the gene specified by i_iter
#   allaij: Numeric vector, all possible values of aij being considered
#   alldelta: Numeric vector, all possible values of delta being considered
#   updateaijmat: Matrix, current estimates of aij
#   updatedelta: Matrix, current estimates of time delay
#   updatebeta: Numeric, continuous parameter
#   GoalDataset: DataFrame, target data for the MCMC optimization
#   InfoDataset: DataFrame, contains timing and cell information for each observation
#   statdiff: Matrix, statistics differential used in likelihood calculations
#   allchooserow: Vector, indices of rows chosen for analysis
#   allstartrow: Vector, indices of the start rows for each segment of data
#   genename: Vector, names of the genes
#   genenum: Integer, number of genes
#   lastvaluelist: List, last observed values for each gene
#   oldcell: Vector, cell identifiers corresponding to each row in InfoDataset
#   oldtime: Vector, timing information corresponding to each row in InfoDataset
#   eachdelta: Numeric, the increment between each possible delta value within the range from 0 to maxdelta
#   alllogprior: Vector, log-prior probabilities associated with each gene's interaction parameters
#   newallhash: List, contains hashed tables of prior interactions
#   logzerodelta: Numeric, log of the prior probability when the interaction strength is zero
#   logdelta: Numeric, log of the prior probability when there is non-zero interaction strength
#   allnewtwo: Matrix, contains all combinations of new possible values for aij and corresponding deltas
#   yesalpha, noalpha: Numeric, log odds parameters for gene expression state changes
# Output:
#   List containing updated rows of the interaction matrix and delta values
rowAij = function(i_iter, newaijindex, allaij, alldelta, updateaijmat, updatedelta, updatebeta,
                     GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow,
                     genename, genenum, lastvaluelist,
                     oldcell, oldtime, eachdelta,
                     alllogprior, newallhash, logzerodelta, logdelta, allnewtwo,
                     yesalpha, noalpha){

  # Loop over each gene index that interacts with the gene specified by i_iter
  for(j_iter in newaijindex[[i_iter]]){
    # Calculate the probability for each potential new interaction strength
    newaij_prob = sapply(allaij, newTwoProb, alldelta, updateaijmat, updatedelta, updatebeta,
                         i_iter, j_iter,
                         GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow,
                         genename, genenum, lastvaluelist,
                         oldcell, oldtime, eachdelta,
                         alllogprior, newallhash, logzerodelta, logdelta,
                         yesalpha, noalpha)

    # Convert log probabilities to probabilities using the transExp function
    newaij_prob = as.numeric(newaij_prob)
    newaij_prob = transExp(newaij_prob)
    
    # Randomly select a new interaction strength based on the calculated probabilities
    choosenum = sample(1:nrow(allnewtwo), 1, FALSE, newaij_prob)
    
    # Update the interaction matrix and delta values based on the selected new interaction strength
    updateaijmat[i_iter, j_iter] = allnewtwo[choosenum, 1]
    updatedelta[i_iter, j_iter] = allnewtwo[choosenum, 2]
  }

  # Compile the updated interaction strengths and deltas into a list for the specified row
  templist = list(rowaij = updateaijmat[i_iter,],
                  rowdelta = updatedelta[i_iter,])
  
  # Return the updated row of interaction strengths and corresponding deltas
  return(templist)
}


# (8) updateTwoMat -------------------------------------------------------
# Define the function updateTwoMat to perform systematic updates to the interaction matrices and delta values
# across all genes in the gene interaction network model.
# Inputs:
#   updateaijmat: Matrix, current estimates of aij
#   updatebeta: Numeric, continuous parameter
#   updatedelta: Matrix, current estimates of time delay
#   maxdelta: Numeric, maximum value for delta in the range to consider
#   eachdelta: Numeric, increment for delta in the range to consider
#   alllogprior: Vector, log-prior probabilities for gene interactions
#   GoalDataset: DataFrame, target data for the MCMC optimization
#   InfoDataset: DataFrame, contains timing and cell information for each observation
#   statdiff: Matrix, statistical differential used in likelihood calculations
#   allchooserow: Vector, indices of rows chosen for analysis
#   allstartrow: Vector, indices of the start rows for each segment of data
#   genenum: Integer, number of genes
#   genename: Vector, names of the genes
#   lastvaluelist: List, last observed values for each gene
#   newallhash: List, contains hashed tables of prior interactions
#   oldcell: Vector, cell identifiers corresponding to each row in InfoDataset
#   oldtime: Vector, timing information corresponding to each row in InfoDataset
#   yesalpha, noalpha: Numeric, log odds parameters for gene expression state changes
#   newaijindex: List, indices of genes that have interactions with each gene
# Output:
#   List containing updated matrices of interaction strengths and delta values
updateTwoMat = function(updateaijmat, updatebeta, updatedelta, maxdelta, eachdelta,
                        alllogprior, GoalDataset, InfoDataset, statdiff,
                        allchooserow, allstartrow, genenum, genename, lastvaluelist,
                        newallhash, oldcell, oldtime,
                        yesalpha, noalpha, newaijindex){

  # Calculate the total number of delta steps to consider
  alltwonum = round(maxdelta / eachdelta) + 1
  
  # Define all possible interaction strengths and their corresponding deltas
  allaij = c(0, 1, -1)  # Possible interaction strength values
  alldelta = seq(0, maxdelta, by = eachdelta)  # Range of delta values
  
  # Create a matrix of all possible new interactions and deltas
  allnewtwo = matrix(c(rep(allaij, each = 6), rep(alldelta, 3)), ncol = 2)
  
  # Log probabilities for zero and non-zero deltas
  logzerodelta = c(log(0.999), rep(log(0.0002), 5))
  logdelta = c(rep(log(1 / 6), 6))
  
  # Apply rowAij function to each gene to update interactions and deltas
  templist = lapply(1:genenum, rowAij, newaijindex, allaij, alldelta, updateaijmat, updatedelta, updatebeta,
                    GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow,
                    genename, genenum, lastvaluelist,
                    oldcell, oldtime, eachdelta,
                    alllogprior, newallhash, logzerodelta, logdelta, allnewtwo,
                    yesalpha, noalpha)
  
  # Extract updated interaction matrices and delta values
  updateaijmat = t(sapply(1:genenum, function(x) templist[[x]][[1]]))
  updatedelta = t(sapply(1:genenum, function(x) templist[[x]][[2]]))
  
  # Compile updated interaction and delta matrices into a list for return
  updatetwo = list()
  updatetwo[[1]] = updateaijmat
  updatetwo[[2]] = updatedelta
  return(updatetwo)
}

# (9) piBeta -------------------------------------------------------------
# Define the function piBeta to calculate the posterior probability of the beta parameter
# Inputs:
#   beta: Numeric, the value of beta to evaluate
#   updateaijmat: Matrix, current estimates of aij
#   updatedelta: Matrix, current estimates of time delay
#   GoalDataset: DataFrame, target data for the MCMC optimization
#   InfoDataset: DataFrame, contains timing and cell information for each observation
#   statdiff: Matrix, statistical differential used in likelihood calculations
#   allchooserow: Vector, indices of rows chosen for analysis
#   allstartrow: Vector, indices of the start rows for each segment of data
#   genename: Vector, names of the genes
#   oldcell: Vector, cell identifiers corresponding to each row in InfoDataset
#   oldtime: Vector, timing information corresponding to each row in InfoDataset
#   newallhash: List, contains hashed tables of prior interactions
#   eachdelta: Numeric, step size for delta calculations
#   genenum: Integer, number of genes
#   lastvaluelist: List, last observed values for each gene
#   yesalpha, noalpha: Numeric, log odds parameters for gene expression state changes
# Output:
#   Numeric, the log of the posterior probability of beta given the data and model
piBeta = function(beta, updateaijmat, updatedelta,
                   GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow, genename,
                   oldcell, oldtime, newallhash, eachdelta,
                   genenum, lastvaluelist, yesalpha, noalpha){

  # Sum the absolute values of interaction strengths for each gene to use in likelihood calculations
  rowaijsum = apply(abs(updateaijmat), 1, sum)
  
  # Compute the total log-likelihood across all genes for the given beta
  allprob = sum(sapply(1:genenum, singleLikelihood, updateaijmat, updatedelta, beta,
                       GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow, genename, genenum, lastvaluelist,
                       oldcell, oldtime, eachdelta,
                       newallhash, rowaijsum,
                       yesalpha, noalpha))
  
  # Compute the log of the prior probability of beta assuming a gamma distribution
  I = dgamma(beta, shape = 100, rate = 100, log = TRUE)
  
  # Return the sum of the log-likelihood and the log prior, representing the log of the posterior probability
  return(allprob + I)
}

# (10) updateBeta ---------------------------------------------------------
# Define the function updateBeta to perform a Metropolis-Hastings update of the beta parameter.
# Inputs:
#   updatebeta: Numeric, current value of beta
#   betasd: Numeric, standard deviation used for the proposal distribution of beta
#   updateaijmat: Matrix, current estimates of aij
#   updatedelta: Matrix, current estimates of time delay
#   GoalDataset: DataFrame, target data for model fitting
#   InfoDataset: DataFrame, contains timing and cell information for each observation
#   statdiff: Matrix, statistical differences used in model evaluation
#   allchooserow: Vector, indices of rows chosen for analysis
#   allstartrow: Vector, indices of the start rows for each data segment
#   genename: Vector, names of the genes
#   genenum: Integer, number of genes
#   lastvaluelist: List, contains last observed values for each gene
#   oldcell: Vector, cell identifiers corresponding to each row in InfoDataset
#   oldtime: Vector, timing information corresponding to each row in InfoDataset
#   eachdelta: Numeric, increment for delta in the model
#   yesalpha, noalpha: Numeric, log odds of observing a gene expression change given the state is unchanged
# Output:
#   List containing the updated beta value and an indicator (0 or 1) of whether the update was accepted
updateBeta = function(updatebeta, betasd, updateaijmat, updatedelta,
                       GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow, genename, genenum, lastvaluelist,
                       oldcell, oldtime, newallhash, eachdelta,
                       yesalpha, noalpha){

  # Propose a new value for beta from a truncated normal distribution limited to positive values
  newbeta = rtruncnorm(1, a = 0, b = Inf, mean = updatebeta, sd = betasd)
  
  # Calculate the log probability ratio of the new beta to the current beta
  ratio = piBeta(newbeta, updateaijmat, updatedelta,
                 GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow, genename,
                 oldcell, oldtime, newallhash, eachdelta,
                 genenum, lastvaluelist, yesalpha, noalpha) -
          piBeta(updatebeta, updateaijmat, updatedelta,
                 GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow, genename,
                 oldcell, oldtime, newallhash, eachdelta,
                 genenum, lastvaluelist, yesalpha, noalpha) +
          log(dtruncnorm(updatebeta, a = 0, b = Inf, mean = newbeta, sd = betasd)) -
          log(dtruncnorm(newbeta, a = 0, b = Inf, mean = updatebeta, sd = betasd))
  
  # Determine if the new beta should be accepted based on the log probability ratio
  accept = c(0, ratio)
  result = list()
  if (log(runif(1)) <= min(accept)) {
    result[[1]] = newbeta  # Accept the new beta value
    result[[2]] = 1  # Indicate acceptance
  } else {
    result[[1]] = updatebeta  # Retain the old beta value
    result[[2]] = 0  # Indicate rejection
  }
  return(result)
}

# (11) piAlpha ------------------------------------------------------------
# Define the function piAlpha to calculate the posterior probability of the alpha parameter
# Inputs:
#   alpha: Numeric, the value of alpha to evaluate
#   updateaijmat: Matrix, current estimates of aij
#   updatedelta: Matrix, current estimates of time delay
#   updatebeta: Numeric, continuous parameter
#   GoalDataset: DataFrame, target data for the MCMC optimization
#   InfoDataset: DataFrame, contains timing and cell information for each observation
#   statdiff: Matrix, statistical differential used in likelihood calculations
#   allchooserow: Vector, indices of rows chosen for analysis
#   allstartrow: Vector, indices of the start rows for each segment of data
#   genename: Vector, names of the genes
#   oldcell: Vector, cell identifiers corresponding to each row in InfoDataset
#   oldtime: Vector, timing information corresponding to each row in InfoDataset
#   newallhash: List, contains hashed tables of prior interactions
#   eachdelta: Numeric, step size for delta calculations
#   genenum: Integer, number of genes
#   lastvaluelist: List, last observed values for each gene
# Output:
#   Numeric, the log of the posterior probability of alpha given the data and model
piAlpha = function(alpha, updateaijmat, updatedelta, updatebeta,
                   GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow, genename,
                   oldcell, oldtime, newallhash, eachdelta,
                   genenum, lastvaluelist){

  # Compute the log odds for the given alpha
  yesalpha = log(1 / (1 + exp(-alpha)))
  noalpha = log(1 / (1 + exp(alpha)))
  
  # Sum the absolute values of interaction strengths for each gene to use in likelihood calculations
  rowaijsum = apply(abs(updateaijmat), 1, sum)
  
  # Compute the total log-likelihood across all genes for the given alpha
  allprob = sum(sapply(1:genenum, singleLikelihood, updateaijmat, updatedelta, updatebeta,
                       GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow, genename, genenum, lastvaluelist,
                       oldcell, oldtime, eachdelta,
                       newallhash, rowaijsum,
                       yesalpha, noalpha))
  
  # Compute the log of the prior probability of alpha assuming a gamma distribution
  I = dgamma(alpha, shape = 1, rate = 10, log = TRUE)
  
  # Return the sum of the log-likelihood and the log prior, representing the log of the posterior probability
  return(allprob + I)
}


# (12) updateAlpha --------------------------------------------------------
# Define the function updateAlpha to perform a Metropolis-Hastings update of the alpha parameter.
# Inputs:
#   updatealpha: Numeric, current value of alpha
#   alphasd: Numeric, standard deviation used for the proposal distribution of alpha
#   updateaijmat: Matrix, current estimates of aij
#   updatedelta: Matrix, current estimates of time delay
#   updatebeta: Numeric, continuous parameter
#   GoalDataset: DataFrame, target data for the MCMC optimization
#   InfoDataset: DataFrame, contains timing and cell information for each observation
#   statdiff: Matrix, statistical differences used in model evaluation
#   allchooserow: Vector, indices of rows chosen for analysis
#   allstartrow: Vector, indices of the start rows for each data segment
#   genename: Vector, names of the genes
#   genenum: Integer, number of genes
#   lastvaluelist: List, contains last observed values for each gene
#   oldcell: Vector, cell identifiers corresponding to each row in InfoDataset
#   oldtime: Vector, timing information corresponding to each row in InfoDataset
#   newallhash: List, contains hashed tables of prior interactions
#   eachdelta: Numeric, increment for delta in the model
# Output:
#   List containing the updated alpha value and an indicator (0 or 1) of whether the update was accepted
updateAlpha = function(updatealpha, alphasd, updateaijmat, updatedelta, updatebeta,
                       GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow, genename, genenum, lastvaluelist,
                       oldcell, oldtime, newallhash, eachdelta){

  # Propose a new value for alpha from a truncated normal distribution limited to positive values
  newalpha = rtruncnorm(1, a = 0, b = Inf, mean = updatealpha, sd = alphasd)
  
  # Calculate the log probability ratio of the new alpha to the current alpha
  ratio = piAlpha(newalpha, updateaijmat, updatedelta, updatebeta,
                  GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow, genename,
                  oldcell, oldtime, newallhash, eachdelta,
                  genenum, lastvaluelist) -
          piAlpha(updatealpha, updateaijmat, updatedelta, updatebeta,
                  GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow, genename,
                  oldcell, oldtime, newallhash, eachdelta,
                  genenum, lastvaluelist) +
          log(dtruncnorm(updatealpha, a = 0, b = Inf, mean = newalpha, sd = alphasd)) -
          log(dtruncnorm(newalpha, a = 0, b = Inf, mean = updatealpha, sd = alphasd))
  
  # Determine if the new alpha should be accepted based on the log probability ratio
  accept = c(0, ratio)
  result = list()
  if (log(runif(1)) <= min(accept)) {
    result[[1]] = newalpha  # Accept the new alpha value
    result[[2]] = 1  # Indicate acceptance
  } else {
    result[[1]] = updatealpha  # Retain the old alpha value
    result[[2]] = 0  # Indicate rejection
  }
  return(result)
}


# (13) piLambda -----------------------------------------------------------
# Define the function piLambda to calculate the posterior probability of the lambda parameter
# Inputs:
#   lambda: Numeric, the value of lambda to evaluate
#   updateaijmat: Matrix, current estimates of aij
#   alllogprior: Vector, log-prior probabilities for gene interactions
# Output:
#   Numeric, the log of the posterior probability of lambda given the data and model
piLambda = function(lambda, updateaijmat, alllogprior){
  # Compute the total log-prior probability across all genes for the given lambda
  allprob = invAijPrior(updateaijmat, alllogprior)
  
  # Compute the prior probability of lambda assuming a gamma distribution
  I = dgamma(lambda, shape = 490, rate = 70)
  
  # Return the sum of the log-prior and the log of the prior probability of lambda
  return(allprob + log(I))
}

# (14) updateLambda -------------------------------------------------------
# Define the function updateLambda to perform a Metropolis-Hastings update of the lambda parameter.
# Inputs:
#   updatelambda: Numeric, current value of lambda
#   lambdasd: Numeric, standard deviation used for the proposal distribution of lambda
#   updateaijmat: Matrix, current estimates of aij
#   alllogprior: Vector, log-prior probabilities for gene interactions
# Output:
#   List containing the updated lambda value and an indicator (0 or 1) of whether the update was accepted
updateLambda = function(updatelambda, lambdasd, updateaijmat, alllogprior){
  
  # Propose a new value for lambda from a truncated normal distribution limited to positive values
  newlambda = rtruncnorm(1, a = 0, b = Inf, mean = updatelambda, sd = lambdasd)
  
  # Calculate the log probability ratio of the new lambda to the current lambda
  ratio = piLambda(newlambda, updateaijmat, alllogprior) -
          piLambda(updatelambda, updateaijmat, alllogprior) +
          log(dtruncnorm(updatelambda, a = 0, b = Inf, mean = newlambda, sd = lambdasd)) -
          log(dtruncnorm(newlambda, a = 0, b = Inf, mean = updatelambda, sd = lambdasd))
  
  # Determine if the new lambda should be accepted based on the log probability ratio
  accept = c(0, ratio)
  result = list()
  if (log(runif(1)) <= min(accept)) {
    result[[1]] = newlambda  # Accept the new lambda value
    result[[2]] = 1  # Indicate acceptance
  } else {
    result[[1]] = updatelambda  # Retain the old lambda value
    result[[2]] = 0  # Indicate rejection
  }
  return(result)
}

# (15) piKa ---------------------------------------------------------------
# Define the function piKa to calculate the posterior probability of the ka parameter
# used in a gene interaction model.
# Inputs:
#   updateka: Numeric, the value of ka to evaluate
#   updateaijmat: Matrix, current estimates of aij
#   genenum: Integer, number of genes
# Output:
#   Numeric, the log of the posterior probability of ka given the data and model
piKa = function(updateka, updateaijmat, genenum){
  
  # Calculate the sum of absolute interaction strengths for each row (gene)
  aijrowsum = apply(abs(updateaijmat), 1, sum)
  
  # For each gene, if it has any interactions (non-zero sum), set the sum to 1; otherwise, keep it 0
  aijrowsum[aijrowsum > 0] = 1
  
  # Calculate the number of genes with at least one interaction (onesum) and no interactions (zerosum)
  onesum = sum(aijrowsum)
  zerosum = genenum - onesum
  
  # Compute the log probability of the interaction states given ka
  allprob = log((updateka^onesum) * ((1 - updateka)^zerosum))
  
  # Compute the prior probability of ka assuming a beta distribution
  I = dbeta(updateka, shape1 = 24, shape2 = 6, log = TRUE)
  
  # Return the sum of the log likelihood and the log prior, representing the log of the posterior probability
  return(allprob + I)
}


# (16) updateKa -----------------------------------------------------------
# Define the function updateKa to perform a Metropolis-Hastings update of the ka parameter.
# Inputs:
#   updateka: Numeric, current value of ka
#   kasd: Numeric, standard deviation used for the proposal distribution of ka
#   updateaijmat: Matrix, current estimates of aij
#   genenum: Integer, number of genes
# Output:
#   List containing the updated ka value and an indicator (0 or 1) of whether the update was accepted
updateKa = function(updateka, kasd, updateaijmat, genenum){
  
  # Propose a new value for ka from a truncated normal distribution limited to the interval [0, 1]
  newka = rtruncnorm(1, a = 0, b = 1, mean = updateka, sd = kasd)
  
  # Calculate the log probability ratio of the new ka to the current ka
  ratio = piKa(newka, updateaijmat, genenum) -
          piKa(updateka, updateaijmat, genenum) +
          log(dtruncnorm(updateka, a = 0, b = 1, mean = newka, sd = kasd)) -
          log(dtruncnorm(newka, a = 0, b = 1, mean = updateka, sd = kasd))
  
  # Determine if the new ka should be accepted based on the log probability ratio
  accept = c(0, ratio)
  result = list()
  if (log(runif(1)) <= min(accept)) {
    result[[1]] = newka  # Accept the new ka value
    result[[2]] = 1  # Indicate acceptance
  } else {
    result[[1]] = updateka  # Retain the old ka value
    result[[2]] = 0  # Indicate rejection
  }
  return(result)
}

# (17) postProb -----------------------------------------------------------
# Define the function postProb to calculate the posterior probability of the model parameters
# used in a gene interaction model.
# Inputs:
#   updateaijmat: Matrix, current estimates of aij
#   updatedelta: Matrix, current estimates of time delay
#   updatebeta: Numeric, continuous parameter
#   updatealpha: Numeric, parameter influencing the probability of gene expression state changes
#   updatelambda: Numeric, parameter influencing the prior probability of gene interactions
#   updateka: Numeric, parameter influencing the probability of gene interactions
#   GoalDataset: DataFrame, target data for the MCMC optimization
#   InfoDataset: DataFrame, contains timing and cell information for each observation
#   statdiff: Matrix, statistical differential used in likelihood calculations
#   allchooserow: Vector, indices of rows chosen for analysis
#   allstartrow: Vector, indices of the start rows for each segment of data
#   genename: Vector, names of the genes
#   genenum: Integer, number of genes
#   lastvaluelist: List, last observed values for each gene
#   oldcell: Vector, cell identifiers corresponding to each row in InfoDataset
#   oldtime: Vector, timing information corresponding to each row in InfoDataset
#   maxdelta: Numeric, maximum value for delta in the range to consider
#   eachdelta: Numeric, increment for delta in the range to consider
#   alllogprior: Vector, log-prior probabilities for gene interactions
#   newallhash: List, contains hashed tables of prior interactions
#   yesalpha, noalpha: Numeric, log odds parameters for gene expression state changes
# Output:
#   Numeric, the log of the posterior probability of the model parameters given the data and model
postProb = function(updateaijmat, updatedelta, updatebeta, updatealpha,
                    updatelambda, updateka, GoalDataset, InfoDataset,
                    statdiff, allchooserow, allstartrow,
                    genename, genenum, lastvaluelist, oldcell, oldtime,
                    maxdelta, eachdelta, alllogprior, newallhash,
                    yesalpha, noalpha){

  # Calculate the sum of absolute interaction strengths for each row (gene)
  allrowaijsum = apply(abs(updateaijmat), 1, sum)
  
  # Compute the total log-likelihood across all genes
  likelihood_prob = sum(sapply(1:genenum, singleLikelihood, updateaijmat, updatedelta, updatebeta,
                               GoalDataset, InfoDataset, statdiff, allchooserow, allstartrow, genename, genenum, lastvaluelist,
                               oldcell, oldtime, eachdelta,
                               newallhash, allrowaijsum,
                               yesalpha, noalpha))
  
  # Compute the log-prior probability of the interaction matrix
  aijmat_prob = invAijPrior(updateaijmat, alllogprior)

  # Compute the log-prior probability of the delta matrix
  delta_prob = invDeltaPrior(updatedelta, updateaijmat, genenum)

  # Compute the log-prior probability of beta assuming a gamma distribution
  beta_prob = dgamma(updatebeta, shape = 100, rate = 100, log = TRUE)

  # Compute the log-prior probability of alpha assuming a gamma distribution
  alpha_prob = dgamma(updatealpha, shape = 1, rate = 10, log = TRUE)

  # Compute the log-prior probability of lambda assuming a gamma distribution
  lambda_prob = dgamma(updatelambda, shape = 490, rate = 70, log = TRUE)

  # Compute the log-prior probability of ka assuming a beta distribution
  ka_prob = dbeta(updateka, shape1 = 24, shape2 = 6, log = TRUE)
  
  # Sum all log-probabilities to get the total log-posterior probability
  allprob = likelihood_prob + aijmat_prob + delta_prob + beta_prob + alpha_prob + lambda_prob + ka_prob

  return(allprob)
}

# (18) realMCMC -----------------------------------------------------------
# Define the function realMCMC to perform MCMC simulations for estimating gene interaction parameters.
# Inputs:
#   subtree: Character, name of the subtree being analyzed
#   p: Numeric, specific parameter value used in setting the random seed
#   mcmcnum: Integer, number of MCMC chains to run in parallel
#   num_maxiters: Integer, maximum number of MCMC iterations
#   iters: Integer, number of iterations per MCMC run
#   betasd, alphasd, lambdasd, kasd: Numeric, standard deviations for proposal distributions of respective parameters
#   maxdelta_time: Numeric, maximum delta time value
#   realpath: Character, path to the directory containing the real data files
#   realresultpath: Character, path to the directory containing the result files
# Output:
#   None, but saves the MCMC results and parameter estimates to files

realMCMC = function(subtree, p, mcmcnum, num_maxiters, iters,
                    betasd, alphasd, lambdasd, kasd,
                    maxdelta_time, realpath, realresultpath){

  # Set the random seed for reproducibility
  set.seed(p + 500)

  # Initialize the total number of iterations
  alliters = 0

  # Initialize parallel processing
  sfInit(parallel = TRUE, cpus = mcmcnum)
  sfLibrary(hash)
  sfLibrary(truncnorm)
  sfLibrary(Matrix)
  sfLibrary(dplyr)
  sfLibrary(tidyr)
  sfLibrary(readr)
  sfLibrary(stringr)
  sfLibrary(purrr)
  sfSource("./code/mainMCMCFun.R")
  sfSource("./code/priorsFun.R")
  sfSource("./code/estimateFun.R")

  # Load necessary data
  load(file = paste(realpath, subtree, "_simu_info.RData", sep = ""))
  Mc = simu_info$Mc
  eachdelta = 1 / (round(Mc) - 1)
  maxdelta = maxdelta_time * eachdelta
  load(paste(realpath, "mer_", subtree, "_", maxdelta_time, "_allchooserow.RData", sep = ""))
  load(paste(realpath, "mer_", subtree, "_", maxdelta_time, "_allstartrow.RData", sep = ""))
  load(paste(realpath, "mer_", subtree, "_", maxdelta_time, "_origenehash.RData", sep = ""))
  newdata = read_csv(file = paste(realpath, "mer_", subtree, "_", maxdelta_time, "_selectdata.csv", sep = ""), show_col_types = FALSE)
  finalinfo = newdata[, 1:6]
  finaldata = newdata[, -c(1:6)]
  oldtime = finalinfo$time
  oldcell = finalinfo$cell
  selectgene = names(finaldata)
  genenum = length(selectgene)
  aijindexlist = lapply(1:genenum, function(x, genenum) list = setdiff(1:genenum, x), genenum)
  names(aijindexlist) = selectgene
  load(file = paste(realpath, "mer_", subtree, "_newaijindex_", 1, "_.RData", sep = ""))
  load(paste(realpath, "mer_", subtree, "_", 1, "_diff", ".RData", sep = ""))
  load(paste(realpath, subtree, "_", p, "_lastvaluelist.RData", sep = ""))
  finaldata = map_dfc(finaldata, function(x) replace_na(x, 0))

  # Function to run the first MCMC chain
  firstsinglemcmc = function(i, p, genenum, lastvaluelist, iters, selectgene, newallhash, aijindexlist,
                             finaldata, finalinfo, statdiff, allchooserow, allstartrow, newaijindex,
                             betasd, alphasd, lambdasd, kasd, maxdelta, eachdelta){

    # Initialize parameters
    updatebeta = rgamma(1, shape = 100, rate = 100)
    updatealpha = rgamma(1, shape = 64, rate = 80)
    updatelambda = rgamma(1, shape = 490, rate = 70)
    updateka = rbeta(1, shape1 = 24, shape2 = 6)

    # Generate initial interaction matrix and delta matrix
    aijmat = geneAij(genenum, updatelambda, updateka, aijindexlist)
    updateaijmat = aijmat[, -(genenum + 1)]
    for(j in 1:nrow(updateaijmat)){
      zeroindex = setdiff(aijindexlist[[j]], newaijindex[[j]])
      updateaijmat[j, zeroindex] = 0
    }
    updatedelta = geneDelta(updateaijmat, maxdelta, eachdelta, genenum)

    # Run the MCMC chain
    chainresult = mainMCMC(iters, genenum, selectgene, lastvaluelist, newallhash,
                           finaldata, finalinfo, statdiff, allchooserow, allstartrow, newaijindex,
                           updateaijmat, updatebeta, updatealpha, updatelambda, updateka, updatedelta,
                           betasd, alphasd, lambdasd, kasd, maxdelta, eachdelta)

    # Extract results from the chain
    chain = chainresult[[1]]
    prob = chainresult[[3]]
    updateaijmat = matrix(chain[iters, 1:(genenum^2)], genenum, genenum)
    updatedelta = matrix(chain[iters, (genenum^2 + 1):(2 * genenum^2)], genenum, genenum)
    updatebeta = chain[iters, 2 * genenum^2 + 1]
    updatealpha = chain[iters, 2 * genenum^2 + 2]
    updatelambda = chain[iters, 2 * genenum^2 + 3]
    updateka = chain[iters, 2 * genenum^2 + 4]

    # Return the results
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

  # Initialize lists to store results
  allchain = list()
  allprob = list()
  allpara = list()

  # Run the first set of MCMC chains
  firstlist = sfLapply(1:mcmcnum, firstsinglemcmc, p, genenum, lastvaluelist, iters, selectgene, newallhash, aijindexlist,
                       finaldata, finalinfo, statdiff, allchooserow, allstartrow, newaijindex,
                       betasd, alphasd, lambdasd, kasd, maxdelta, eachdelta)

  # Process the results of the first set of MCMC chains
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
    allchain[[i]] = sparseMatrix(i = 1, j = 1, x = 1, dims = c(iters, 2 * genenum^2 + 4))
    allchain[[i]][, 1:(genenum^2)] = chain[, 1:(genenum^2)]
    allchain[[i]][, (genenum^2 + 1):(2 * genenum^2)] = chain[, (genenum^2 + 1):(2 * genenum^2)]
    allchain[[i]][, 2 * genenum^2 + 1] = chain[, 2 * genenum^2 + 1]
    allchain[[i]][, 2 * genenum^2 + 2] = chain[, 2 * genenum^2 + 2]
    allchain[[i]][, 2 * genenum^2 + 3] = chain[, 2 * genenum^2 + 3]
    allchain[[i]][, 2 * genenum^2 + 4] = chain[, 2 * genenum^2 + 4]
    allprob[[i]] = thisprob
  }

  # Save the results of the first set of MCMC chains
  save(allchain, file = paste(realresultpath, subtree, "_", p, "_result", ".RData", sep = ""))
  save(allprob, file = paste(realresultpath, subtree, "_", p, "_allprob", ".RData", sep = ""))
  choosepara = maxEstimate(allchain, allprob, genenum)
  preaij = choosepara[[1]]
  predelta = choosepara[[2]]
  preparas = choosepara[[3]]
  write.csv(preaij, file = paste(realresultpath, subtree, "_", p, "_preaij.csv", sep = ""), row.names = FALSE)
  write.csv(predelta, file = paste(realresultpath, subtree, "_", p, "_predelta.csv", sep = ""), row.names = FALSE)
  write.csv(preparas, file = paste(realresultpath, subtree, "_", p, "_preparas.csv", sep = ""), row.names = FALSE)
  alliters = alliters + iters

  # Function to run subsequent MCMC chains
  eachsinglemcmc = function(i, p, allpara, genenum, lastvaluelist, iters, selectgene, newallhash,
                            finaldata, finalinfo, statdiff, allchooserow, allstartrow, newaijindex,
                            betasd, alphasd, lambdasd, kasd, maxdelta, eachdelta){

    # Initialize parameters from previous run
    updatebeta = allpara[[i]][["updatebeta"]]
    updatealpha = allpara[[i]][["updatealpha"]]
    updatelambda = allpara[[i]][["updatelambda"]]
    updateka = allpara[[i]][["updateka"]]
    updateaijmat = allpara[[i]][["updateaijmat"]]
    updatedelta = allpara[[i]][["updatedelta"]]

    # Run the MCMC chain
    chainresult = mainMCMC(iters, genenum, selectgene, lastvaluelist, newallhash,
                           finaldata, finalinfo, statdiff, allchooserow, allstartrow, newaijindex,
                           updateaijmat, updatebeta, updatealpha, updatelambda, updateka, updatedelta,
                           betasd, alphasd, lambdasd, kasd, maxdelta, eachdelta)

    # Extract results from the chain
    chain = chainresult[[1]]
    prob = chainresult[[3]]
    updateaijmat = matrix(chain[iters, 1:(genenum^2)], genenum, genenum)
    updatedelta = matrix(chain[iters, (genenum^2 + 1):(2 * genenum^2)], genenum, genenum)
    updatebeta = chain[iters, 2 * genenum^2 + 1]
    updatealpha = chain[iters, 2 * genenum^2 + 2]
    updatelambda = chain[iters, 2 * genenum^2 + 3]
    updateka = chain[iters, 2 * genenum^2 + 4]

    # Return the results
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

  # Run subsequent MCMC chains until the maximum number of iterations is reached
  while(alliters < num_maxiters){
    alliters = alliters + iters
    eachlist = sfLapply(1:mcmcnum, eachsinglemcmc, p, allpara, genenum, lastvaluelist, iters, selectgene, newallhash,
                        finaldata, finalinfo, statdiff, allchooserow, allstartrow, newaijindex,
                        betasd, alphasd, lambdasd, kasd, maxdelta, eachdelta)
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
      subchain[[i]] = sparseMatrix(i = 1, j = 1, x = 1, dims = c(iters, 2 * genenum^2 + 4))
      subchain[[i]][, 1:(genenum^2)] = chain[, 1:(genenum^2)]
      subchain[[i]][, (genenum^2 + 1):(2 * genenum^2)] = chain[, (genenum^2 + 1):(2 * genenum^2)]
      subchain[[i]][, 2 * genenum^2 + 1] = chain[, 2 * genenum^2 + 1]
      subchain[[i]][, 2 * genenum^2 + 2] = chain[, 2 * genenum^2 + 2]
      subchain[[i]][, 2 * genenum^2 + 3] = chain[, 2 * genenum^2 + 3]
      subchain[[i]][, 2 * genenum^2 + 4] = chain[, 2 * genenum^2 + 4]
      allchain[[i]] = rbind(allchain[[i]], subchain[[i]])
      allprob[[i]] = c(allprob[[i]], thisprob)
    }

    # Save the results of the current set of MCMC chains
    save(allchain, file = paste(realresultpath, subtree, "_", p, "_result", ".RData", sep = ""))
    save(allprob, file = paste(realresultpath, subtree, "_", p, "_allprob", ".RData", sep = ""))
    choosepara = maxEstimate(allchain, allprob, genenum)
    preaij = choosepara[[1]]
    predelta = choosepara[[2]]
    preparas = choosepara[[3]]
    write.csv(preaij, file = paste(realresultpath, subtree, "_", p, "_preaij.csv", sep = ""), row.names = FALSE)
    write.csv(predelta, file = paste(realresultpath, subtree, "_", p, "_predelta.csv", sep = ""), row.names = FALSE)
    write.csv(preparas, file = paste(realresultpath, subtree, "_", p, "_preparas.csv", sep = ""), row.names = FALSE)
  }

  # Stop the parallel processing
  sfStop()
}
