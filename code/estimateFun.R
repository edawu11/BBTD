# (1) catchTable ---------------------------------------------------------
# Define the function catchTable to count occurrences of a specific numeric value in a vector
# Inputs:
#   a: Numeric vector, the data in which occurrences are to be counted
#   num: Numeric, the specific number whose occurrences are to be counted
# Output:
#   Numeric, the count of occurrences of 'num' in vector 'a'
catchTable = function(a, num){
  # Convert the input vector to numeric type to ensure proper handling
  a = as.numeric(a)
  
  # Convert the numeric vector into a factor, specifying possible levels to capture all potential values even if absent
  a = factor(a, levels = c(-1, 0, 1), labels = c(-1, 0, 1))
  
  # Generate a frequency table of the factor levels
  temp = table(a)
  
  # Identify the index of the desired number 'num' within the names of the table, ensuring it matches as character
  namenum = which(names(temp) == as.character(num))
  
  # Return the count of occurrences of 'num'; if 'num' is not present, handle potential NA by converting to numeric
  return(as.numeric(temp)[namenum])
}


# (2) maxEstimate --------------------------------------------------------
# Define the function maxEstimate to find parameter estimates with the highest posterior probabilities from MCMC chains
# Inputs:
#   allchain: List of matrices, each containing parameter chains for different MCMC runs
#   allprob: List of numeric vectors, each containing posterior probabilities corresponding to the parameter chains
#   genenum: Integer, number of genes, used to reshape the parameter vectors into matrices
# Output:
#   List containing matrices of aij, delta, and other parameter estimates (beta, alpha, lambda, ka)
maxEstimate = function(allchain, allprob, genenum){
  # Calculate the number of MCMC runs and the length of each chain
  mcmcnum = length(allprob)
  eachlen = nrow(allchain[[1]])

  # Discard the first half of each chain to allow for burn-in
  for(i in 1:mcmcnum){
    allchain[[i]] = allchain[[i]][(eachlen/2 + 1):eachlen, , drop = FALSE]
    allprob[[i]] = allprob[[i]][(eachlen/2 + 1):eachlen]
  }

  # Initialize a vector to store the maximum probabilities from each chain
  listvalue = c()
  for(i in 1:mcmcnum){
    listvalue = c(listvalue, max(allprob[[i]]))
  }
  
  # Identify the index of the chain with the highest overall probability
  maxlist_index = which.max(listvalue)
  chooseprob = allprob[[maxlist_index]]
  
  # Find the index of the maximum probability within the selected chain
  maxprob_index = which.max(chooseprob)
  choosevalue = allchain[[maxlist_index]][maxprob_index,]

  # Reshape the chosen parameters back into matrices and scalars for interpretation
  chooseaijmat = matrix(choosevalue[1:(genenum^2)], genenum, genenum)
  choosedelta = matrix(choosevalue[((genenum^2) + 1):(2 * genenum^2)], genenum, genenum)
  choosebeta = choosevalue[2 * genenum^2 + 1]
  choosealpha = choosevalue[2 * genenum^2 + 2]
  chooselambda = choosevalue[2 * genenum^2 + 3]
  chooseka = choosevalue[2 * genenum^2 + 4]

  # Organize the scalar parameters into a data frame for easy viewing
  finalvalue = data.frame(beta = choosebeta,
                          alpha = choosealpha,
                          lambda = chooselambda,
                          ka = chooseka)

  # Package the matrices and final parameter values into a list to return
  result = list()
  result[[1]] = chooseaijmat
  result[[2]] = choosedelta
  result[[3]] = finalvalue
  return(result)
}

# (3) evaAijDelta ----------------------------------------------------------------
# Define the function evaAijDelta to evaluate the agreement between two gene interaction matrices and their delta values
# Inputs:
#   a_gene: Matrix, ground truth gene interaction matrix
#   b_gene: Matrix, predicted gene interaction matrix
#   a_delta: Matrix, delta values associated with the first gene interaction matrix
#   b_delta: Matrix, delta values associated with the second gene interaction matrix
#   genenum: Integer, number of genes, used to iterate over the matrices
# Output:
#   Data frame containing evaluation metrics such as TPR, PPR, TNR, PNR, and ACC
evaAijDelta = function(a_gene, b_gene, a_delta, b_delta, genenum){
  # Initialize a list to store the evaluation metrics
  pingjia = list()
  
  # Calculate TPR: True Positive Rate
  choosenum = 0
  allnum_cujin = catchTable(a_gene, 1)  # Count all '1' values in a_gene
  if(allnum_cujin == 0){
    pingjia$TPR = 99  # Assign 99 if there are no positive values to avoid division by zero
  } else {
    for(i in 1:genenum){
      for(j in 1:genenum){
        if((a_gene[i, j] == 1) & (b_gene[i, j] == 1)){
          choosenum = choosenum + 1
        }
      }
    }
    pingjia$TPR = choosenum / allnum_cujin
  }
  
  # Calculate PPR: Positive Predictive Rate
  choosenum = 0
  allnum_cujin = catchTable(b_gene, 1)  # Count all '1' values in b_gene
  if(allnum_cujin == 0){
    pingjia$PPR = 99
  } else {
    for(i in 1:genenum){
      for(j in 1:genenum){
        if((a_gene[i, j] == 1) & (b_gene[i, j] == 1)){
          choosenum = choosenum + 1
        }
      }
    }
    pingjia$PPR = choosenum / allnum_cujin
  }
  
  # Calculate TNR: True Negative Rate
  choosenum = 0
  allnum_yizhi = catchTable(a_gene, -1)  # Count all '-1' values in a_gene
  if(allnum_yizhi == 0){
    pingjia$TNR = 99
  } else {
    for(i in 1:genenum){
      for(j in 1:genenum){
        if((a_gene[i, j] == -1) & (b_gene[i, j] == -1)){
          choosenum = choosenum + 1
        }
      }
    }
    pingjia$TNR = choosenum / allnum_yizhi
  }
  
  # Calculate PNR: Positive Negative Rate
  choosenum = 0
  allnum_yizhi = catchTable(b_gene, -1)  # Count all '-1' values in b_gene
  if(allnum_yizhi == 0){
    pingjia$PNR = 99
  } else {
    for(i in 1:genenum){
      for(j in 1:genenum){
        if((a_gene[i, j] == -1) & (b_gene[i, j] == -1)){
          choosenum = choosenum + 1
        }
      }
    }
    pingjia$PNR = choosenum / allnum_yizhi
  }
  
  # Calculate ACC: Accuracy for delta values
  a_delta = round(a_delta, 3)
  b_delta = round(b_delta, 3)
  choosenum = 0
  allnum_nozero = catchTable(b_gene, 1) + catchTable(b_gene, -1)  # Count all non-zero values in b_gene
  if(allnum_nozero == 0){
    pingjia$ACC = 1
  } else {
    for(i in 1:genenum){
      for(j in setdiff(1:genenum, i)){  # Exclude diagonal elements (self-interactions)
        if(a_delta[i, j] == b_delta[i, j] & b_gene[i, j] != 0){
          choosenum = choosenum + 1
        }
      }
    }
    pingjia$ACC = choosenum / allnum_nozero
  }
  
  # Convert the list of metrics to a data frame for easier viewing and further analysis
  pingjia = as.data.frame(pingjia)
  return(pingjia)
}