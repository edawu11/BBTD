# (1) singleAijPrior --------------------------------------------------------
# Define the function singleAijPrior to calculate log-prior probabilities for a_ij
# Inputs:
#   genenum: Integer, the total number of genes in the network
#   lambda: Numeric, hyper-parameter
#   ka: Numeric, hyper-parameter
# Output:
#   priorhash: A hash table where each key is a number of possible connections (0 to m)
#              and each value is the log-prior probability of that number of connections
singleAijPrior = function(genenum, lambda, ka){
  # Calculate m, the maximum possible connections a gene can have (total genes - 1)
  m = genenum - 1
  
  # Create a sequence from 1 to m, representing the possible number of connections
  z = 1:m
  
  # Compute the normalization constant
  cc = sum(exp(-lambda * z))
  
  # Define a nested function to compute log-prior for a row of a_ij
  row_logprior = function(x, m, lambda, ka, cc){
    if(x == 0){
      # Calculate log probability of having zero connections
      log_p = log(1 - ka)
    }
    else{
      # Calculate log probability of having x connections using a detailed formula
      log_p = log(ka) + (-lambda * x) - log(cc) - log(choose(m, x)) - x * log(2)
    }
    return(log_p)
  }
  
  # Apply the row_logprior function to each value from 0 to m
  allprior = sapply(0:m, row_logprior, m, lambda, ka, cc)
  
  # Create a hash table to store each number of connections and its corresponding log-prior probability
  priorhash = hash(0:m, allprior)
  
  # Return the hash table containing log-prior probabilities
  return(priorhash)
}

# (2) invAijPrior ---------------------------------------------------------
# Define the function invAijPrior to compute the cumulative log-prior probabilities for a matrix of gene connections
# Inputs:
#   aijmat: Numeric matrix, where each row represents the connection status from one gene to all other genes
#   alllogprior: A hash or list where keys are the number of non-zero connections and values are the corresponding log-prior probabilities
# Output:
#   The sum of all log-prior probabilities for the connections described in aijmat
invAijPrior = function(aijmat, alllogprior){
  # Define a nested function to compute the log-prior for a single row of connection statuses
  rowmatrix_fun = function(rowmatrix, alllogprior){
    # Count the number of non-zero entries in the row (indicative of existing connections) and convert to character
    nozeronum = as.character(sum(abs(rowmatrix)))
    # Retrieve the log-prior probability for this number of connections from alllogprior
    return(alllogprior[[nozeronum]])
  }
  
  # Apply the nested function to each row of aijmat, resulting in a vector of log-prior probabilities
  temp = apply(aijmat, 1, rowmatrix_fun, alllogprior)
  
  # Sum up the log-prior probabilities from all rows to get the total log-prior for the matrix
  return(sum(temp))
}

# (3) invDeltaPrior -------------------------------------------------------
# Define the function invDeltaPrior to compute the cumulative log-prior probabilities for delta updates in a gene regulatory network
# Inputs:
#   updatedelta: Numeric matrix, representing time delay
#   updateaijmat: Numeric matrix, representing a_ij
#   genenum: Integer, the total number of genes in the network
# Output:
#   The sum of log-prior probabilities of updatedelta
invDeltaPrior = function(updatedelta, updateaijmat, genenum){
  # Convert the update matrices to numeric vectors for easier indexing and manipulation
  allaij = as.numeric(updateaijmat)
  alldelta = as.numeric(updatedelta)
  
  # Find indices where there are connections (non-zero entries in the adjacency matrix) and calculate their log-prior probabilities
  nozeroindex = which(allaij != 0)
  nozeroprob = length(nozeroindex) * log(1 / 6)  # Assuming uniform probability for the existence of connections
  
  # Find indices where there are no connections and exclude self-connections (diagonal entries in the matrix)
  allzeroindex = which(allaij == 0)
  selfindex = seq(1, genenum^2, by = genenum + 1)  # Indices of diagonal entries (self-connections)
  yeszeroindex = setdiff(allzeroindex, selfindex)  # Non-self zero entries
  
  # Calculate probabilities for zeros in delta matrix where the corresponding entries in the adjacency matrix are zero
  zeronum = length(which(alldelta[yeszeroindex] == 0.0))
  if(zeronum == length(yeszeroindex)){
    # If all zeros in delta correspond to zeros in adjacency, apply a high probability
    zeroprob = zeronum * log(0.999)
  }
  else{
    # Otherwise, apply different probabilities for zero and non-zero deltas
    temp1 = zeronum * log(0.999)
    temp2 = (length(yeszeroindex) - zeronum) * log(0.0002)
    zeroprob = temp1 + temp2
  }
  
  # Return the sum of log-prior probabilities for non-zero and zero updates
  return(nozeroprob + zeroprob)
}

# (4) geneAijRow -------------------------------------------------------
# Define the function geneAijRow to simulate a row of gene regulatory connections based on given probabilities
# Inputs:
#   aijindex: Numeric vector, indices where connections can potentially exist
#   lambda: Numeric, hyper-parameter
#   ka: Numeric, hyper-parameter
#   genenum: Integer, the total number of genes in the network
# Output:
#   bengenrow: Numeric vector representing a single row of an adjacency matrix for gene regulatory interactions
geneAijRow = function(aijindex, lambda, ka, genenum){
  # Initialize the output vector for the gene row with all zeros (no connections)
  bengenrow = rep(0, genenum)
  
  # Decide randomly based on probability ka if this gene will have any connections at all
  ifzero = rbinom(1, 1, prob = ka)
  
  if(ifzero == 0){
    # If the gene has no connections, return the vector filled with zeros
    return(bengenrow)
  }
  else{
    # If the gene does have connections, proceed with determining the specifics of those connections
    
    # Determine the number of potential connection positions based on aijindex
    m = length(aijindex)
    z = 1:m
    
    # Calculate the sum of exponential decays for normalization
    sumexp = sum(exp(-lambda * z))
    
    # Calculate the probability of connecting to each possible position using the exponential decay
    p = exp(-lambda * z) / sumexp
    
    # Randomly choose a number of positions to connect based on probabilities p
    choosenum = sample(z, 1, replace = FALSE, prob = p)
    
    # Randomly select the actual positions for connections from aijindex
    choosepos = sample(aijindex, choosenum, replace = FALSE)
    
    # Set the connections in the output vector; randomly decide on positive or negative influence for each connection
    bengenrow[choosepos] = sample(c(-1, 1), choosenum, replace = TRUE, prob = c(0.5, 0.5))
    
    # Return the completed row vector representing the gene's regulatory influences
    return(bengenrow)
  }
}

# (5) geneAij ----------------------------------------------------------
# Define the function geneAij to generate an adjacency matrix for gene regulatory networks
# Inputs:
#   genenum: Integer, the total number of genes in the network
#   lambda: Numeric, hyper-parameter
#   ka: Numeric, hyper-parameter
#   aijindexlist: List of numeric vectors, each vector containing indices where connections can potentially exist for each gene
# Output:
#   bengen: Matrix, the adjacency matrix for gene regulatory interactions, with gene indices appended as the last column
geneAij = function(genenum, lambda, ka, aijindexlist){
  # Apply the geneAijRow function to each set of indices in aijindexlist to simulate the connections for each gene
  bengenlist = lapply(aijindexlist, geneAijRow, lambda, ka, genenum)
  
  # Unlist the list of vectors into a single vector
  bengen = unlist(bengenlist)
  
  # Reshape the unlisted vector into a matrix with one row per gene
  bengen = matrix(bengen, genenum, genenum, byrow = TRUE)
  
  # Append a column to the matrix that contains the indices of the genes
  bengen = cbind(bengen, 1:genenum)
  
  # Return the completed adjacency matrix with gene indices
  return(bengen)
}

# (6) subDeltaPrior -------------------------------------------------------
# Define the function subDeltaPrior to sample delta values for updates based on connection presence
# Inputs:
#   subaij: Integer, indicating whether a regulatory connection exists (0 for no connection, non-zero for a connection)
#   maxdelta: Numeric, the maximum possible value for time delay in regulation
#   eachdelta: Numeric, the increment between each possible delta value within the range from 0 to maxdelta
# Output:
#   delta: Numeric, the sampled delta value corresponding to the aij
subDeltaPrior = function(subaij, maxdelta, eachdelta){
  if(subaij == 0){
    # If there is no existing connection, sample a delta value from a sequence, favoring a delta of zero
    delta = sample(seq(0, maxdelta, by = eachdelta), 1, replace = FALSE, 
                   prob = c(0.999, rep(0.0002, 5)))
  }
  else{
    # If there is an existing connection, sample a delta value uniformly from the possible values
    delta = sample(seq(0, maxdelta, by = eachdelta), 1, replace = FALSE, 
                   prob = rep(1 / 6, 6))
  }
  # Return the sampled delta value
  return(delta)
}


# (7) geneDelta -----------------------------------------------------------
# Define the function geneDelta to generate a matrix of delta values for updates in a gene regulatory network
# Inputs:
#   updateaijmat: Numeric matrix, representing the adjacency matrix of gene connections (1 for a connection, 0 for no connection)
#   maxdelta: Numeric, the maximum possible value for time delay in regulation
#   eachdelta: Numeric, the increment between each possible delta value within the range from 0 to maxdelta
#   genenum: Integer, the total number of genes in the network
# Output:
#   Matrix, a matrix of delta values
geneDelta = function(updateaijmat, maxdelta, eachdelta, genenum){
  # Convert the adjacency matrix to a numeric vector for easier manipulation
  aijnum = as.numeric(updateaijmat)
  
  # Apply the subDeltaPrior function to each element in the adjacency matrix
  # This function samples delta values based on whether a regulatory connection exists
  alldelta = sapply(aijnum, subDeltaPrior, maxdelta, eachdelta)
  
  # Reshape the vector of delta values back into a matrix with the same dimensions as the adjacency matrix
  return(matrix(alldelta, genenum, genenum))
}