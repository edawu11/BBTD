# (1) preScreen -----------------------------------------------------------
# Define the function preScreen to pre-screen aij based on statistical p-values
# Inputs:
#   aijindexlist: List, indices of potential interactions for each gene
#   selectgene: Vector, list of selected gene identifiers
#   GoalDataset: DataFrame, contains target data for analysis
#   finalinfo: List, contains final processed data from previous analysis steps
#   allchooserow: Vector, indices of rows chosen for analysis
#   realaijmat: Matrix, actual gene interaction data for comparison
#   maxdelta_time: Numeric, max delta time
#   Mc: Numeric, number of data points for interpolation
#   yuzhipoint: Numeric, threshold for filtering based on p-values
#   ifallpvalue: Logical, indicates whether to calculate or load p-values
#   ifreal: Logical, indicates whether to compare with real data matrix
#   simupath: String, directory path for saving/loading results
#   p: Numeric, a parameter that might be used to distinguish different runs or settings
#   bosscell: String, identifier for the main cell or dataset being analyzed
#   thishash: Hash or List, previously calculated hash tables or lists containing data
# Output:
#   None directly; function saves results of the screening or loads previous screening data
preScreen = function(aijindexlist, selectgene, GoalDataset, finalinfo, allchooserow, realaijmat,
                                       maxdelta_time, Mc, yuzhipoint, ifallpvalue, ifreal, simupath, p, bosscell, thishash){
  
  # Calculate delta values based on Mc and maxdelta_time
  eachdelta = 1 / (round(Mc) - 1)
  maxdelta = maxdelta_time * eachdelta
  genenum = length(selectgene)
  
  # Decide whether to calculate all p-values or load them from file
  if(ifallpvalue == TRUE){
    all_pvalue = aijPvalue(aijindexlist, selectgene, genenum, GoalDataset, finalinfo, allchooserow,
                                     maxdelta, eachdelta, thishash)
    save(all_pvalue, file = paste(simupath, "mer_", bosscell, "_all_pvalue_", p, "_", "_.RData", sep = ""))
  } else {
    load(paste(simupath, "mer_", bosscell, "_all_pvalue_", p, "_", "_.RData", sep = ""))
  }
  
  # Filter the aij indices based on the p-value threshold
  newaijindex = aijindexlist
  for(i in 1:length(all_pvalue)){
    thepvalue = all_pvalue[[i]]
    tempindex = which(thepvalue < yuzhipoint)
    newaijindex[[i]] = aijindexlist[[i]][tempindex]
  }
  
  # If not comparing with real data, calculate initial results based on the new indices
  if(ifreal == FALSE){
    initialresult = list()
    initialresult[[1]] = length(unlist(newaijindex))
    initialmat = matrix(0, nrow = genenum, ncol = genenum)
    for(i in 1:length(newaijindex)){
      initialmat[i, newaijindex[[i]]] = 1
    }
    
    # Compare predicted interaction matrix with the real interaction matrix
    colnowgene = factor(as.numeric(initialmat), levels = c(0, 1), labels = c(0, 1))
    conmat = table(Prediction = factor(as.numeric(initialmat), levels = c(0, 1), labels = c(0, 1)),
                   Reference = factor(as.numeric(as.matrix(abs(realaijmat[, 1:genenum]))), levels = c(0, 1), labels = c(0, 1)))
    conmat[1, 1] = conmat[1, 1] - genenum
    
    # Calculate false positive rate
    if(sum(conmat[, 2]) == 0){
      FPR = 0
    } else {
      FPR = conmat[1, 2] / sum(conmat[, 2])
    }
    
    # Compile confusion matrix and other results
    conf_matrix = confusionMatrix(conmat)
    initialresult[[2]] = conf_matrix
    initialresult[[3]] = genenum * (genenum - 1) - initialresult[[1]]
    initialresult[[4]] = FPR
    initialresult[[5]] = conmat[, 2]
    
    # Save initial results
    save(initialresult, file = paste(simupath, "mer_", bosscell, "_initial_result_", p, "_.RData", sep = ""))
  }
  
  # Save the new aij index list
  save(newaijindex, file = paste(simupath, "mer_", bosscell, "_newaijindex_", p, "_.RData", sep = ""))
}

# (2) aijPvalue -----------------------------------------------------------
# Define the function aijPvalue to calculate p-values for gene interaction indices
# Inputs:
#   aijindexlist: List, indices of potential interactions for each gene
#   selectgene: Vector, list of selected gene identifiers
#   genenum: Integer, number of genes
#   GoalDataset: DataFrame, contains target data for analysis
#   finalinfo: DataFrame, contains final processed data from previous analysis steps
#   allchooserow: Vector, indices of rows chosen for analysis
#   maxdelta: Numeric, the maximum possible value for time delay in regulation
#   eachdelta: Numeric, the increment between each possible delta value within the range from 0 to maxdelta
#   thishash: Hash or List, previously calculated hash tables or lists containing data
# Output:
#   List of vectors, each containing p-values for potential interactions of each gene
aijPvalue = function(aijindexlist, selectgene, genenum, GoalDataset, finalinfo, allchooserow,
                     maxdelta, eachdelta, thishash){
  # Retrieve cell and time information from finalinfo
  oldcell = finalinfo$cell
  oldtime = finalinfo$time
  # Generate a sequence of delta values
  alldelta = seq(0, maxdelta, by = eachdelta)
  
  # Define a sub-function to calculate p-values for one gene against all others
  subPvalue = function(numindex, aijindexlist, selectgene, genenum, allchooserow, oldcell, oldtime, alldelta, thishash){
    thisgenename = selectgene[numindex]
    
    # Retrieve future gene values from the GoalDataset and replace NA with 0
    colfuturegene = GoalDataset[[numindex]][allchooserow]
    colfuturegene = replace_na(colfuturegene, 0)
    
    # Define another sub-function to calculate the minimum p-value for interactions of a specific gene
    eachMin = function(eachindex, alldelta, thisgenename, eachdelta, colfuturegene,
                       finalinfo, allchooserow, oldcell, oldtime, thishash){
      origenename = selectgene[eachindex]
      # Calculate values for each delta delay
      origenenowvalue = sapply(alldelta, orivalueFun, origenename, thisgenename, eachdelta,
                               finalinfo, allchooserow, oldcell, oldtime, thishash)
      # Compute p-values for the past values compared to actual values
      thisgenepvalue = apply(origenenowvalue, 2, inGeneCol, colfuturegene)
      return(min(thisgenepvalue))
    }
    
    # Calculate p-values for all potential interactions of this gene
    thisgeneallpvalue = sapply(aijindexlist[[numindex]], eachMin, alldelta, thisgenename, eachdelta, colfuturegene,
                               finalinfo, allchooserow, oldcell, oldtime, thishash)
    
    return(thisgeneallpvalue)
  }
  
  # Apply the subPvalue function to each gene to compute all p-values
  allpvalue = sapply(1:genenum, subPvalue, aijindexlist, selectgene, genenum,
                     allchooserow, oldcell, oldtime, alldelta, thishash, simplify = FALSE)
  names(allpvalue) = selectgene
  
  return(allpvalue)
}

# (3) orivalueFun -------------------------------------------------------------
# Define the function orivalueFun to compute original gene values based on specified time delay.
# Inputs:
#   choosedelta: Numeric, the time delta for backtracking the original values
#   origenename: String, the name of the originating gene for which values are being computed
#   futgenename: String, the name of the future gene (not used in this function, may be for context)
#   eachdelta: Numeric, each delta step size (not directly used in this function)
#   finalinfo: DataFrame, contains 'lastrow' and timing information
#   allchooserow: Vector, indices of the chosen rows for analysis
#   oldcell: Vector, cell identifiers corresponding to each row in finalinfo
#   oldtime: Vector, timing information corresponding to each row in finalinfo
#   thishash: List of hash tables, contains historical data for each gene indexed by cell
# Output:
#   Numeric vector, contains the predicted past values of the gene based on choosedelta
orivalueFun = function(choosedelta, origenename, futgenename, eachdelta,
                        finalinfo, allchooserow, oldcell, oldtime, thishash){
  
  # Initialize a vector to hold the predicted values for each choice row
  alllastvalue = rep(0, length(allchooserow))
  lastrow = finalinfo$lastrow
  
  # Define a sub-function to calculate the value for a single row based on the specified delta
  singleGeneOrivalue = function(k, allchooserow, choosedelta, origenename, lastrow){
    oldrow = allchooserow[k]
    thistime = oldtime[oldrow]
    thiscell = oldcell[oldrow]
    ini_point = floor(thistime)  # Initial time point for the current observation
    lasttime = round(thistime - choosedelta, 3)  # Calculate the time point to look up
    
    # Determine if the time point goes before the start of the observed time
    if(lasttime < ini_point){
      lastcell = str_sub(thiscell, 1, -2)  # Assume previous cell if crossing time boundary
      thish = thishash[[origenename]][[lastcell]]
      # Check hash for the previous value; default to 0 if not found
      if(is.null(thish) || length(intersect(keys(thish), as.character(lasttime))) == 0){
        lastvalue = 0
      } else {
        lastvalue = values(thish, round(lasttime, 3))
      }
      
    } else {
      thish = thishash[[origenename]][[thiscell]]
      if(is.null(thish) || length(intersect(keys(thish), as.character(lasttime))) == 0){
        lastvalue = 0
      } else {
        lastvalue = values(thish, round(lasttime, 3))
      }
    }
    return(as.numeric(lastvalue))
  }
  
  # Apply the sub-function to each row index, calculating past values
  alllastvalue = sapply(1:length(allchooserow), singleGeneOrivalue, allchooserow,
                        choosedelta, origenename, lastrow)
  
  return(alllastvalue)
}

# (4) inGeneCol -----------------------------------------------------------
# Define the function inGeneCol to calculate a p-value using Fisher's Exact Test for gene comparison
# Inputs:
#   colnowgene: Numeric vector, current gene states (0 or 1)
#   colfuturegene: Numeric vector, future gene states (0 or 1), to be compared against current states
# Output:
#   Numeric, the p-value from Fisher's Exact Test indicating the statistical significance of the association between states
inGeneCol = function(colnowgene, colfuturegene){
  # Convert numeric gene state vectors into factors with levels specified to ensure all cases are covered
  colnowgene = factor(colnowgene, levels = c(0, 1), labels = c(0, 1))
  colfuturegene = factor(colfuturegene, levels = c(0, 1), labels = c(0, 1))

  # Create a contingency table of current and future gene states
  ktab = table(colnowgene, colfuturegene)
  
  # Reconvert factors to numeric for Fisher test compatibility
  colnowgene = as.numeric(as.character(colnowgene))
  colfuturegene = as.numeric(as.character(colfuturegene))
  
  # Check if any elements in the contingency table are zero and adjust by adding missing levels if necessary
  if(any(ktab == 0)){
    if(ktab[1] == 0){
      colnowgene = c(colnowgene, 0)  # Add false positive
      colfuturegene = c(colfuturegene, 0)
    }
    if(ktab[2] == 0){
      colnowgene = c(colnowgene, 1)  # Add true negative
      colfuturegene = c(colfuturegene, 0)
    }
    if(ktab[3] == 0){
      colnowgene = c(colnowgene, 0)  # Add false negative
      colfuturegene = c(colfuturegene, 1)
    }
    if(ktab[4] == 0){
      colnowgene = c(colnowgene, 1)  # Add true positive
      colfuturegene = c(colfuturegene, 1)
    }
    
    # Conduct Fisher's Exact Test and return the p-value
    result = fisher.test(matrix(c(colnowgene, colfuturegene), ncol = 2))
    return(result$p.value)
    
  } else {
    # Conduct Fisher's Exact Test and return the p-value when no adjustment is necessary
    result = fisher.test(matrix(c(colnowgene, colfuturegene), ncol = 2))
    return(result$p.value)
  }
}
