# (1) expressData ---------------------------------------------------------
# Define the function expressData to process and export experimental data based on specified criteria
# Inputs:
#   fileindex: Integer, index to select specific files from filename lists
#   oldfilename: Character vector, list of old filenames for reading initial data
#   hufilename: Character vector, list of filenames containing onset information
#   newfilename: Character vector, list of new filenames for saving the processed data
#   hufilepath: String, path to the directory containing onset files
#   oldfilepath: String, path to the directory containing old data files
#   expdatapath: String, path to the directory where processed data will be saved
# Output:
#   None directly; function saves processed data files to specified path
expressData = function(fileindex, oldfilename, hufilename, newfilename,
                        hufilepath, oldfilepath, expdatapath){
  # Read an Excel file containing onset information, specifically from the third sheet
  hufile = read_excel(path = paste(hufilepath, hufilename[fileindex], sep = ""), sheet = 3)
  
  # Process 'cell' column by removing any quotation marks
  choosecell = str_replace_all(hufile$cell, "\"", "")
  
  # Process 'time' column by converting to numeric after removing any quotation marks
  choosetime = as.numeric(str_replace_all(hufile$time, "\"", "")) 
  
  # Read a CSV file containing old experimental data
  oldfile = read_csv(file = paste(oldfilepath, oldfilename[fileindex], sep = ""), show_col_types = FALSE)
  
  # Select relevant columns from the old data
  oldfile = oldfile %>% select(cell, time, blot)
  
  # Get unique values from the 'cell' column
  oldcell = unique(oldfile$cell)
  
  # Filter the old data based on the first selected 'cell' from the onset file
  choosedata = oldfile %>% filter(cell %in% str_subset(oldcell, paste("^", choosecell[1], sep = "")))
  
  # Loop through the remaining cells from the onset file to filter and combine data
  for(k in 2:length(choosecell)){
    tempdata = oldfile %>% filter(cell %in% str_subset(oldcell, paste("^", choosecell[k], sep = "")))
    choosedata = choosedata %>% bind_rows(tempdata)
  }
  
  # Add a column to indicate presence of 'blot' data initially setting all to 1
  choosedata = choosedata %>% mutate(ifblot = 1)
  
  # Update 'ifblot' values based on comparison of 'time' in old data to 'time' in onset file
  for(k in 1:length(choosecell)){
    subdata = choosedata %>% filter(cell == choosecell[k]) %>% 
      mutate(ifblot = ifelse(time >= choosetime[k], 1, 0))
    chooseindex = which(choosedata$cell == choosecell[k])
    choosedata[chooseindex,] = subdata
  } 
  
  # Write the processed data to a new CSV file in the specified export path
  write_csv(choosedata, file = paste(expdatapath, newfilename[fileindex], sep = ""))
}

# (2) cellnumStat ---------------------------------------------------------
# Define the function cellnumStat to aggregate cell data based on specified subtree criteria and save the results
# Inputs:
#   subtree: Character, the specific subtree identifier (e.g., "E") for which data is to be processed
#   usepath: String, the directory path where the input CSV files are stored
#   savepath: String, the directory path where processed files will be saved
# Output:
#   None directly; function writes processed data files to the specified save path
cellnumStat = function(subtree, usepath, savepath){
  
  # Adjust the subtree pattern based on specific criteria
  if(subtree == "E"){
    usesubtree = paste("E", "[^M]|E", sep="")
  }
  else{
    usesubtree = subtree
  }
  
  # List all CSV files in the specified directory
  allgenefile = list.files(path = usepath, pattern = ".csv")
  
  # Initialize an empty list to store data
  alldata = list()
  
  # Extract the gene name from the first file name by removing the last 4 characters (.csv)
  genename = str_sub(allgenefile[1], start = 1, end = nchar(allgenefile[1])-4)
  
  # Read the first CSV file and filter data based on the subtree pattern
  genedata = read_csv(file = paste(usepath, allgenefile[1], sep = ""), show_col_types = FALSE)
  allchoosedata = genedata %>% filter(cell %in% str_subset(genedata$cell, paste("^", usesubtree, sep = ""))) %>% 
    group_by(cell) %>% summarise(cellname = sum(ifblot))
  names(allchoosedata) = c("cell", genename)
  
  # Process remaining files in the same way, join each file's data with the aggregated data frame
  for(i in allgenefile[-1]){
    genedata = read_csv(file = paste(usepath, i, sep = ""), show_col_types = FALSE)
    genename = str_sub(i, start = 1, end = nchar(i)-4)
    choosedata = genedata %>% filter(cell %in% str_subset(genedata$cell, paste("^", usesubtree, sep = ""))) %>% 
      group_by(cell) %>% summarise(cellname = sum(ifblot))
    names(choosedata) = c("cell", genename)
    allchoosedata  = allchoosedata %>% full_join(choosedata, by = "cell")
  }
  
  # Write the fully joined data to a CSV file
  write_csv(allchoosedata, file = paste(savepath, "all_", subtree, "_cell.csv", sep = ""))
  
  # Filter out columns with all NA values and write the result to another CSV file
  allchoosedata = allchoosedata %>% select_if(~!all(is.na(.)))
  write_csv(allchoosedata, file = paste(savepath, "all_", subtree, "_nona_cell.csv", sep = ""))
}

# (3) geneCopyFun ---------------------------------------------------------
# Define the function geneCopyFun to process and save gene data based on cell lineage information
# Inputs:
#   subtree: Character, the specific subtree identifier for which data is to be processed
#   usedata: DataFrame, contains data to be used for generating lineage information
#   allfilename: List of strings, contains filenames which include identifiers and gene names
#   savepath: String, the directory path where the R data file will be saved
# Output:
#   None directly; function saves a comprehensive summary of gene data and lineage information as an R data file
geneCopyFun = function(subtree, usedata, allfilename, savepath){
  # Extract copyname and genename from filenames, removing the last 4 characters (.csv)
  copyname = str_sub(allfilename[[3]], 1, nchar(allfilename[[3]]) - 4)
  genename = allfilename[[4]]
  
  # Create a hash to map copyname to genename
  genecopy_h = hash(copyname, genename)
  
  # Sort the data by 'cell' column
  newdata = usedata %>% arrange(cell)
  
  # Extract the 'cell' column
  allcellname = newdata$cell
  
  # Initialize a vector to track the number of unique tree numbers, starting with 1 for the first cell
  alltreenum = rep(1, nrow(newdata))
  nowmother = allcellname[1]
  
  # Iterate through each cell to determine if it belongs to the same lineage (tree number remains the same) or a new one
  for(i in 2:nrow(newdata)){
    if(str_sub(allcellname[i], 1, nchar(nowmother)) == nowmother){
      alltreenum[i] = alltreenum[i-1]
    }
    else{
      nowmother = allcellname[i]
      alltreenum[i] = alltreenum[i-1] + 1
    }
  }
  
  # Create a hash to map each cell name to its corresponding tree number
  tree_h = hash(allcellname, alltreenum)
  
  # Initialize a list to store simulation information
  simu_info = list()
  simu_info[["genecopy"]] = genecopy_h
  simu_info[["copyname"]] = keys(genecopy_h)
  simu_info[["genename"]] = sort(as.character(unique(values(genecopy_h))))
  simu_info[["cellname"]] = allcellname
  simu_info[["Mc"]] = median(unlist(newdata[,-1]), na.rm = TRUE)  # Calculate the median of all non-cell-name data
  
  # Save the simulation information to an R data file
  save(simu_info, file = paste(savepath, subtree, "_simu_info.RData", sep = ""))
}

# (4) interData -----------------------------------------------------------
# Define the function interData to interpolate experimental data and save processed results
# Inputs:
#   cellname: Vector, list of cell identifiers
#   copyname: Vector, list of gene copy identifiers
#   Mc: Numeric, the number of data points used to interpolate
#   expdatapath: String, path to the directory containing experimental data files
#   subtree: String, identifier for the subtree used in naming saved files
#   maxdelta_time: Numeric, maximum delta time
#   savepath: String, directory path where processed files will be saved
# Output:
#   None directly; function saves interpolated and processed data to specified path
interData = function(cellname, copyname, Mc, expdatapath, subtree, maxdelta_time, savepath){
  # Calculate delta parameters based on Mc and maxdelta_time
  M = round(Mc)
  eachdelta = 1 / (M - 1)
  maxdelta = maxdelta_time * eachdelta
  
  # Initialize lists to hold functions and data for each copy
  allfun = list()
  copynum = length(copyname)
  length(allfun) = copynum
  names(allfun) = copyname
  
  # Initialize list to hold data for each cell
  allcelllist = list()
  length(allcelllist) = length(cellname)
  names(allcelllist) = cellname
  
  # Assign an empty list to each gene copy in allfun
  for(i in copyname){
    allfun[[i]] = allcelllist
  }

  # Read data for each gene copy, process it, and create spline functions for interpolation
  for(i in copyname){
    genefile = read.csv(paste(expdatapath, i, ".csv", sep = ""))
    for(j in cellname){
      choosedata = genefile %>% filter(cell == j)
      if(nrow(choosedata) == 0){
        next
      }
      if(nrow(choosedata) < 3 & nrow(choosedata) > 0){
        # Handle cases with insufficient data by appending data from the parent cell
        lastcell = str_sub(j, 1, nchar(j) - 1)
        lastcelldata = genefile %>% filter(cell == lastcell)
        alldata = lastcelldata %>% add_row(choosedata) %>% 
                  mutate(gap = nchar(cell)) %>%
                  mutate(scale_time = gap + (time - min(time)) / (max(time) - min(time) + 1)) %>%
                  mutate(newblot = blot - lag(blot)) %>% na.omit()
        allfun[[i]][[j]] = splinefun(x = alldata$scale_time, y = alldata$newblot)
      }
      else{
        # Regular case: create spline functions directly from the data
        choosedata = choosedata %>% mutate(gap = nchar(cell)) %>%
                    mutate(scale_time = gap + (time - min(time)) / (max(time) - min(time) + 1))  %>%
                    mutate(newblot = blot - lag(blot)) %>% na.omit()
        allfun[[i]][[j]] = splinefun(x = choosedata$scale_time, y = choosedata$newblot)
      }
    }
  }
  
  # Save all functions to an RData file
  save(allfun, file = paste(savepath, subtree, "_allfun.RData", sep = ""))
  
  # Create a new data frame to hold interpolated data
  newdata = matrix(0, nrow = length(cellname) * (M - 1), ncol = 2)
  newdata = data.frame(newdata)
  names(newdata) = c("cell", "time")
  newdata$cell = rep(cellname, each = M - 1)
  newdata = newdata %>% mutate(gap = nchar(cell)) %>% group_by(cell) %>% 
            mutate(time = seq(min(gap), min(gap) + 1, length.out = M)[-M]) %>% ungroup() %>% select(-gap)
  
  infodata = newdata

  # Define a nested function to interpolate data for each cell using precomputed spline functions
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
  
  # Define a nested function to interpolate data for all cells under a single gene copy
  singlecopyinter_fun = function(thiscopyname,cellname,infodata,allfun){
    thiscopyvalue = unlist(sapply(cellname,singlecellinter_fun,infodata,thiscopyname,allfun,simplify = F))
    return(thiscopyvalue)
  }

  # Apply the single copy interpolation function across all gene copies
  adddata = sapply(copyname,singlecopyinter_fun,cellname,infodata,allfun)
  # Combine the interpolated data back with the original data structure
  newdata = newdata %>% bind_cols(adddata)
  # Update column names to reflect cell, time, and each gene copy
  names(newdata) = c("cell", "time", copyname)

  # Calculate lineage tree indices based on cell names
  allcellname = newdata$cell

  alltreenum = rep(1, nrow(newdata))
  nowmother = allcellname[1]
  for(i in 2:nrow(newdata)){
    if(str_sub(allcellname[i], 1, nchar(nowmother)) == nowmother){
      alltreenum[i] = alltreenum[i-1]
    } else {
      nowmother = allcellname[i]
      alltreenum[i] = alltreenum[i-1] + 1
    }
  }

  # Create a hash to map cell names to their tree indices
  tree_h = hash(allcellname, alltreenum)
  
  # Assign tree indices to the data frame
  newdata$treeindex = 0
  for(i in 1:nrow(newdata)){
    newdata$treeindex[i] = tree_h[[newdata$cell[i]]]
  }

  # Arrange the data by cell and time, and compute additional lineage-based indices
  newdata = newdata %>% arrange(cell, time) %>%  
    mutate(originrow = 1:n(), lastrow = originrow - 1) %>% 
    group_by(treeindex) %>% mutate(ifchange = as.numeric(lag(cell) != cell)) %>%
    select(cell, time, originrow, lastrow, treeindex, ifchange, everything()) %>% ungroup() %>% 
    mutate(ifchange = replace_na(ifchange, 9))

  # Correct lastrow indices based on lineage changes
  changeindex = which(newdata$ifchange == 1)
  c_cell = newdata$cell
  lastrow = newdata$lastrow
  for(i in changeindex){
    thiscell = c_cell[i]
    lastcell = substr(thiscell, 1, nchar(thiscell) - 1)
    choosedata = newdata %>% filter(cell == lastcell) %>% select(originrow)
    temprow = choosedata$originrow
    lastrow[i] = temprow[length(temprow)]
  }
  newdata$lastrow = lastrow
  
  # Save the processed data with detailed cell lineage and interpolation information
  write_csv(newdata, file = paste(savepath, subtree, "_", maxdelta_time, "_condata.csv", sep = ""))
  
}

# (5) hashData ------------------------------------------------------------
# Define the function hashData to interpolate experimental data and store in hash tables for quick lookup
# Inputs:
#   allfun: List, contains spline functions for interpolation, indexed by copyname and cellname
#   realdata: List of DataFrames, contains real experimental data indexed by gene copy
#   copyname: Vector, list of gene copy identifiers
#   cellname: Vector, list of cell identifiers
#   maxdelta_time: Numeric, maximum delta time
#   Mc: Numeric, median count used to calculate delta
#   savepath: String, directory path where processed hash files will be saved
#   subtree: String, identifier for the subtree used in naming saved files
# Output:
#   None directly; function saves hash tables with interpolated data to specified path
hashData = function(allfun, realdata, copyname, cellname, maxdelta_time, Mc, savepath, subtree){
  # Calculate delta parameters based on Mc and maxdelta_time
  M = round(Mc)
  eachdelta = 1 / (M - 1)
  maxdelta = maxdelta_time * eachdelta

  # Initialize a list of hash tables for each gene copy and cell
  allhash = list()
  copynum = length(copyname)
  length(allhash) = copynum
  names(allhash) = copyname

  # Populate the hash list with empty hash tables for each cell under each copy
  for(i in copyname){
    for(j in cellname){
      allhash[[i]][[j]] = hash()
    }
  }
  
  # Generate a sequence of time points for interpolation
  allseq = seq(0, maxdelta, by = eachdelta)
  
  # Process each gene copy
  for(i in copyname){
    thiscopydata = realdata[[i]]
    # Skip processing if all data points are NA
    if(all(is.na(thiscopydata))){
      next
    }
    # Process each cell
    for(j in cellname){
      choosedata = realdata %>% filter(cell == j)
      # Check if there's any non-NA data for this cell and gene copy
      if(any(!is.na(choosedata[[i]]))){
        oldtime = choosedata$time
        # Filter out NA times for processing
        if(any(is.na(choosedata[[i]]))){
          nonaindex = which(!is.na(choosedata[[i]]))
          oldtime = oldtime[nonaindex]
        }
        ini_point = floor(min(oldtime))

        newx = oldtime[1]
        newtime = c(sort(newx), oldtime)
        
        # Handle cases where interpolated times are earlier than initial points
        if(any(newtime < ini_point)){
          chooindex = which(newtime < ini_point)
          newtime_1 = newtime[-chooindex]
          newtime_2 = newtime[chooindex]
          lastcell = str_sub(j, 1, nchar(j) - 1)
          
          fcg1 = allfun[[i]][[j]]
          allhash[[i]][[j]] = hash(keys = round(newtime_1, 3), values = fcg1(newtime_1))
          
          fcg2 = allfun[[i]][[lastcell]]
          .set(allhash[[i]][[lastcell]], keys = round(newtime_2, 3), values = fcg2(newtime_2))
        } else {
          # Store interpolated values in hash table
          fcg = allfun[[i]][[j]]
          allhash[[i]][[j]] = hash(keys = round(newtime, 3), values = fcg(newtime))
        }
      }
    }
  }
  
  # Save the hash tables to an RData file
  save(allhash, file = paste(savepath, subtree, "_", maxdelta_time, "_allhash.RData", sep = ""))
}

# (6) computeM ------------------------------------------------------------
# Define the function computeM to calculate the median of all interpolated values stored in hash tables
# Inputs:
#   hashdata: List of hash tables, contains interpolated data indexed by gene copy and cell
#   subtree: String, identifier for the subtree used in naming saved files
#   maxdelta_time: Numeric, maximum delta time
#   savepath: String, directory path where the median values will be saved
# Output:
#   None directly; function saves the computed median values to specified path
computeM = function(hashdata, subtree, maxdelta_time, savepath){
  
  # Define a helper function to extract and convert hash table values to numeric for a single hash table
  firstlist_fun = function(subfirst){
    tempvalue = as.numeric(values(subfirst))
    return(tempvalue)
  }
  
  # Define another helper function to apply the first helper function to each sublist (each cell in a copy)
  # and calculate the median of all numeric values retrieved
  secondlist_fun = function(sublist){
    tempvalue = sapply(sublist, firstlist_fun)  # Apply firstlist_fun to each cell's data
    return(median(unlist(tempvalue)))  # Calculate median of all values across cells
  }
  
  # Apply the second helper function to each gene copy's data in hashdata and compute the median value
  M = sapply(hashdata, secondlist_fun)
  
  # Save the computed median values to an RData file, naming it according to the subtree and maxdelta_time
  save(M, file = paste(savepath, subtree, "_M_", maxdelta_time, "_.RData", sep = ""))
}



# (7) mergeCopyFun --------------------------------------------------------
# Define the function mergeCopyFun to process, merge, and save genetic data based on biovalue thresholding and other criteria
# Inputs:
#   allhash: List of hash tables containing interpolated gene expression data indexed by gene copy and cell
#   genecopy_h: Hash table mapping gene copies to gene names
#   M: Numeric vector, contains median values for each gene copy
#   realdata: List of data frames, each containing real data for a gene copy
#   GeneDataset: Data frame, contains gene expression data
#   InfoDataset: Data frame, contains cell metadata and time points
#   subtree: String, identifier for the subtree used in naming saved files
#   deletecellrate: Numeric, threshold for deleting cells based on some criterion
#   deletegenerate: Numeric, threshold for generation-based deletion criteria
#   copyname: Vector, list of gene copy identifiers
#   cellname: Vector, list of cell identifiers
#   mergerate: Numeric, threshold for merging gene data
#   savepath: String, directory path where processed files will be saved
# Output:
#   None directly; function saves processed genetic data and related statistics to specified path
mergeCopyFun = function(allhash, genecopy_h, M, realdata, GeneDataset, InfoDataset, subtree, deletecellrate,
                        deletegenerate, copyname, cellname, mergerate, savepath){
  
  # Process each gene copy's hash table data to threshold values based on median M
  for(i in 1:length(allhash)){
    M_i = M[i]
    if(is.na(M_i)){
      next
    }
    for(j in cellname){
      if(length(keys(allhash[[i]][[j]])) == 0){
        next
      }
      thekey = keys(allhash[[i]][[j]])
      thevalue = as.numeric(values(allhash[[i]][[j]]))
      biovalue = ifelse(thevalue > M_i, 1, 0)
      allhash[[i]][[j]] = hash(keys = thekey, values = biovalue)
    }
  }
  
  # Save the processed hash data to an RData file
  save(allhash, file = paste(savepath, "mer_", subtree, "_", maxdelta_time, "_bioallhash.RData", sep = ""))

  # Aggregate all gene names from gene copies
  allcopygene = c()
  for(i in copyname){
    thiscopygene = genecopy_h[[i]]
    allcopygene = c(allcopygene, thiscopygene)
  }
  uniquegene = sort(unique(allcopygene))

  # Initialize a new list for storing merged hash data
  newallhash = list()
  for(i in uniquegene){
    colindex = which(allcopygene == i)
    copynum = length(colindex)
    nacopy = 0
    for(q in colindex){
      thiscopy = copyname[q]
      if(all(is.na(GeneDataset[[thiscopy]]))){
        nacopy = nacopy + 1
        next
      }
    }
    if(nacopy == copynum){
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
        if(length(thiskeys) == 0){
          next
        }
        tempdata = data.frame(time = as.numeric(thiskeys), thisvalues)
        names(tempdata) = c("time", copyname[j])
        mergedata = mergedata %>% full_join(tempdata, by = "time")
      }
      if(nrow(mergedata) == 1){
        next
      }
      mergedata = mergedata %>% filter(time > 0) %>% arrange(time) %>%
        mutate(blot = rowMeans(select(., contains(i)), na.rm = TRUE)) %>%
        mutate(blot = ifelse(blot >= mergerate, 1, 0))
      newallhash[[i]][[k]] = hash(keys = mergedata$time, values = mergedata$blot)
    }
  }
  
  # Save the new merged hash data to an RData file
  save(newallhash, file = paste(savepath, "mer_", subtree, "_", maxdelta_time, "_origenehash.RData", sep = ""))

  # Select names of genes from the newly merged hash table
  selectgene = names(newallhash)

  # Create a new dataset initialized to zeros, with rows equal to the GeneDataset and columns to selected genes
  bioGeneDataset = matrix(0, nrow = nrow(GeneDataset), ncol = length(selectgene))
  bioGeneDataset = as_tibble(bioGeneDataset)
  names(bioGeneDataset) = selectgene

  # Define a function to extract gene values for each cell from the hash tables
  eachgenevalue = function(eachgene, newallhash, InfoDataset, cellname){
    thisvale = rep(0, nrow(InfoDataset))
    for(i in cellname){
      chooseindex = which(InfoDataset$cell == i)
      choosetime = round(InfoDataset$time[chooseindex], 3)
      if(is.null(newallhash[[eachgene]][[i]])){
        thisvale[chooseindex] = rep(NA, length(chooseindex))
      } else {
        thisvale[chooseindex] = as.numeric(values(newallhash[[eachgene]][[i]], keys = choosetime))
      }
    }
    return(thisvale)
  }

  # Apply the function to each gene to populate the bioGeneDataset
  bioGeneDataset[,] = sapply(selectgene, eachgenevalue, newallhash, InfoDataset, cellname)

  # Write the populated dataset to a CSV file
  write_csv(bioGeneDataset, file = paste(savepath, "mer_", subtree, "_", maxdelta_time, "_taskdata.csv", sep = ""))

  # Calculate statistics for each gene and cell combination, then reshape and save the data
  newstatdata = InfoDataset[,1] %>% add_column(bioGeneDataset) %>%
    gather(key = "genename", value = "num", -cell) %>% group_by(cell, genename) %>%
    summarise(count = sum(!is.na(num)), .groups = "drop") %>% spread(key = genename, value = count)
  write_csv(newstatdata, file = paste(savepath, "mer_", subtree, "_", maxdelta_time, "_statdata.csv", sep = ""))

  # Filter cells and genes based on defined rates and generate a final selected dataset
  selectinfo = delectFun(newstatdata, deletecellrate, deletegenerate)
  selectcell = selectinfo$selectcell
  selectgene = selectinfo$selectgene

  selectindex = which(names(bioGeneDataset) %in% selectgene)
  newbioGeneDataset = bioGeneDataset[, selectindex]
  mergedata = InfoDataset %>% add_column(newbioGeneDataset) %>%
    filter(cell %in% selectcell)

  # Assign tree indices and handle cell changes across generations
  allcellname = mergedata$cell
  alltreenum = rep(1, nrow(mergedata))
  nowmother = allcellname[1]
  for(i in 2:nrow(mergedata)){
    if(str_sub(allcellname[i], 1, nchar(nowmother)) == nowmother){
      alltreenum[i] = alltreenum[i - 1]
    } else {
      nowmother = allcellname[i]
      alltreenum[i] = alltreenum[i - 1] + 1
    }
  }

  tree_h = hash(allcellname, alltreenum)
  mergedata$treeindex = 0
  for(i in 1:nrow(mergedata)){
    mergedata$treeindex[i] = tree_h[[mergedata$cell[i]]]
  }

  newmergedata = mergedata %>% arrange(cell, time) %>%
    mutate(originrow = 1:n(), lastrow = originrow - 1) %>%
    group_by(treeindex) %>% mutate(ifchange = as.numeric(lag(cell) != cell))  %>%
    select(cell, time, originrow, lastrow, treeindex, ifchange, everything()) %>% ungroup() %>%
    mutate(ifchange = replace_na(ifchange, 9))

  # Handle changes at the start of new cell generations and save the processed data
  changeindex = which(newmergedata[["ifchange"]] == 1)
  c_cell = newmergedata$cell
  lastrow = newmergedata$lastrow
  for(i in changeindex){
    thiscell = c_cell[i]
    lastcell = substr(thiscell, 1, nchar(thiscell) - 1)
    choosedata = newmergedata %>% filter(cell == lastcell) %>% select(originrow)
    temprow = choosedata$originrow
    lastrow[i] = temprow[length(temprow)]
  }
  newmergedata$lastrow = lastrow

  write_csv(newmergedata, file = paste(savepath, "mer_", subtree, "_", maxdelta_time, "_selectdata.csv", sep = ""))

  # Calculate starting rows for cell generations and save the information for reference
  allstartrow = which(newmergedata[["ifchange"]] == 9)
  allstartrow = as.numeric(sapply(allstartrow, function(x) seq(from = x, length.out = 5)))
  allchooserow = setdiff(1:nrow(newmergedata), allstartrow)

  save(allstartrow, file = paste(savepath, "mer_", subtree, "_", maxdelta_time, "_allstartrow.RData", sep = ""))
  save(allchooserow, file = paste(savepath, "mer_", subtree, "_", maxdelta_time, "_allchooserow.RData", sep = ""))

}

# (8) delectFun -----------------------------------------------------------
# Define the function delectFun to filter cells and genes based on specific criteria
# Inputs:
#   newstatdata: Data frame, contains statistics for each gene in each cell
#   deletecellrate: Numeric, threshold rate for deleting cells with a high proportion of zeros
#   deletegenerate: Numeric, threshold rate for deleting genes with a high proportion of zeros
# Output:
#   List containing vectors of selected cell names and gene names after filtering
delectFun = function(newstatdata, deletecellrate, deletegenerate){
  # Extract gene names (excluding the first column which contains cell names)
  genename = names(newstatdata)[-1]
  # Extract cell names from the first column
  cellname = newstatdata[[1]]
  # Convert the data frame to a matrix excluding the first column of cell names
  tempdata = as.matrix(newstatdata[,-1])
  
  # Calculate the proportion of zeros in each row (cell) and identify rows exceeding the deletion threshold
  rowstat = apply(tempdata, 1, function(x) sum(x == 0) / length(x))
  derowindex = which(rowstat > deletecellrate)
  if(length(derowindex) > 0){
    selectcell = cellname[-derowindex]
  } else {
    selectcell = cellname
  }

  # Filter the dataset to include only the selected cells
  newstatdata = newstatdata %>% filter(cell %in% selectcell)
  tempdata = as.matrix(newstatdata[,-1])
  
  # Calculate the proportion of zeros in each column (gene) and identify columns exceeding the deletion threshold
  colstat = apply(tempdata, 2, function(x) sum(x == 0) / length(x))
  decolindex = which(colstat > deletegenerate)
  if(length(decolindex) > 0){
    selectgene = genename[-decolindex]
  } else {
    selectgene = genename
  }

  # Select the data for the remaining genes after filtering
  selectdata = newstatdata %>% select(all_of(selectgene))
  
  # Return a list containing the names of selected cells and genes
  return(list(selectcell = selectcell, selectgene = selectgene))
}
