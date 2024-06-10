# (1) allrelationFun ------------------------------------------------------
# Define the function allrelationFun to process gene interaction data and generate a connection table.
# Inputs:
#   subtree: Character, name of the subtree being analyzed
#   maxdelta_time: Numeric, maximum delta time value
#   p: Numeric, specific parameter value used in file naming
#   realpath: Character, path to the directory containing the real data files
#   realresultpath: Character, path to the directory containing the result files
# Output:
#   DataFrame, containing the connection table with columns from, to, preaij, and predelta
allrelationFun = function(subtree, maxdelta_time, p, realpath, realresultpath){
  
  # Read the selected data from the CSV file
  newdata = read_csv(file=paste0(realpath, "mer_", subtree, "_", maxdelta_time, "_selectdata.csv"), show_col_types = FALSE)
  
  # Extract final information and data
  finalinfo = newdata[, 1:6]
  finaldata = newdata[, -c(1:6)]
  genename = names(finaldata)
  
  # Initialize the connection data frame
  connect = data.frame(from = "a", to = "b", preaij = 0, predelta = 0)
  
  # Load the simulation information
  load(file=paste0(realpath, subtree, "_simu_info.RData"))
  Mc = simu_info$Mc
  eachdelta = 1 / (round(Mc) - 1)
  
  # Read the precomputed interaction matrices
  preaij = read.csv(file=paste(realresultpath, subtree, "_", p, "_preaij.csv", sep=""))
  predelta = read.csv(file=paste(realresultpath, subtree, "_", p, "_predelta.csv", sep=""))
  
  # Loop through each gene to find interactions and update the connection table
  for(i in 1:length(genename)){
    nozero = which(preaij[i, ] != 0)
    if(length(nozero) == 0){
      next
    } else {
      for(j in nozero){
        if(preaij[i, j] == -0.1){
          preaij[i, j] = -1
        }
        tempconnect = data.frame(from = genename[j],
                                 to = genename[i],
                                 preaij = preaij[i, j],
                                 predelta = round(predelta[i, j] / eachdelta))
        connect = rbind(connect, tempconnect)
      }
    }
  }
  
  # Remove the initial dummy row from the connection table
  connect = connect[-1, ]
  
  # Print the unique delta values
  print(sort(unique(connect$predelta)))
  
  # Return the connection table
  return(connect)
}

# (2) singleNetFun ----------------------------------------------------------
# Define the function singleNetFun to generate a circular network plot for gene interactions.
# Inputs:
#   subtree: Character, name of the subtree being analyzed
#   connect: DataFrame, containing the connection table with columns from, to, preaij, and predelta
#   maxdelta_time: Numeric, maximum delta time value
#   realpath: Character, path to the directory containing the real data files
#   realresultpath: Character, path to the directory containing the result files
# Output:
#   Plot, circular network plot for gene interactions
singleNetFun = function(subtree, connect, maxdelta_time, realpath, realresultpath){
  
  # Read the selected data from the CSV file
  newdata = read_csv(file=paste0(realpath, "mer_", subtree, "_", maxdelta_time, "_selectdata.csv"), show_col_types = FALSE)
  
  # Extract final information and data
  finalinfo = newdata[, 1:6]
  finaldata = newdata[, -c(1:6)]
  genename = names(finaldata)
  
  # Create a vertices data frame with unique gene names and their counts
  c(as.character(connect$from), as.character(connect$to)) %>%
    as_tibble() %>%
    group_by(value) %>%
    summarize(n = n()) -> vertices
  colnames(vertices) <- c("name", "n")
  
  # Filter the connection table to include only those vertices present in the vertices data frame
  connect <- connect %>%
    filter(from %in% vertices$name) %>%
    filter(to %in% vertices$name) %>%
    left_join(vertices, by = c('from' = 'name'))
  
  # Calculate angles and positions for the circular layout
  number_of_bar = nrow(vertices)
  vertices$id = seq(1, nrow(vertices))
  angle = 360 * (vertices$id - 0.5) / number_of_bar
  vertices$hjust = ifelse(angle > 180, 1, 0)
  vertices$angle = ifelse(angle > 180, 90 - angle + 180, 90 - angle)
  
  # Create the graph object
  mygraph <- graph_from_data_frame(connect, vertices = vertices, directed = TRUE)
  
  # Define colors for the plot
  mycolor = brewer.pal(8, 'Dark2')[c(1:6)]
  
  # Generate the circular network plot
  g1 = ggraph(mygraph, layout = "linear", circular = TRUE) +
    geom_edge_arc(aes(edge_colour = factor(predelta, levels = c(0, 1, 2, 3, 4, 5)),
                      linetype = factor(preaij, levels = c(1, -1))), edge_alpha = 0.8, 
                  edge_width = 0.3,
                  arrow = arrow(length = unit(1.5, "mm"), type = "open", ends = "last"),
                  start_cap = square(1.5, 'mm'),
                  end_cap = circle(1.5, 'mm')) +
    scale_edge_linetype_manual(values = c("solid", "dashed"), labels = c("Positive", "Negative"), drop = FALSE) +
    scale_edge_color_manual(values = mycolor,
                            labels = c("0", expression(Delta * t), expression(2 * Delta * t),
                                       expression(3 * Delta * t), expression(4 * Delta * t),
                                       expression(5 * Delta * t)), drop = FALSE) +
    geom_node_text(aes(x = x * 1.16, y = y * 1.16, label = name, angle = angle, hjust = hjust), size = 3, fontface = 'bold') +
    geom_node_point(aes(shape = 1), fill = "#ffa400", shape = 21, size = 2,
                    color = '#ffa400', alpha = 0.6) +
    scale_size_continuous(range = c(0.5, 10)) +
    scale_fill_manual(values = mycolor) +
    scale_color_manual(values = mycolor) +
    expand_limits(x = c(-1.6, 1.6), y = c(-1.6, 1.6)) +
    coord_fixed() +
    theme_minimal() +
    labs(edge_colour = "Time delay",
         edge_linetype = "Regulation types") +
    theme(
      legend.title = element_text(colour = 'black', size = 14, face = "bold"),
      legend.text = element_text(margin = margin(l = 10), colour = "black", size = 12, face = "bold", hjust = 0),
      legend.direction = 'vertical',
      legend.box.background = element_rect(fill = NA, color = "black", linetype = 1),
      panel.grid = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "null"),
      panel.spacing = unit(c(0, 0, 0, 0), "null")) +
    annotate("text", x = 0, y = 0,
             label = paste(subtree),
             colour = "black", size = 6, fontface = "bold")
  
  # Return the generated plot
  return(g1)
}

# (3) multiNetFun ----------------------------------------------------------
# Define the function multiNetFun to generate a composite plot of circular network plots for multiple subtrees.
# Inputs:
#   allsubtree: Character vector, names of the subtrees being analyzed
#   p: Numeric, specific parameter value used in file naming
#   maxdelta_time: Numeric, maximum delta time value
#   readypath: Character, path to the directory containing the real data files for all subtrees
#   realresultmainpath: Character, path to the directory containing the result files for all subtrees
#   figurepath: Character, path to the directory where the resulting figure will be saved
# Output:
#   None, but saves a composite plot of circular network plots as a PDF file
multiNetFun = function(allsubtree, p, maxdelta_time, readypath, realresultmainpath, figurepath){
  
  # Initialize a data frame to store connection data for all subtrees
  alltreedata = data.frame(from = "a", to = "b", preaij = 1, predelta = 1, subtree = "AB")
  
  # Loop through each subtree to generate individual network plots
  for(i in 1:5){
    subtree = allsubtree[i]
    realpath = paste(readypath, subtree, "/", sep = "")
    realresultpath = paste(realresultmainpath, subtree, "/", sep = "")
    
    # Get the connection data for the current subtree
    connect = allrelationFun(subtree, maxdelta_time, p, realpath, realresultpath)
    
    # Generate the network plot for the current subtree and store it in a variable
    assign(paste("g", i, sep = ""), singleNetFun(subtree, connect, maxdelta_time, realpath, realresultpath) + 
             theme(legend.position = "none"))
  }
  
  # Extract the legend from the first plot
  legend <- g_legend(g1 + theme(legend.position = 'right'))
  
  # Create a composite plot with the individual network plots and the legend
  plot_grid(g1, g2, g3, g4, g5, legend, labels = c('(a)', '(b)', '(c)', '(d)', '(e)'), nrow = 2, label_size = 14, hjust = 0, scale = 1)
  
  # Save the composite plot as a PDF file
  ggsave(file = paste0(figurepath, "network.pdf"), width = 240, height = 180, units = "mm")
}
