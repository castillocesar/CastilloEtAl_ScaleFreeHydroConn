################################################################################
# This script is provided as part of the methods section for the manuscript
# entitled "Scale-free structure of surface-water connectivity within a lowland
# river floodplain" by Cesar R Castillo, Inci Güneralp, Billy Hales, and Burak
# Güneralp. The corresponding author is Cesar R Castillo and can be contacted at
# castillocesar@tamu.edu.
#
# This script uses output from the
# CastilloEtAl_ScaleFreeHydroConn_DetermineFloodExtent_GenerateNetwork.py
# Python script that performs the GIS analysis to create the graph data
# structures and perform the topological and algebraic analysis for each
# iteration of the netowrk. This code will iterate over a set of tables in
# comma-separated format will columns that will be used to create and edge-list
# that will be the basis for creating the graph objects over which this analysis
# is based. The outputs from this script include a table of the topological and
# algebraic metrics for all network states that were iterated over and a degree
# distribution table for each network iteration.
#
# It is assumed that all inputs and outputs will be within a source directory
# that is specified below (variable is "prj_dir"). The user can specify within
# the directories and subdirectory where outputs will be sent
#
# This code was designed to be used with the "igraph" library within R.
# igraph_1.2.x was the version that was used in the analysis.
################################################################################


# tracking the speed of the script
start.time <- Sys.time()
# importing libraries
library(igraph)

################################################################################
# specify the filename for the output csv file that will be created
out.filename.csv <- 'FILENAME'
i.start <- 1

# specify the directory for the project
prj.dir <- 'DIRECTORY'
# specify the subdirectory where the tables containing the columns that 
# will be used as the edge-list for constructing the graph objects
wd <- paste(prj.dir,
            '/',
            'SUBDIRECTORY',
            sep='')
# setting the working directory for the project
setwd(wd)

# constructing list of comma-separated files with edge-listings
# the tables directly exported from the python script that does
# the GIS analysis will be comma-separated ".txt" files and
# thus the code below specifically calls these types of files
files <- list.files(pattern='*.txt$')

# creating a vector that depicts the stage height at gaging station
# here we use values that span the historical range for the 08189500
# USGS gaging station on Mission River 
# change as needed for one's own respective analysis
stage <- c(seq(from=1,
               to=11.5,
               by=0.5),
           11.63)

# creating an empty matrix that will hold the sets of topological and
# algebraic graph metrics
conn.udir.smet <- matrix(ncol=11,
                         nrow=1)
# UNDIRECTED GRAPH WITH NO WEIGHTS BEING USED
#for(i in 1:5){
for(i in i.start:length(files)){
#for(i in 1:length(files)){
  # reading in table that contains columns to construct edge-list for the
  # respective iteration
  tbl <- read.csv(files[i],
                     header=TRUE,
                     sep=',')
  # identifying the stage associated with each flow so that it can
  # be included in the filename of the degree matrix that will be
  # outputed
  stage_val <- ifelse(test=(nchar(x=format(x=stage[i])) < 4),
                      yes=(paste('0',
                                 format(x=stage[i]*100),
                                 sep='')),
                      no=(format(x=stage[i]*100)))
  # creating a filename for the degree matrix that will be outputed
  filename <- paste(prj.dir,
                    '/',
                    'graph_analysis_degreematrices',
                    '/',
                    'graph_analysis_degreematrices_',
                    stage_val,
                    '.csv',
                    sep='')
  # removing the rows that have an na in them this shouldn't
  # be an issue if the input data has been cleaned but this
  # is in there cause it can blow things up
  tbl <- na.omit(tbl)
  # creating the edge-list by identifying the source and
  # neighboring patches
  el <- as.matrix(cbind(as.character(tbl$SRC_PID),
                        as.character(tbl$NBR_PID)))
  #converting initial edge-list vector to an undirected igraph/graph object
  gph <- graph.edgelist(el=el,
                        directed=FALSE)
  # finding the adjacency matrix of the graph
  adj <- get.adjacency(graph=gph,
                       type='both',
                       edges=FALSE)
  # finding the lapacian matrix of the graph
  lap <- graph.laplacian(graph=gph,
                         normalized=FALSE,
                         weights=NULL,
                         sparse=getIgraphOpt('sparsematrices'))
  
  # finding the network order
  n <- nrow(adj)
  
  # finding the network size
  m <- gsize(graph=gph)
  # finding the min number of edges for the current number of vertices
  m.min <- n - 1
  # finding the max number of edges for the current number of vertices
  m.max <- (n^2 - n) / 2
  # finding the normalized network size
  m.norm <- (m - m.min) / (m.max - m.min)
  
  # find the degree matrix for the graph
  deg <- degree(graph=gph)
  # finding the mean vertex degree
  deg.mean <- mean(deg)
  
  # finding the mean geodesic path (l)
  l <- mean_distance(graph=gph,
                     directed=FALSE,
                     unconnected=TRUE)
  
  # find the transitivity of the graph (clustering coefficient)
  tr <- transitivity(graph=gph,
                     type='undirected')
  
  # finding the graph level centrality index score
  cb <- centr_betw(graph=gph,
                   directed=FALSE,
                   nobigint=TRUE,
                   normalized=TRUE)
  
  # decomposing the adjacency matrix to find the spectral radius (largest eigenvalue)
  sr <- max(eigen(x=adj,
                  symmetric=TRUE)$values)
  # maximum theoretical value for spectral radius for a given n
  sr.max <- n - 1
  # theoretical upper bound for the spectral radius for a given n and m
  sr.upper <- ((2 * m * (n - 1)) / n)^(0.5)
  # computing the ratios between sr and its theoretical bounds
  sr.srmax.ratio <- sr / sr.max
  sr.srupper.ratio <- sr / sr.upper
  
  # finding the relative contribution to sr from graph connectivity
  rel.cont.conn <- (sr.max - sr.upper) / (sr.max - sr)
  # finding the relative contribution from the wiring
  rel.cont.wir <- 1 - rel.cont.conn
  
  # decomposing the laplacian matrix to find the algebraic connectivity (2nd smallest eigenvalue)
  ac <- sort(eigen(x=lap,
                   symmetric=TRUE)$values,
             partial=n-n+2)[n-n+2]
  # finding the min theoretical value for ac
  gph.diameter <- diameter(graph=gph,
                           directed-FALSE,
                           unconnected=TRUE)
  ac.min <- 4 / (n * gph.diameter)
  # finding the max theoretical value for ac
  ac.max <- vertex.connectivity(graph=gph)
  # computing the normalized ac value
  ac.norm <- (ac - ac.min) / (ac.max - ac.min)
  
  # creating the degree matrix data frame and adding column labels
  # this is what will be exported to its own csv file
  deg <- cbind(as.numeric(labels(deg)),
               as.numeric(deg))
  deg <- data.frame(deg)
  colnames(deg) <- c('id', 'degree')
  
  # putting together the matrix that will hold the values for each 
  # iteration of the network
  ifelse(i <= i.start,
         grph.met <- cbind(n, m.norm, deg.mean,
                           l, tr, cb,
                           sr.srupper.ratio, sr.srmax.ratio, rel.cont.conn, rel.cont.conn, ac.norm),
         graph.met <- rbind(gph.met, cbind(n, m.norm, deg.mean,
                                           l, tr, cb,
                                           sr.srupper.ratio, sr.srmax.ratio, rel.cont.conn, rel.cont.conn, ac.norm))
         )

  
  # writing the degree data frame for each iteratoin to file as a csv
  write.table(x=deg,
              file=filename,
              sep=",",
              row.names=FALSE,
              col.names=TRUE)
  
  # deleted variables in order to free up memory for next iteration in loop
  remove(el, gph, adj, lap,
         n, m, m.min, m.min,
         m.norm, deg.mean, l, tr,
         cb, sr, sr.max, sr.upper,
         sr.srupper.ratio, sr.srmax.ratio, rel.cont.wir, rel.cont.conn,
         ac, ac.min, ac.max, ac.norm)
  
  cat('finished with iteration', i, 'of', length(files),
      'from the graph analysis at', 
      as.character(Sys.time()))
  cat('\n')
  
}

gph.met <- cbind(stage, gph.met)
gph.met <- data.frame(gph.met)
gph.met <- round(x=gph.met,
                 digits=3)
colnames(conn.udir.smet) <- c('stage', 'n', 'm_norm',
                              'deg_mean', 'l', 'C',
                              'Cbtw', 'sr_srupper_ratio', 'sr_srmax_ratio',
                              'rel_cont_conn', 'rel_cont_wir', 'ac_norm')

filename <- paste(prj.dir,
                  '/',
                  filename.csv,
                  sep='')
write.table(gph.met,
            file=filename,
            sep=',',
            row.names=FALSE,
            col.names=TRUE)

end.time <- Sys.time()
cat('graph analysis started at', as.character(start.time))
cat('\n')
cat('graph analysis ended at', as.character(end.time))
cat('\n')
cat('processing time for graph analysis is', end.time - start.time)
cat('\n')
remove(start.time, end.time)

# end of the script

