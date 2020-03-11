################################################################################
# This script is provided as part of the methods section for the manuscript
# entitled "Scale-free structure of surface-water connectivity within a lowland
# river floodplain" by Cesar R Castillo, İnci Güneralp, Billy Hales, and Burak
# Güneralp. The corresponding author is Cesar R Castillo and can be contacted at
# castillocesar@tamu.edu.
#
# This script uses the degree distributions outputed from the R script
# CastilloEtAl_ScaleFreeHydroConn_GenerateGraphObjects_CalculateConnectivityMetrics_GenerateDegreeDistributions.R
# to create a degree sequence that can be used as an input into the Python code
# that performs the empirical analysis for determing the scale-freeness of a
# degree distribution. This code in particular, loops through a set of degree
# distributions to generate the seperate degree distributions and outputs them
# as csv files in another directory.
#
# It is assumed that all inputs are within the same directory and all outputs will 
# be saved in a separate directory
################################################################################
# defining a function that creates a degree sequence from a given degree matrix
create.degree.seq <- function(degree.vector){
  deg <- table(degree.vector)
  deg <- cbind(as.numeric(names(deg)),
               as.numeric(deg))
  deg <- as.data.frame(deg)
  cnames <- c('xvalue',
              'counts')
  colnames(deg) <- cnames
  deg.vec <- seq(from=1,
                 to=max(deg$xvalue),
                 by=1)
  counts.vec <- rep(x=0,
                    times=length(deg.vec))
  deg.seq <- cbind(deg.vec, counts.vec)
  deg.seq <- as.data.frame(deg.seq)
  colnames(deg.seq) <- cnames
  for(i in 1:nrow(deg.seq)){
    for(j in 1:nrow(deg)){
      if(deg.seq$xvalue[i] == deg$xvalue[j]){
        deg.seq$counts[i] <- deg$counts[j]
      }
    }
  }
  return(deg.seq)
}

################################################################################
# specifying the directories and filenames

# specifying the general directory where the inputs and outputs will be held
gen.dir <- 'DIRECTORY'
# specifying the subdirectory within the general directory where the inputs
# are being stored
deg.mat.dir <- paste(gen.dir,
                     '/',
                     'SUBDIRECTORY',
                     sep='')
# creating a list of the files that contains the degree distributions
# the degree distributions provided by the 
# GenerateGraphObjects_CalculateConnectivityMetrics_GenerateDegreeDistributions
# will be in csv format
deg.mat.files <- list.files(path=deg.mat.dir,
                            pattern='*.csv')
# specify the subdirectory where the outputs will be sent
deg.seq.dir <- paste(gen.dir,
                     '/',
                     'SUBDIRECTORY',
                     sep='')

####################################################################################
# looping through the set of degree distribution files and outputing the degree
# sequences into a different subdirectory
for(i in 1:length(deg.mat.files)){
  deg.mat.file <- paste(deg.mat.dir,
                        '/',
                        deg.mat.files[i],
                        sep='')
  # reading in the degree matrix
  dat <- read.csv(file=deg.mat.file,
                  header=TRUE,
                  sep=',')
  dat.seq <- create.degree.seq(dat$degree)
  # adding the stage in the degree distribution filename to the
  # degree sequence files
  # here we use values that span the historical range for the 08189500
  # USGS gaging station on Mission River 
  # change as needed for one's own respective analysis
  stage <- substr(x=strsplit(x=deg.mat.files[i],
                             split='_')[[1]][4],
                  start=1,
                  stop=4)
  dat.seq.filename <- paste(deg.seq.dir,
                            '/',
                            'degree_sequence_',
                            stage,
                            '.csv',
                            sep='')
  write.table(x=dat.seq,
              file=dat.seq.filename,
              col.names=TRUE,
              row.names=FALSE,
              sep=',')
}



