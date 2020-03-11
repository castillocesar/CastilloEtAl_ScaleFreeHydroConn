# Castillo <i>et.al.</i> Scale Free Hydrologic Connectivity
This set of codes is provided as part of the methods section for the manuscript entitled "Scale-free
structure of surface-water connectivity within a lowland river floodplain" by Cesar R Castillo,
İnci Güneralp, Billy Hales, and Burak Güneralp. The corresponding author is Cesar R Castillo
and can be contacted at castillocesar@tamu.edu.

This set of scripts is designed to be used in four steps using the Python and R computing languages
and brief descriptions of each are listed below. More detailed/technical descriptions can be found
within the documentation/comments of the scripts themselves. The data that can be used to reproduce
the analysis performed by Castillo et al. can be found at the following URL 
https://dataverse.tdl.org/dataverse/CastilloEtAl_ScaleFreeHydroConn. Data can also be provided by
the corresponding author upon request.

1. The first script used in the analysis is
CastilloEtAl_ScaleFreeHydroConn_DetermineFloodExtent_GenerateNetwork.py and it forms the basis of the
GIS analysis performed in our manuscript.This Python script uses tif rasters of inundation depth
exported from HEC-RAS to determine the surface-water connections between landscape patches in the
river-floodplain landscape. The landscape patches in this particular case are composed of an irregular
tessalation that are represented by a polygon shapefile/feature class with a number of specific fields
in the attribute table. The connectivity analysis is only really focused on portions of the landscape
that are connected to the main channel by surface water and this script also requires a polyline
shapefile/feature class the centerline or thalwag for the main channel to filter out inundation that
is not connected to the main channel. This code was designed to be used with Python 3.6 and ArcGIS
Pro 2.x.

2. The second script in the analysis is 
CastilloEtAl_ScaleFreeHydroConn_GenerateGraphObjects_CalculateConnectivityMetrics_GenerateDegreeDistributions.R
and it forms the basis for generating the graph/network structure and it generates/computes the degree 
distributions and topological/algebraic metrics that we use in our analysis. This R script uses output
from the CastilloEtAl_ScaleFreeHydroConn_DetermineFloodExtent_GenerateNetwork.py Python script that
performs the GIS analysis to create the graph data structures and perform the topological and
algebraic analysis for each iteration of the netowrk. This code will iterate over a set of tables in
comma-separated format will columns that will be used to create and edge-list that will be the basis
for creating the graph objects over which this analysis is based. The outputs from this script include
a table of the topological and algebraic metrics for all network states that were iterated over and a
degree distribution table for each network iteration. This code was designed to be used with the "igraph"
library within R. igraph_1.2.x was the version that was used in the analysis.

3. The third script in the analysis is
CastilloEtAl_ScaleFreeHydroConn_CreateDegreeSequenceFromDegreeDistribution_ForBroidoAndClauset2019Code.r
and the R scripts uses the degree distributions outputed from the R script
CastilloEtAl_ScaleFreeHydroConn_GenerateGraphObjects_CalculateConnectivityMetrics_GenerateDegreeDistributions.R
to create a degree sequence that can be used as an input into the Python code that performs the
empirical analysis for determing the scale-freeness of a degree distribution. This code in particular,
loops through a set of degree distributions to generate the seperate degree distributions and outputs
them as csv files in another directory.

4. The fourth script in the analysis is 
CastilloEtAl_ScaleFreeHydroConn_ScaleFreeTestingOfDegreeDistribution_CombiningBroidoAndClauset2019Code.py
and the Python script uses the degree sequences created by the 
CastilloEtAl_ScaleFreeHydroConn_CreateDegreeSequenceFromDegreeDistribution_ForBroidoAndClauset2019Code
code and subjects it to code provided by Broido and Clauset in their 2019 paper published in
Nature Communications under the title "Scale-free networks are rare". This script combines all their
code and reads in the degree sequence and conducts an empirical analysis on the degree distribution
that is used to determine how scale-free a network is. The original Broido and Clauset 2019 code can
be found at the following URL https://github.com/adbroido/SFAnalysis/tree/master/code. This code was
designed to be used in Python 3.6 with the most up-to-date version of the following Python libraries:
pandas, numpy, scipy, time, mpmath, norm, chi2, and os.
