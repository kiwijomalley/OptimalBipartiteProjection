# OptimalBipartiteProjection
R and Python Code for paper under review at the Journal of Data Science. The current version of the scripts are setup to run directly on the data set for our motivating example. This contains a lot of supurfluous information from administrative claims data that is not used in the analyses in the current paper. It is anticipated that the script will in time be simplified to a barebones version. Note that the data wrangling and data analysis are performed in distinct steps because the raw data is patient confidential and can only be analyzed behind a firewall.

Python steps

Save your network in a data set 'Network.txt'. The key columns are: <br/>
 Col1 = BeneID = Beneficiary ID <br/>
 Col2 = SFROMDT = Date of patient-physician encounter <br/>
 Col3 = PRFNPI = Physician ID <br/>
 Col4 = WorkRVU = Measure of average amount of physician effort related to encounter 
 
Step 1: Run Script1_FormNetworks.py to form the 20 directed and 20 undirected weighted networks. This generates 3 output files.
Step 2: Run Script2_BinaryFilter.py to obtain the binary versions of the networks computed in Step 1
Step 3: Run Script3_BatchSubNetworks.py to break the network into components corresponding to study units (in our case, the 4800 hospitals) and combines the weighted and the binary projections together (80 in total). A separate file that shows the individuals (in our case, the physicians) in each study unit is needed. If you only have a single network change nbatches from 160 to 1.
Step 4: Run Script4_SummaryMeasures.py to compute the 16 network summary measures for each network and each of the 80 projections

R steps

The networks formed under the 80 different methods are separately analyzed with respect to Discrimimatory Power and Predictive Power. Within these, two R scripts are run. The first yields estimates of network-wide and measure-specific measures of Discrimimatory Power and Predictive Power, respectively. A second script is then run in each case to summarize the results in a variety of ways.  

Discriminatory Power

Step 1: Run DiscPowerModels.r to obtain overall network and measure specific summary measures of Discriminatory Power
Step 2: Run NetDiscPowerAnal.r to summarize the measures of Discriminatory Power in various ways, including those reported in the paper.

Predictive Power

Step 1: Run PredPowerModels.r to obtain estimates of Predictive Power under each projection method and for the model with the base design's measures added as additional predictors
Step 2: Run NetPredPowerAnal.r to obtain absolute and percentage improvement in predictive power for each projection method compared to the base design.  

The names of the network summary measures formed in the Pythin script and analyzed in the R script are: <br/>
"Nodes" = Size of network <br/>
"Density" = Network density <br/>
"Specialization" = Variance of (undirected) degree distribution <br/>
"CorDeg_InOut" = Correlation between actor in-degree and receiver out-degree <br/>
"Assort_InIn" = Correlation between sender in-degree and receiver in-degree <br/>
"Assort_InOut" = Correlation between sender in-degree and receiver out-degree <br/>
"Assort_OutIn" = Correlation between sender out-degree and receiver in-degree <br/>
"Assort_OutOut" = Correlation between sender out-degree and receiver out-degree <br/>
"Reciprocity" = Reciprocity of edge status <br/>
"Transitivity" = Transitivity coefficient <br/>
"Ave_Clustering" = Average of the local (egocentric) clustering coefficient <br/>
"Ncomponents" = Number of components in the network <br/>
"Diameter" = Diameter of the network <br/>
"sumTransitivity" = Total number of transitive triads in the directed network <br/>
"sum3Cycle" = Total number of closed 3-cycles in the directed network
