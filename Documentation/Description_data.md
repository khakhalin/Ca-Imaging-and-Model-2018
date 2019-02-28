Comments on how data is stored
============

# Main analysis

### all_cells_latencies.csv

* Contains: info on cell position and response latency, for space-related plots.

### avamps_allcells_allbrains stable.csv

* Contains: response amplitudes for every cell recorded in every experiment.

TODO Q: Is it truly all cells (before we restrict analysis to the largest component), or is it post-restriction?

### caimg_network_measures.csv

* Contains: global network measures (Figure 4)
* Produced by:
* Used by: caimg_network_analysis.r

Columns: 

1. type - one of two groups of values:
  * 1actual - actual measurements on a reconstructed graph. Numbers before the name are so that R would sort them properly by default.
  * 2rewire - an average of the same value, across 100 degree-preserving rewires of the reconstructed graph
1. stage - tadpole stage (46 or 49)
1. name - experiment code (YYMMDD + a letter if more than one experiment in a given day)
1. eff - network efficiency
1. a_oi - Out-In assortativity
1. a_io - In-Out assortativity
1. a_oo - Out-Out assortativity
1. a_ii - In-In assortativity
1. modul - Spectral modularity
1. clust - Global clustering coefficient
1. flow - Hierarchical flow (see the text)

### caimg_network_summary.csv

* Contains: Basic low-level information about each reconstructed network
* Produced by:
* Used by: caimg_experimens_assorted_figures.r

Columns:

1. name - experiment name
1. stage - tadpole stage
1. nCells - total number of cells quantified
1. nSweeps - total number of sweeps (stimulus presentations, across all 3 types of stimuli) in this experiment
1. ampF - average response amplitude for Flashes
1. ampS - same for Scrambled
1. ampC - same for Looming (C stands for "Crash", which is an old internal term for Looming) 
1. nEdges - number of edges detected
1. largest - size of the largest weakly connected component
1. encoding - how well stimulus identity can be reconstructed from this recording (see the text)
1. gamma - degree distribution power, for a scale-free fit (see the text)

### sel_allcells_allbrains stable.csv

* Contains: selectivity measurements for every cell
* Produced by:
* Used by:

Columns:

1. fc - Flash-Looming selectivity
1. fs - Flash-Scrambled selectivity
1. sc - Scrambled-Looming selectivity
1. ibrain - experiment number
1. stage - developmental stage for this experiment

### sel_centrality_allcells.csv

* Contains: Centrality measurements for each cell, as well as cell selectivity
* Produced by: graph_structure_analyzer.m (but only if doSelectivity flag is set to TRUE). 
* Used by: caimg_centrality_analysis.r

Columns: 

1. nexp - experiment number
1. sel - Flash-Looming selectivity
1. indegree - in-degree of this node (cell)
1. katz - Katz centrality of this node (cell)
1. spiking - average amplitude of spiking response in this cell
1. insel - average selectivity of all cells connected to this cell (in-connections)

Currently does not contain cell ids for two reasons. One, it is the only secondary data file that contains information on per-cell level. Second, only nodes from the largest weakly connected component are included in this analysis.

First 14 experiments are stage 46; the reast are stage 49.


# Original dataset info

For each experiment, we originally had the following set of data files:

1.	Excel file with fluorescence signals
2.	TIFF containing ROIs – Not necessary, as ROI (x,y) coordinates are also stored in same Excel data file.
3.	A print-screen of ROIs (for further processing in imageJ)
4.	.mat file with a Matlab structure (S) containing all data from this Excel file, and also spikes reconstructed
