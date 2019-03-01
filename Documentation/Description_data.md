Comments on how data is stored
============

# Main analysis

### all_cells_latencies.csv

* Contains: info on position and response latency for every cell. Used for space-related plots. Contains all cells, before any network-based subsetting (as used for network analysis; see below).
* Produced by:
* Used by:

Columns:

1. ibrain - experiment number
1. stage - developmental stage
1. cellid - cell number (arbitrarily assigned at quantification stage), used as id
1. x - cell position, x
1. y - cell position, y
1. dist - distance from the "most optimal" retinotopic center for the looming response (see methods in the main text for the algorithm)
1. lat - latency
1. ac - amplitude looming
1. af - same for flash
1. as - same for scramble
1. selfc - selectivity for Looming stimulus, compared to Flash
1. selsc - same for Looming compared to Scramble

### avamps_allcells_allbrains stable.csv

* Contains: response amplitudes for every cell recorded in every experiment.

TODO Q: Most probably  truly all cells (before we restrict analysis to the largest component), or is it post-restriction?

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
* Produced by: manually, by combining different console outputs of graph_structure_analyzer.m into one meaningful table
* Used by: caimg_experimens_assorted_figures.R, caimg_network_analysis.R

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
1. selassort - weighted assortativity for cell selectivity (whether similarly selective cells tended to be connected to each other)

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
* Produced by: graph_structure_analyzer.m ; only if doSelectivity flag is set to TRUE.
* Used by: caimg_centrality_analysis.r

Columns: 

1. nexp - experiment number
1. ncell - cell number in the original list (not all numbers are present, as we are only looking at the largest weakly connected component here)
1. selFC - Flash-Looming selectivity (default for most analyses)
1. selSC - Scrambled-Looming selectivity
1. indegree - in-degree of this node (cell)
1. outdegree - out-degree of this node
1. katz - Katz centrality of this node: a measure of information sink
1. revkaz - Katz centrality calculated on a reversed graph: a feature of information (activation) sources
1. clust - local clustering coefficient for this node
1. spiking - average amplitude of spiking response in this cell
1. insel - average selectivity of all cells connected to this cell (in-connections)

First 14 experiments are stage 46; the reast are stage 49.


# Original dataset info

For each experiment, we originally had the following set of data files:

1.	Excel file with fluorescence signals
2.	TIFF containing ROIs – Not necessary, as ROI (x,y) coordinates are also stored in same Excel data file.
3.	A print-screen of ROIs (for further processing in imageJ)
4.	.mat file with a Matlab structure (S) containing all data from this Excel file, and also spikes reconstructed
