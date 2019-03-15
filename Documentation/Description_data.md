Comments on how data is stored
============

# Main analysis

### all_cells_latencies.csv

* Contains: info on position and response latency for every cell. Used for space-related plots. Contains all cells, before any network-based subsetting (as used for network analysis; see below).
* Produced by: caimaging_pca.m
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

* Contains: response amplitudes for every cell recorded in every experiment. (Small amplitude cells aren't removed)
* Produced by: caimaging_pca.m

### sel_allcells_allbrains stable.csv

* Contains: selectivity measurements for every cell. (Small amplitude cells aren't removed)
* Produced by: caimaging_pca.m
* Used by:

Columns:

1. fc - Flash-Looming selectivity
1. fs - Flash-Scrambled selectivity
1. sc - Scrambled-Looming selectivity
1. ibrain - experiment number
1. stage - developmental stage for this experiment

### caimg_experiment_summary.csv

* Contains: Basic low-level information about each experiment and each reconstructed network. This file contains 30 rows; one row for each experiment.
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
1. rLatSel - Pearson correlation coefficient between cell latency and its selectivity, for this experiment
1. rLatDist - Pearson correlation coeff between cell latency and distance from the estimated "retinotopy center" for looming stimuli
1. rSelDist - Pearson correlation coeff between cell FL selectivity and distance from the "retintopy center"
1. nEns - number of detected ensembles
1. maxModul - maximal modularity achieved (an estimation of how well the ensembles were resolved)
1. clustCompact - a measure of ensemble spatial locality (see text)
1. clustPref - a measure of preferential connections between cells within the same ensemble (see text)


### caimg_network_measures.csv

* Contains: global network measures (Figure 4). The biggest difference between this file and the previous one `caimg_experiment_summary` is that here each experimenet is represented by 2 different rows (so 60 rows in total): one row for normal analysis, and another row with all same measurements, but from a randomly rewired graph. That's why this file and `caimg_experiment_summary` cannot be combined in one.
* Produced by: manually, by adding columns of interest
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

### Experiments used in figures

* Figure 1: visual field: from an early multisensory experiment 140128 (stage 49) that wasn't included in analysis
* Figure 1 signals: 140722 - stage 49
* Figure 2 selectivity map, selectivity correlations: 140718b
* Figure 3 components: not sure, but not 140709
* Figure 3: PCA field, correlation matrix, modularity plot, ensemble field: 140709 - stage 49.
* FIgure 4 Corr, TE, and W matrices: 140722 - stage 49
* Figure 4 connectivity graph: not sure

# Model analysis

### Model Analysis summaries

Files:

* modelAnalysis181020 1 slide looming.csv - main set of experiments
* modelAnalysis181020 2 slide vis.csv - experiments in which looming stimuli were replaced with random visual transitions (that include looms as a subset)
* modelAnalysis181020 3 slide rand.csv - experiments in which the model was exposed to random noise
* modelAnalysis181023 4 slide looming nointrinsic.csv - looming stimuli, but weak intrinsic plasticity (see Methods)
* modelAnalysis181023 5 slide looming Hebb.csv - looming stimuli, but Hebbian plasticity instead of STDP
* modelAnalysis181023 6 decay looming.csv - looming stimuli, but synaptic decay instead of synaptic competition

Each file contains analysis results for 50 modeling experiments. For each experiments, we have values from 5 developmental points, taken at equal spaces across simulation.

Each file contains the following columns:

1. file - name of the matching file that contains simulation results (technical column)
1. type - experiment type (should be the same within each csv file)
1. competition - competition type: either "slide" (default; synaptic competition) or "decay" (synaptic decay)
1. stage - developmental stage, from 1 (naive network, one step after initialization), to 5 (final network)
1. rewire - whether the network was randomly rewired before analysis (should be "original")
1. fullBrainSel - full brain selectivity to Loom compared to Flash (first responses of all cells are added, then selectivity is calculated)
1. meanSel - first selectivity of each cell is calculated, then selectivity values are averaged
1. shareSelCells - share of cells that respond to looming more strangly than to flash
1. sel90perc - Flash-Looming selectivity of top 90% most selective cells
1. sel90m50 - a difference between top 90% selectivity (above) and median selectivity
1. bestPredict - best possible prediction of stimulus identity achieved from cell response amplitudes (see methods)
1. fullBrainSel_SC - same as fullBrainSel above, but for Scramble-Looming rather than Flash-looming selectivity
1. meanSel_SC - ditto
1. shareSelCells_SC - ditto
1. sel90perc_SC - ditto
1. sel90m50_SC - ditto
1. rPosSel - correlation between node position (distance from the center) and node selecttivity
1. rDirWei - correlation nbetween edge direction off-center, and its weight
1. mDistWei - a ratio of weighted average distance between cells (weighted by edge weight) and average distance between all cells. The smaller the value, the more local are the connections.
1. rSelClu - correlation between node selectivity and its local clustering coefficient
1. rCluSpk - correlation between node clustering and its activity (average spiking)
1. rSelNet - correlation between node selectivity and its Katz rank
1. rSelRnt - same, but for Katz rank on a reveersed graph
1. rSelGth - correlation between node selectivity and an adjusted ratio of (1+ weighted in-degree) devided by (1+ weighted out-degree) (abandoned analysis)
1. rSelIns - correlation between node selectivity and its in-degree
1. selAssort - assortativity of node selectivities
1. shESelGrow - share of edges that increase selectivity (sel grows along the edge)
1. selEGrowth - weighted average selectivity growth across the edge (weighted by edge weight)
1. gammaIn - power of in-degree distribution
1. gammaOu - power of out-degree distribution
1. deg0 - frequency of nodes with 0 degree
1. deg12 - frequency of nodes with either degree 1 or 2
1. deg5p - frequency of nodes with degree of 5 or higher
1. nPCAto80 - number of PCA components to achieve 80% reconstructon of total response variability (total response of each cell, analyzed across all stimuli, of all 3 types)
1. nRichTo80F - same, but only for a subset of responses to "Flashes". Quantifies intrinsic network variability, rather than its representation of different stimuli
1. nRichTo80C - same, but for looming stimuli
1. nClusters - number of response clusters (modularity-based, see Methods)
1. clusterPreference - preference of connectedness within a cluster; see Methods
1. clusterCompactness - geographical compactness of clusters (quantification, see methods)
1. clustPrefPval - a p-value associated with cluster preference (whether weights within a cluster were statistically different from weights between clusters)
1. eff - network efficiency
1. modul - network modularity
1. clust - global clustering coefficient
1. flow - global hierarchical flow
1. revFlow - same, on a reversed graph
1. cycl - measure of graph cyclicity (abandoned analysis)
1. recip - share of reciprocal connections
1. rSelSpk - correlation of selectivity and sliking
1. rSelfcSelfs - correlation between Flash-Loom and Flash-Scrambled selectivities
1. rSelfcSelsc - correlation between Flash-Loom and Scrambled-Loom selectivities


# Original dataset info

Described here for historical purposes (not presented in this repo). For each experiment, we originally had the following set of data files:

1.	Excel file with fluorescence signals
2.	TIFF containing ROIs – Not necessary, as ROI (x,y) coordinates are also stored in same Excel data file.
3.	A print-screen of ROIs (for further processing in imageJ)
4.	.mat file with a Matlab structure (S) containing all data from this Excel file, and also spikes reconstructed
