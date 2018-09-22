Ca imaging project – Master File
================================

#Protocols

File: ca protocol calculations.xls for a comparison of concentrations that different people from different labs use in their experiments.

## Early processing, comments

### Intensity measurement

### Data export and transformation

Beware: by default NES exports to Excel not full data, but only a brief report on data, which looks like a nice square table, but is limited to 74 first ROIs. To get all data, switch export model to "full data", in which case you won't get a square table, but rather information for each ROI given one under another. An unfortunate surprise you should be prepared for is that in "Report" mode time axis is given in seconds, while in "Full data" mode it is given in ms. I spent about 3 hours troubleshooting my program, trying to understand why it works on old format data, but does not work on new, before I realized that the time axis is all wrong.

# Processing programs

## Main pipeline, parts that should not be repeated

### caimaging_read

Obsolete; see caimaging_read2

### caimaging_read2

Reads all original Excel files (should be explicitly coded in the file). Saves them as mat files. The first step in processing. This is a fairly slow process.

### caimaging_basic_batch.m

Loads data in mat form, applies non-negative deconvolution (the one by Vogelstein), and saves the data back again, but now with the spike-vector in it. All math is performed by the external module caimaging_basic (see below). Input files are hard-coded, and need to be commented/uncommented manually to make the thing work. Also relies on caimaging_find_good_cell to find a strong cell; this cell number is then passed along to the caimaging_basic processor, and is used to estimate deconvolution parameters. This is a fairly slow process as well.

###caimaging_basic.m

Receives data in stricture form (S) as an input, applies [Vogelstein] fast_oopsi.m deconvolution, and updates the S structure.

###caimaging_find_good_cell

Automatically finds a cell with strong signal. In practice, removes biases, looks at 10 cells with biggest time-variances (strongest signals), and then returns the number of the one in the middle. So it should be 5th strongest cell in the set.

###caimaging_process.m

Appears to be abandoned and obsolete. First attempt to build a comprehensive data analyzer. Was superseded by caimaging_structure.

###caimaging_structure.m

Batch connectivity reconstruction routine. Goes through mat files one by one (file names are hard-coded in the header of the program), and calls transfer entropy on them; then saves the result in a separate file (different folder). Both in/out folders and file names are hard-coded in the body. This process is very slow, and takes somewhere between 10 and 24 hours, for the full set of data (50 brains). TE calculation can be disabled by the doTE flag (row 80). Has some logic on when to update and when not to update the data, so that it would not ruin expensive TE calculations by mistake. Also calculates correlations using reshuffle_corr, averages the data in several different ways and concatenates it together. Saves all results in a structure that is then passed to graph_structure_analyzer.

Formally it calculates selectivity for every cell, but in fact this selectivity was calculated in a very simplistic way that is now considered obsolete. So it needs to be recalculated (by caimaging_pca described below), and included into any analyses of graph from there; NOT from the file that contains the graph.

Also builds a "mugshot" picture of each dataset, based on the timing of spiking, and PCA analysis. This functionality is deprecated however, and moved to caimaging_pca.

The routine returns all matrices (te, corr, p) in the "direct", unflipped format, so w_ij is a connection from i to j. It means that G(w) makes sense, but to run a process on this matrix you'd have to flip it: s = w'*s.

##Main pipeline, late stages of analysis that can be repeated

###caimaging_pca

Despite the name, the logic here is that it contains all analysis that doesn't involve the connectivity reconstruction graph. Selectivity, PCA, nice mugshot, a hunt for cell types, activation trajectories. Hard-codes which brains are read, so this should be changed manually by commenting and uncommenting individual rows. Contains a bunch of flags and switches (around row 70) that turn different types of analysis and figures on and off.

###graph_structure_analyzer

Graph analysis. Relies on the data files created by caimaging_structure (file names and data folders are hard-coded in the header of the program). For each file then calls an internal processing routine (first function, located in the same m-file, below the main function). It means that all constants and switches that turn different parts of the code on and off are defined in the 1st subroutine (approx. row 70 and below), not in the main function header. Most analysis blocks can be turned on and on with flags (constants), but some also contain some lazy hard-coded switches (if(1) or if(0)), which is particularly true for troubleshooting figures. The main output goes into console with fprintf, which means that in most cases two functional blocks cannot be active at the same time, as it will make lines intercolate. Better get one table first, then another table later, etc.

##Important dependencies 

###Borrowed dependencies: 

several subroutines (mostly network and centrality measures) from this set of tools:
https://sites.google.com/site/bctnet/measures/list 

###transfer_entropy

Calculates transfer entropy. Gets its main inputs as parameters; the only part that is hard-coded is whether the histogram is based on quantiles (default), or on linear binning. Called by caimaging_process. Returns w as a "straight w", where w_ij goes from i to j.

Includes a built-in tester, which also calls and tests reshuffle_corr, by courtesy.

### reshuffle_corr

Calculates reshuffle-adjusted correlations, similarly to how we do it for TE. Called by caimaging_process. By default calculates "connectivity reconstruction correlations" that shift one of the traces by one frame. For "shifted" correlations, returns a "straight w" with w_ij going from i to j. Performs 20 reshuffles (hard-coded in the header).

### spectralClustering

Performs adaptive spectral clustering, by Ng Jordan Weiss algorithms, followed by adaptive k-means clustering with elbow method. Tests itself.

### network_rewire

Performs degree-preserving rewiring. Get an adjacency matrix and N shuffles, returns a new adjacency matrix.

### selectivity_graph

Draws the connectivity graph, with nodes colored according to some feature (such as selectivity), both in original coordinates and optimized coordinates. Is not heavily used, but seems to be up to date.

### myCentrality

A collection of different centrality measures I coded, and some graph statistics. Take the adjacency matrix, returns either a column of centrality values, or one characteristic measure for the entire graph (open the program, and see what options are there).

### myCyclycity

My attempts to build a useful cyclicity measure that for now failed. So it may not be included in the final product; we'll see.

### kmeans_opt

Performs optimal k-means clustering using the elbow method (looks at the point when the improvement in variance explained becomes worse). It was not my program originally (it has attribution in the program file itself), but I changed the way it reports the share of explained variance.

## Helper programs

### caimaging_extra_figures1

Buiilds all kinds of hard-coded figures for the paper.

### caimaging_graph_given_degree
Tries to build a graph with a given degree distribution, for illustration purposes. Doesn't always do a good job. Can work both on directed and undirected graphs.

### graph_structure_tester

A toy model (build a graph manually, then gradually destroy it) that can be used to assess the validity of different network measures. 
Abandoned and untested (potentially obsolete) utilities

### network_flow

Obsolete. Used to calculate the network flow (those "new" measures I introduced). Could be called from the analysis program, or from the model. Includes a built-in tester. Now however these measures are re-implemented in myCentrality and maybe, to some extent, myCyclicity, so this realization is abandoned.

### caimaging_mat2csv_batch

Technical file that reads mat files and dumps them into csv. Inputs are hard-coded. Is only used when data needs to be passed to students. Before it can be reused, needs to be seriously revisited, as I'm not sure what data it's even saving (is it only the spike-series?).

### caimaging_browser

Interactive browser of cells; can be called by passing the data structure S. Is based on my general script scatterbrowser. Currently is not used.

### caimaging_draw_spikes

Not used; abandoned.

## Model

### model_stdp_multisens_1

Main program to generate data. Set the constants and flags in the header, then run it from the console in a cycle. It would generate a bunch of files, name them according to the YYMMDDhhmmss convention, and save them in a certain folder.

### model_stdp_tester

Main analysis program. Reads files from the model results folder, and analyzes them one by one. Calculates lots of different values, and tried to plot them all, as well as save them all in a separate file. The file name changes every day, so from day to day it won't overwrite the previous results, but WITHIN a day it WOULD overwrite the results, so be careful and rename the output files if you are analyzing stuff in batches. The file name is "modelAnalysisYYMMDD", in the "Modeling" folder.

May be fairly slow, if there are lots of datasets to process.

### model_stdp_curve_plotter

Reads the curves of how different network measures changed in time, as the model_stdp_multisens_1 program ran; from the model result files created by this program. Then plots them all on the same plot, creating a draft figure that can then be ported to Illustrator. The way it reads the data is very similar to model_stdp_tester, as the program was forked from it when the functionalities diverged.

### shotgun_testing

Creates a bunch of modular networks, then gradually destroys them, every time randomly subsampling, and comparing network measurements on a full network to that on a subsampling network. Sadly, some most interesting measures break down on networks of about 5000 nodes, so I cannot check realistic numbers for the tectum. Also it seems to be quite dependent on the structure of the underlying network (for example, a dense torus of random modules seems to have too much flow, so my "hierarchy" measure doesn't change as it is randomized. (or maybe the measure is wrong?)

Depends on create_modular_network

### create_modular_network

Creates a bunch of non-random networks.

## Data files

A structure of the original dataset:

1.	Excel file with fluorescence signals
2.	TIFF containing ROIs – Not necessary, as ROI (x,y) coordinates are also stored in same Excel data file.
3.	A print-screen of ROIs (for further processing in imageJ)
4.	.mat file with a Matlab structure (S) containing all data from this Excel file, and also spikes reconstructed

### 131213 Latency calibration

Old-school (ephys paper –compatible) stimuli, recorded by the camera, and compared. Stored in the main folder, not in respective Data folder.
 
Based on this data (average lightness of the fiber, negative thresholdoing at a 99.9% value), average latencies compared to a flash (in ms) are: [0 225 172 468] for flash, crash, grid and back respectively. The reason for latency adjustment for "grid" and "crash" is that for these stimuli a one-pixel-size object (in the very beginning of the stimulus) would not actually be shown on the screen, even though the program would try to produce it (due to computer->VGA->NTSC->LCD conversion). Therefore actual stimulus is shown a tiny bit later than Matlab starts generating it. For "back" stimulus the situation is even more dire, as the stimulus starts off-limits of the fiber, and only some time later enters the fiber image.

Note that this approach differs from what was used in "Tectal processing" paper, as there we assumed "flash", "grid", "ramp" and "crash" having equal latencies. Check out the file "crampf - Eye recordings.xls", tab "Stimuli" for comparison.

I think it would make sense to introduce discontinuity with the previous paper (ephys), and re-adjust stimuli, so that they start and end simultaneously (otherwise dynamics comparisons will be really problematic). Some pilot experiments however (131209 and 131213) still use old stimuli, and thus require a heavy latency adjustments.