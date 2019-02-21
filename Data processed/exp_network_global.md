# exp_network_global

Description of the file.

Columns: 

1. stage - tadpole stage (46 or 49)
1. name - experiment code (YYMMDD + a letter if more than one experiment in a given day)
1. type - one of three groups of values:
  1. actual - actual measurements on a reconstructed graph
  2. rewired - an average of the same value, across 100 degree-preserving rewires of the reconstructed graph
  3. erdos - an average of the same value, across 100 complete random rewires of the reconstructed graph (but with weights of individual edges still preserved)
1. eff - network efficiency
1. a_oi - Out-In assortativity
1. a_io - In-Out assortativity
1. a_oo - Out-Out assortativity
1. a_ii - In-In assortativity
1. modul - Spectral modularity
1. clust - Global clustering coefficient
1. flow - Hierarchical flow (see the text)