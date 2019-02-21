# caimg_network_measures

Data file description

Columns: 

1. type - one of three groups of values:
  * 1actual - actual measurements on a reconstructed graph. Numbers before the name are so that R would sort them properly by default.
  * 2rewire - an average of the same value, across 100 degree-preserving rewires of the reconstructed graph
  * 3random - an average of the same value, across 100 complete random rewires of the reconstructed graph (but with weights of individual edges still preserved)
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