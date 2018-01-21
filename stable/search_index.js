var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures documentation",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#NetalignMeasures-documentation-1",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures documentation",
    "category": "section",
    "text": "Network alignment measures including S3, DS3, WEC, DWEC, as well as node similarities for network alignment including GDV similarity, degree similarity, etc."
},

{
    "location": "index.html#Installation-1",
    "page": "NetalignMeasures documentation",
    "title": "Installation",
    "category": "section",
    "text": "NetalignMeasures can be installed as follows.Pkg.clone(\"https://github.com/vvjn/NetalignMeasures.jl\")"
},

{
    "location": "index.html#NetalignMeasures.NetalignMeasure",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.NetalignMeasure",
    "category": "Type",
    "text": "abstract type NetalignMeasure end\nabstract type NetalignScore end\n\nSub-types of NetalignMeasure are network alignment measures. A type instance will contain information such that given the type instance and an alignment, we are able to calculate the alignment score.\n\nSub-types of NetalignScore contain network alignment score information. The type instance will contain information of a particular network alignment measure and an alignment. The score field in the type instance will contain the network alignment score. score measures network similarity and will be between 0 and 1.\n\nTwo networks\n\nNamely, let S3Measure <: NetalignMeasure be the S3 network alignment measure. It will be instantiated as     meas = S3Measure(G1,G2) where G1 and G2 are the adjacency matrices of the two relevant networks such that n1 <= n2, where n1 and n2 are the number of nodes in G1 and G2, respectively.\n\nThen, given an alignment (described using a vector of Ints, f, such that when node i in graph G1 is mapped to node j in graph G2, f[i] = j), we can calculate the alignment score using the measure and score functions.\n\nGiven a type NAM <: NetalignMeasure, the interface functions such a type under NetalignMeasure include:\n\nmeasure(meas::NAM, f::Vector{Int}) -> NAS\n\nReturns a struct of type NAS <: NetalignScore, which corresponds to type NAM, containing information regarding the alignment score of f with respect the alignment measure NAM.\n\nNamely, given an instance of NAM, for example, meas = NAM(G1,G2), then     x = measure(meas,f) x.score will contain the alignment score. x might also contain other information such as the numerator and denominator of the alignment score, or the number of conserved and non-conserved edges, etc.\n\nscore(meas::NAM, f::Vector{Int}) -> Float64\n\nReturns the alignment score of f with respect the alignment measure NAM.\n\nMultiple networks\n\nNot all NetalignMeasures work for more than two networks. One that does it the LCCS (largest common conserved sub-graph) measure. This is an example with 3 networks but it is easily extended to k networks. Namely, let LCCSMeasure <: NetalignMeasure be the LCCS network alignment measure. It will be instantiated as     meas = LCCSMeasure(G1,G2,G3) where G1,G2,G3 are the adjacency matrices of the four relevant networks such that n1 <= n2 <= n3, where n1,n2,n3 are the number of nodes in G1,G2,G3, respectively.\n\nThen, given an alignment we can calculate the alignment score using the measure and score functions. An alignment between k networks, G[1],...,G[k], is described using a tuple, fs, of k-1 vectors of Ints. For a given vector of Ints, fs[l], fs[l] is 1-1 mapping from graph G[l] to graph G[k]. And similar to the two network case, each fs[l] is such that when node i in graph G[l] is mapped to node j in graph G[k], fs[l][i] = j).\n\nGiven a type NAM <: NetalignMeasure, the interface functions such a type under NetalignMeasure include:\n\nmeasure(meas::NAM, fs::Tuple{Vector{Int}}) -> NAS\nmeasure(meas::NAM, fs::Vector{Int}...) -> NAS\n\nReturns a struct of type NAS <: NetalignScore, which corresponds to type NAM, containing information regarding the alignment score of fs with respect the alignment measure NAM.\n\nscore(meas::NAM, fs::Tuple{Vector{Int}}) -> Float64\nscore(meas::NAM, f::Vector{Int}...) -> Float64\n\nReturns the alignment score of f with respect the alignment measure NAM.\n\nscore(x::NAS) -> Float64\n\nReturns x.score.\n\ndim(meas::NAM) -> Int\n\nReturns the number of nodes in the largest network among the relevant networks. Note: this is always the last network since the networks _must_ be ordered w.r.t. node size.\n\ndim : NAM x Int -> Int\n\nReturns the number of nodes in the kth network.\n\nOther functions\n\nscore(x::NAS) -> Float64\n\nReturns x.score.\n\ndim(meas::NAM) -> Int\n\nReturns the number of nodes in the largest network among the relevant networks. Note: this is always the last network since the networks _must_ be ordered w.r.t. node size.\n\ndim : NAM x Int -> Int\n\nReturns the number of nodes in the kth network.\n\nNote: The input adjacency matrices _must_ be symmetric.\n\n\n\n"
},

{
    "location": "index.html#Overview-1",
    "page": "NetalignMeasures documentation",
    "title": "Overview",
    "category": "section",
    "text": "CurrentModule = NetalignMeasuresNetalignMeasure"
},

{
    "location": "index.html#NetalignMeasures.S3Measure",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.S3Measure",
    "category": "Type",
    "text": "S3Measure(G1::SparseMatrixCSC,G2::SparseMatrixCSC)\n\nS3 edge conservation measure.\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.DS3Measure",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.DS3Measure",
    "category": "Type",
    "text": "DS3Measure(G1::SparseMatrixCSC,G2::SparseMatrixCSC)\n\nDS3 event conservation measure.\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.WECMeasure",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.WECMeasure",
    "category": "Type",
    "text": "WECMeasure(G1::SparseMatrixCSC,G2::SparseMatrixCSC,S::AbstractMatrix)\n\nWeighted edge conservation from WAVE.\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.DWECMeasure",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.DWECMeasure",
    "category": "Type",
    "text": "DWECMeasure(G1::SparseMatrixCSC{Events},G2::SparseMatrixCSC{Events},S::AbstractMatrix)\n\nDynamic weighted edge conservation from DynaWAVE.\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.LCCSMeasure",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.LCCSMeasure",
    "category": "Type",
    "text": "LCCSMeasure(G::SparseMatrixCSC...)\n\nLargest common conserved sub-graph, from multiMAGNA++ paper. Works for >2 networks.\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.GhostMeasure",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.GhostMeasure",
    "category": "Type",
    "text": "GHOST signature and node similarity\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.ConvexCombMeasure",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.ConvexCombMeasure",
    "category": "Type",
    "text": "ConvexCombMeasure(S::A,T::B,alpha::Float64) where {A<:NetalignMeasure,B<:NetalignMeasure}\n\nConvex combination of two NetalignMeasures. 0 <= alpha <= 1.\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.NodeSimMeasure",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.NodeSimMeasure",
    "category": "Type",
    "text": "Node conservation measure\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.NullMeasure",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.NullMeasure",
    "category": "Type",
    "text": "Null measure Always returns 0.0 score\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.NetalMeasure",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.NetalMeasure",
    "category": "Type",
    "text": "NetalMeasure(G1::SparseMatrixCSC,G2::SparseMatrixCSC)\n\nMeasure from NETAL paper.\n\n\n\n"
},

{
    "location": "index.html#Measures-1",
    "page": "NetalignMeasures documentation",
    "title": "Measures",
    "category": "section",
    "text": "S3Measure\nDS3Measure\nWECMeasure\nDWECMeasure\nLCCSMeasure\nGhostMeasure\nConvexCombMeasure\nNodeSimMeasure\nNullMeasure\nNetalMeasure\n"
},

{
    "location": "index.html#Dynamic-Networks-1",
    "page": "NetalignMeasures documentation",
    "title": "Dynamic Networks",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#NetalignMeasures.Events",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.Events",
    "category": "Type",
    "text": "An dynamic edge in a dynamic network is defined as a set of events. An event is an interaction between two nodes from start time t0 to stop time t1, t0 <= t1, and is represent as a tuple (t0,t1). The set of events is represented as a Vector, sorted wrt start time t0, and if two start times are equal, then stop time t1.  Note: Unless specified, it is assumed that two events in an event set will not overlap, i.e. there exists no two tuples (t0,t1) and (s0,s1) s.t. intersect([t0,t1],[s0,s1]) != null\n\nA dynamic network can be represented using a SparseMatrixCSC{Events} structure where an element in the matrix (i.e. an edge) is represented using the following structure\n\n\n\n"
},

{
    "location": "index.html#Types-1",
    "page": "NetalignMeasures documentation",
    "title": "Types",
    "category": "section",
    "text": "Events"
},

{
    "location": "index.html#NetalignMeasures.cet_ncet",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.cet_ncet",
    "category": "Function",
    "text": "Calculates CET and NCET of the first n elements in a sorted array of events As in, the events from two node pairs are combined together and tis is the combined events Sortedness is not checked\n\n\n\nCalculates CET and NCET of the events in tsl and tsr Assumes the events in each array are non-overlapping Sortedness is not checked\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.networkactivity",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.networkactivity",
    "category": "Function",
    "text": "Sum of active duration over all events in G\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.snapshots",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.snapshots",
    "category": "Function",
    "text": "Convert to snapshots with window size t_w\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.nodeactivity",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.nodeactivity",
    "category": "Function",
    "text": "Total time during which events in event set are active\ni.e. active duration of an edge\n\n\n\nInput: adj. matrix of a dynamic network Output: adj. matrix of static network containing corresponding  duration of each edge in the dynamic network\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.fixevents",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.fixevents",
    "category": "Function",
    "text": "Given a vector of timestamps (sorted wrt the start time), with possibly overlapping durations, return a vector of timestamps with no overlapping durations\n\n\n\nGiven two vectors of timestamps (sorted wrt the start time), with possibly overlapping durations, return a vector of timestamps with no overlapping durations\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.fixevents!",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.fixevents!",
    "category": "Function",
    "text": "Given a vector of timestamps (sorted wrt the start time), with possibly overlapping durations, modify vector of timestamps s.t. it has no overlapping durations\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.snapshot",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.snapshot",
    "category": "Function",
    "text": "If event intersects with [t_s,t_e] then add edges to snapshot\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.widenevents!",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.widenevents!",
    "category": "Function",
    "text": "Make each event in network wider by adding making each event occur pre time ahead and post time later\n\n\n\n"
},

{
    "location": "index.html#Functions-1",
    "page": "NetalignMeasures documentation",
    "title": "Functions",
    "category": "section",
    "text": "cet_ncet\nnetworkactivity\nsnapshots\nnodeactivity\nfixevents\nfixevents!\nsnapshot\nwidenevents!"
},

{
    "location": "index.html#Static-Networks-1",
    "page": "NetalignMeasures documentation",
    "title": "Static Networks",
    "category": "section",
    "text": ""
},

{
    "location": "index.html#NetalignMeasures.degree",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.degree",
    "category": "Function",
    "text": "degree vector of an adjacency graph \n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.dependency",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.dependency",
    "category": "Function",
    "text": "probability that any of i's neighbors will be conserved (netal)\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.expci",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.expci",
    "category": "Function",
    "text": "approx expected number of conserved interactions (netal)\n\n\n\nexpci(..)[i,j] = I[i,j] = expected # of conserved edges if i is aligned to j\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.adjnodes",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.adjnodes",
    "category": "Function",
    "text": "Given an adj. matrix of a graph, return index of nodes  adjacent to node i\n\n\n\n"
},

{
    "location": "index.html#NetalignMeasures.adjnzval",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.adjnzval",
    "category": "Function",
    "text": "Given an adj. matrix of a graph, return values of nodes  adjacent to node i\n\n\n\n"
},

{
    "location": "index.html#Functions-2",
    "page": "NetalignMeasures documentation",
    "title": "Functions",
    "category": "section",
    "text": "degree\ndependency\nexpci\nadjnodes\nadjnzval"
},

{
    "location": "index.html#NetalignMeasures.pca",
    "page": "NetalignMeasures documentation",
    "title": "NetalignMeasures.pca",
    "category": "Function",
    "text": "Input: X : m x n matrix Each column of X is a feature vector Get the first dimensions containing k fraction of variance\n\nOutput:     PCA object     with mean, projection matrix, and variances     U is projection matrix     projection is T = U' * X     U*T =~ X, where X is mean centered\n\n\n\n"
},

{
    "location": "index.html#Other-functions-1",
    "page": "NetalignMeasures documentation",
    "title": "Other functions",
    "category": "section",
    "text": "pca"
},

]}
