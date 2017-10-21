__precompile__()

module NetalignMeasures

export NetalignMeasure, NetalignScore, measure, score, dim

"""
Types under NetalignMeasure describes a network alignment measure.
A type structure will contain information such that given the
type structure and an alignment, we are able to calculate
the alignment score.

So, a typical measure NAM <: NetalignMeasure
will be instantiated as
    meas = NAM(G1,G2)
where G1 and G2 are the adjacency matrices of
two relevant networks such that n1 <= n2,
where n1 and n2 are the number of nodes in G1 and G2 respectively

Then, given an alignment (described using a permutation,
specifically a vector of Ints of length n2),
we can calculate the alignment score using the score function.

The interface functions to types under NetalignMeasure include:

 measure : NAM x f -> Float
Returns object containing information regarding the alignment score
    of f with respect the alignment measure NAM.
Given
    x = measure(NAM,f)
x.score will contain the alignment score.
x might also contain other information such as the numerator
and denominator of the alignment score, or the number of conserved
and non-conserved edges, etc.

 dim : NAM -> Int
Returns the number of nodes in the largest network
Note: this is always the last network since the networks
_must_ be ordered wrt node size.

 dim : NAM x Int -> Int
Returns the number of nodes in the kth network.

Note: The input adjacency matrices _must_ be symmetric.
"""
abstract type NetalignMeasure end
function dim end

abstract type NetalignScore end
score(x::NetalignScore) = x.score

function measure end

include("mergesorted.jl")
include("networks.jl")
include("dynamicnetworks.jl")
include("measures.jl")

end
