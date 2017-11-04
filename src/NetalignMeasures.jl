__precompile__()

module NetalignMeasures

using MatrixNetworks
using Distances
using StatsBase
using Distributions
import Distributions: dim

export NetalignMeasure, NetalignScore, measure, score, dim

"""
    abstract type NetalignMeasure end
    abstract type NetalignScore end

Sub-types of `NetalignMeasure` are network alignment measures.
A type instance will contain information such that given the
type instance and an alignment, we are able to calculate
the alignment score.

Sub-types of `NetalignScore` contain network alignment score information.
The type instance will contain information of a particular network alignment measure and
an alignment. The `score` field in the type instance will contain the network alignment score.
`score` measures network similarity and will be between `0` and `1`.

### Two networks

Namely, let `S3Measure <: NetalignMeasure` be the S3 network alignment measure.
It will be instantiated as
    meas = S3Measure(G1,G2)
where `G1` and `G2` are the adjacency matrices of the two relevant
networks such that `n1 <= n2`, where `n1` and `n2` are the number of nodes
in `G1` and `G2`, respectively.

Then, given an alignment (described using a vector of `Int`s, `f`, such that when
node `i` in graph `G1` is mapped to node `j` in graph `G2`, `f[i] = j`),
we can calculate the alignment score using the `measure` and `score` functions.

Given a type `NAM <: NetalignMeasure`, the interface functions such a type
under `NetalignMeasure` include:

    measure(meas::NAM, f::Vector{Int}) -> NAS
Returns a struct of type `NAS <: NetalignScore`, which corresponds to type `NAM`,
containing information regarding the alignment score
of `f` with respect the alignment measure `NAM`.

Namely, given an instance of `NAM`, for example, `meas = NAM(G1,G2)`, then
    x = measure(meas,f)
`x.score` will contain the alignment score.
`x` might also contain other information such as the numerator
and denominator of the alignment score, or the number of conserved
and non-conserved edges, etc.

    score(meas::NAM, f::Vector{Int}) -> Float64
Returns the alignment score of `f` with respect the alignment measure `NAM`.

### Multiple networks

Not all `NetalignMeasure`s work for more than two networks. One that does it the
LCCS (largest common conserved sub-graph) measure.
This is an example with 3 networks but it is easily extended to `k` networks.
Namely, let `LCCSMeasure <: NetalignMeasure` be the LCCS network alignment measure.
It will be instantiated as
    meas = LCCSMeasure(G1,G2,G3)
where `G1,G2,G3` are the adjacency matrices of the four relevant
networks such that `n1 <= n2 <= n3`, where `n1,n2,n3` are the number of nodes
in `G1,G2,G3`, respectively.

Then, given an alignment we can calculate the alignment score using the `measure` and `score` functions.
An alignment between `k` networks, `G[1],...,G[k]`, is described using a tuple, `fs`,
of `k-1` vectors of `Int`s.
For a given vector of `Int`s, `fs[l]`, fs[l] is 1-1 mapping from graph `G[l]` to graph `G[k]`.
And similar to the two network case, each `fs[l]` is such that when
node `i` in graph `G[l]` is mapped to node `j` in graph `G[k]`, `fs[l][i] = j`).

Given a type `NAM <: NetalignMeasure`, the interface functions such a type
under `NetalignMeasure` include:

    measure(meas::NAM, fs::Tuple{Vector{Int}}) -> NAS
    measure(meas::NAM, fs::Vector{Int}...) -> NAS
Returns a struct of type `NAS <: NetalignScore`, which corresponds to type `NAM`,
containing information regarding the alignment score
of `fs` with respect the alignment measure `NAM`.

    score(meas::NAM, fs::Tuple{Vector{Int}}) -> Float64
    score(meas::NAM, f::Vector{Int}...) -> Float64
Returns the alignment score of `f` with respect the alignment measure `NAM`.

    score(x::NAS) -> Float64
Returns `x.score`.

    dim(meas::NAM) -> Int
Returns the number of nodes in the largest network among the relevant networks.
Note: this is always the last network since the networks
_must_ be ordered w.r.t. node size.

    dim : NAM x Int -> Int
Returns the number of nodes in the kth network.

### Other functions

    score(x::NAS) -> Float64
Returns `x.score`.

    dim(meas::NAM) -> Int
Returns the number of nodes in the largest network among the relevant networks.
Note: this is always the last network since the networks
_must_ be ordered w.r.t. node size.

    dim : NAM x Int -> Int
Returns the number of nodes in the kth network.

Note: The input adjacency matrices _must_ be symmetric.
"""
abstract type NetalignMeasure end
function dim end

abstract type NetalignScore end
score(x::NetalignScore) = x.score
score(m::NetalignMeasure,fs::Vector{<:Integer}...) = score(measure(m,fs...))

function measure end
measure(meas::NetalignMeasure, fs::Tuple{Vector{<:Integer}}) =
    error("measure function not implemented for $(typeof(meas))")
measure(meas::NetalignMeasure, fs::Vector{<:Integer}...) = measure(meas,fs)

include("mergesorted.jl")
include("networks.jl")
include("dynamicnetworks.jl")
include("measures.jl")
include("nodemeasures.jl")
include("aln.jl")

end
