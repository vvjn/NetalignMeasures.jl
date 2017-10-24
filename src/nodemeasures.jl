export NodeSimMeasure, NodeSimScore

using MatrixNetworks
using Distances
using MultivariateStats

"""
Node conservation measure
"""
immutable NodeSimMeasure{T<:AbstractMatrix} <: NetalignMeasure
    S :: T
    function NodeSimMeasure{T}(S::T) where {T<:AbstractMatrix}
        if size(S,1) > size(S,2) error("Matrix dims") end
        new(S)
    end
end
NodeSimMeasure(S::T) where {T<:AbstractMatrix} =  NodeSimMeasure{T}(S)

immutable NodeSimScore <: NetalignScore
    score :: Float64
end

function measure(m::NodeSimMeasure,f::Vector{Int})
    x = 0.0
    S = m.S
    n = size(S,1)
    for i = 1:n
        x += S[i,f[i]]
    end
    NodeSimScore(x/n)
end

dim(m::NodeSimMeasure,d::Int) = size(m.S,d)
dim(m::NodeSimMeasure) = dim(m,2)

NodeSimMeasure(sym::Symbol, args...) =  NodeSimMeasure(Val{sym}(), args...)

"GDV similarity, from the original paper"
function NodeSimMeasure(::Val{:gdvs}, gdv1::AbstractMatrix{T},gdv2::AbstractMatrix{T}) where {T}
    weights = T[1.0, 0.838444533, 0.838444533, 0.838444533, 0.743940642, 0.676889065, 0.743940642, 0.743940642, 0.676889065, 0.743940642, 0.676889065, 0.676889065, 0.676889065, 0.676889065, 0.743940642, 0.676889065, 0.582385175, 0.624879821, 0.676889065, 0.624879821, 0.582385175, 0.582385175, 0.676889065, 0.676889065, 0.676889065, 0.624879821, 0.546456463, 0.676889065, 0.582385175, 0.582385175, 0.546456463, 0.676889065, 0.582385175, 0.582385175, 0.582385175, 0.624879821, 0.582385175, 0.546456463, 0.546456463, 0.624879821, 0.546456463, 0.582385175, 0.546456463, 0.582385175, 0.624879821, 0.624879821, 0.582385175, 0.515333598, 0.546456463, 0.582385175, 0.582385175, 0.515333598, 0.582385175, 0.487881285, 0.624879821, 0.582385175, 0.676889065, 0.582385175, 0.582385175, 0.546456463, 0.515333598, 0.582385175, 0.582385175, 0.515333598, 0.546456463, 0.582385175, 0.546456463, 0.546456463, 0.515333598, 0.624879821, 0.582385175, 0.582385175, 0.676889065]
    if size(gdv1,2) != size(gdv2,2) || size(gdv1,2) != length(weights) error("GDV sizes don't match") end
    @inline kernel(oai,obi,tai,tbi,wi) = wi * abs(oai-obi) / max(tai,tbi)
    D = zeros(T,size(gdv1,1),size(gdv2,1))
    olog1 = log.(gdv1 .+ one(T))
    tlog1 = log.(gdv1 .+ 2one(T));
    olog2 = log.(gdv2 .+ one(T))
    tlog2 = log.(gdv2 .+ 2one(T));    
    D = zeros(T,size(gdv1,1),size(gdv2,1))
    for i = 1:length(weights), v = 1:size(gdv2,1), u = 1:size(gdv1,1)
        @inbounds D[u,v] += kernel(olog1[u,i], olog2[v,i], tlog1[u,i], tlog2[v,i], weights[i])
    end
    D ./= sum(weights)
    D .= 1 .- D
    NodeSimMeasure(D)
end

"GDV similarity from the L-GRAAL paper"
function NodeSimMeasure(::Val{:lgraalgdvs}, gdv1::AbstractMatrix{T},gdv2::AbstractMatrix{T}) where {T}
    if size(gdv1,2) != size(gdv2,2) error("GDV sizes don't match") end
    S = zeros(T,size(gdv1,1),size(gdv2,1))
    @inbounds for i = 1:size(gdv1,2), v = 1:size(gdv2,1), u = 1:size(gdv1,1)
        a,b = minmax(gdv1[u,i], gdv2[v,i])
        S[u,v] += ifelse(b>0, a/b, 0.0)
    end
    S ./= size(gdv1,2)
    NodeSimMeasure(S)
end

"The normalized GDV similarity from the H-GRAAL paper"
function NodeSimMeasure(::Val{:hgraalgdvs}, G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                        gdv1::AbstractMatrix{T},gdv2::AbstractMatrix{T},alpha::T=0.5) where {T}
    gm = NodeSimMeasure(:gdvs, gdv1,gdv2)
    deg1 = vec(full(sum(x->T(x),G1,1)))
    deg2 = vec(full(sum(x->T(x),G2,1)))
    maxdeg = maximum(deg1) + maximum(deg2)
    S = gm.S
    @inbounds for v = 1:size(S,2), u = 1:size(S,1)
        S[u,v] = alpha * S[u,v] + (one(T) - alpha) * (deg1[u]+deg2[v]) / maxdeg
    end
    gm
end

"The GDV similarity normalized using PCA from Yuriy's dynamic graphlets paper"
function NodeSimMeasure(::Val{:pcagdvs}, gdv1::AbstractMatrix,gdv2::AbstractMatrix)
    X = hcat(gdv1',gdv2') # switch to make columns feature vectors    
    X = X[vec(var(X,2)) .>= 1e-8,:] # remove features w/ no variance
    X .-= mean(X,2)
    X ./= std(X,2) # standardize
    res = fit(PCA,X; pratio=0.99, mean=0)
    Y = sqrt.(res.prinvars + 1e-8) .\ res.proj' * X # project and whiten    
    R = 1 .- pairwise(CosineDist(), Y[:,1:size(gdv1,1)], Y[:,size(gdv1,1)+1:end]) ./ 2
    NodeSimMeasure(R)
end

"Relative similarity between v1 and v2: 1 - abs(a-b)/max(a,b)"
function NodeSimMeasure(::Val{:relative}, v1::AbstractVector,v2::AbstractVector)
    m = length(v1)
    n = length(v2)
    S = zeros(Float64,m,n)
    for i = 1:m, j = 1:n
        S[i,j] = 1 - abs(v1[i] - v2[j])/max(v1[i],v2[j])
    end
    NodeSimMeasure(S)
end

"Degree similarity"
NodeSimMeasure(::Val{:degree}, G1::SparseMatrixCSC,G2::SparseMatrixCSC) =
    NodeSimMeasure(:relative, vec(sum(G1,1)), vec(sum(G2,1)))

"Relative clustering coefficient"
NodeSimMeasure(::Val{:clustering}, G1::SparseMatrixCSC,G2::SparseMatrixCSC) =
    NodeSimMeasure(:relative, clustercoeffs(G1), clustercoeffs(G2))

"Eccentricity similarity"
function NodeSimMeasure(::Val{:eccentricity}, G1::SparseMatrixCSC,
                        G2::SparseMatrixCSC)
    D1 = hcat(map(i->dijkstra(G1,i)[1],1:size(G1,1))...)
    D2 = hcat(map(i->dijkstra(G2,i)[1],1:size(G2,1))...)
    NodeSimMeasure(:relative, maximum(D1,1), maximum(D2,1))
end

"Performs the migraal confidence transform"
function NodeSimMeasure(::Val{:migraal}, Xs::AbstractMatrix...)
    C = zeros(Float64,size(Xs[1]))
    for X in Xs
        for i = 1:size(X,1)
            for j = 1:size(X,2)
                sless = 0
                for k = 1:size(X,2)
                    sless += Int(X[i,k] < X[i,j])
                end
                C[i,j] = sless/size(X,2)
            end
        end
    end
    C ./= length(Xs)
    NodeSimMeasure(C)
end

"Transforms blast e-values to nodes similarities"
function NodeSimMeasure(::Val{:evalues}, E::AbstractMatrix, clamp=1e-101)
    # -log(e-value) transform
    B = -log.(E)
    # If -log(E-value) above -log(clampval), ie if E-value below clampval,
    # then clamp it
    # Otherwise we get issues with -log(E-value) being Inf or some
    # ridiculously huge number
    clamp!(B, 0, -log(clampval))
    # Normalize by dividing by maximum value
    B ./= maximum(B)
    NodeSimMeasure(B)    
end

"Transforms blast e-values to node similarities"
function NodeSimMeasure(::Val{:evalues}, E::SparseMatrixCSC, clamp=1e-101)
    B = copy(E)
    # Assume structural sparsity, ie if B[i,j] is not in the structural non-zeros
    # of B, then there is no given E-value between i and j
    # But if B[i,j] is, then take the -log(E-value) and clamp the result
    @. B.nzval = -log(B.nzval) 
    clamp!(B.nzval, 0, -log(clampval))
    B.nzval ./= maximum(B.nzval)
    NodeSimMeasure(B)
end
