export NodeSimMeasure, NodeSimScore

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

type PCA
    mean::AbstractArray
    proj::AbstractMatrix       # projection
    prinvars::AbstractVector # variance of projX*projX' instead of cov(projX) like fit(PCA)
end
"""
Input:
X : m x n matrix
Each column of X is a feature vector
Get the first dimensions containing k fraction of variance

Output:
    PCA object
    with mean, projection matrix, and variances
    U is projection matrix
    projection is T = U' * X
    U*T =~ X, where X is mean centered
"""
function pca(X::AbstractMatrix,k::Real;mu=mean(X,2))
    m,n = size(X)
    if mu!=0
        X = X .- vec(mu)
    else
        mu = zeros(0,0)
    end
    if m <= n # feature size is less than number of samples
        C = X*X'
        res = eigfact!(C)
        d = res.values
        U = res.vectors
        ix = sortperm(d,rev=true)
        d = d[ix]
        clamp!(d, 0.0, Inf)
    else
        C = X'*X
        res = eigfact!(C)
        d = res.values
        V = res.vectors
        ix = sortperm(d,rev=true)
        d = d[ix]
        clamp!(d, 0.0, Inf)
    end
    if !(0.0 <= k <= 1.0) error("k needs to be between 0 and 1") end
    k = findfirst(cumsum(d)./sum(d) .>= k)
    d = d[1:k]
    ix = ix[1:k]
    if m <= n
        U = U[:,ix]
    else
        V = V[:,ix]
        U = X*V # convert evecs from X'*X to X*X'. the evals are the same.
        s = sqrt.(d)
        U = U ./ s'
    end
    PCA(mu, U, d)
end

"The GDV similarity calculated using PCA from Yuriy's dynamic graphlets paper"
function NodeSimMeasure(::Val{:pcagdvs}, gdv1::AbstractMatrix,gdv2::AbstractMatrix)
    X = hcat(gdv1',gdv2') # switch to make columns feature vectors
    X = X[vec(var(X,2)) .>= 1e-5,:] # remove features w/ no variance
    X .-= mean(X,2)
    X ./= std(X,2) # standardize
    res = pca(X,0.99,mu=0) #fit(PCA,X; pratio=0.99, mean=0)
    Y = sqrt.(res.prinvars + 1e-5) .\ res.proj' * X # project and whiten
    R = 1 .- pairwise(CosineDist(), Y[:,1:size(gdv1,1)], Y[:,size(gdv1,1)+1:end]) ./ 2
    NodeSimMeasure(R)
end

"The GDV similarity calculated using CCA"
function NodeSimMeasure(::Val{:ccagdvs}, gdv1::AbstractMatrix,gdv2::AbstractMatrix)
    X = gdv1'
    Y = gdv2'
    res = fit(CCA,X,Y)
    R = 1 .- pairwise(CosineDist(), xtransform(res,X), ytransform(res,Y)) ./ 2
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

"""Transforms blast e-values to nodes similarities.
# Arguments
- `clamp = 1e-101` : Clamp E-values to a minimum of `clamp`."""
function NodeSimMeasure(::Val{:evalues}, E::AbstractMatrix, clamp=1e-101)
    # -log(e-value) transform
    B = -log.(E)
    # Clamp E-value to a minimum of `clamp`
    # That is, clamp -log(E-value) to max of -log(clampval)
    # Otherwise we get issues with -log(E-value) being Inf or some
    # ridiculously huge number
    clamp!(B, 0.0, -log(clamp))
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
    clamp!(B.nzval, 0.0, -log(clamp))
    B.nzval ./= maximum(B.nzval)
    NodeSimMeasure(B)
end

"GHOST node similiarity"
function NodeSimMeasure(::Val{:ghost}, G1::SparseMatrixCSC,
                        G2::SparseMatrixCSC,k::Integer=4;verbose=true)
    m = size(G1,1)
    n = size(G2,1)
    gmeas = GhostMeasure(G1,G2,k,verbose=verbose)
    verbose && println("Calculating GHOST node similarities")
    R = zeros(Float64,m,n)
    for v = 1:n, u = 1:m
        R[u,v] = score(meas,u,v)
    end
    NodeSimMeasure(R)
end

# interaction matrix = min(G1.expci[i],G2.expci[j])/max(G1.maxdeg,G2.maxdeg)
function NodeSimMeasure(::Val{:interactionscore}, G1::SparseMatrixCSC,G2::SparseMatrixCSC)
    NodeSimMeasure(expci(G1,G2))
end
