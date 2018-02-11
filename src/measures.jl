export
ConvexCombMeasure, ConvexCombScore,
S3Measure, S3Score, DS3Measure, DS3Score,
WECMeasure, WECScore, DWECMeasure, DWECScore,
NetalMeasure, NetalScore, NullMeasure, NullScore,
LCCSMeasure, LCCSScore, CIQMeasure, CIQScore,
GOCMeasure, GOCScore,
GhostMeasure, GhostScore

"""
Null measure
Always returns 0.0 score
"""
immutable NullMeasure <: NetalignMeasure
    G
    function NullMeasure(G::SparseMatrixCSC...)
        if length(G) < 2 error("Need 2 of more networks") end
        for i = 2:length(G)
            if size(G[i-1],1) > size(G[i],1) error("Network size order") end
        end
        new(G)
    end
end

immutable NullScore <: NetalignScore
    score :: Float64
end

measure(m::NullMeasure,fs::Tuple{Vector{<:Integer}}) = NullScore(0.0)
measure(m::NullMeasure,fs::Vector{<:Integer}...) = NullScore(0.0)
dim(m::NullMeasure,d::Int) = size(m.G[end],1)
dim(m::NullMeasure) = size(m.G[end],length(m.G))

"""
    ConvexCombMeasure(S::A,T::B,alpha::Float64) where {A<:NetalignMeasure,B<:NetalignMeasure}

Convex combination of two `NetalignMeasure`s. `0 <= alpha <= 1`.
"""
immutable ConvexCombMeasure{A<:NetalignMeasure,B<:NetalignMeasure} <: NetalignMeasure
    S :: A
    T :: B
    alpha :: Float64 # alpha S + (1-alpha) T
    function ConvexCombMeasure{A,B}(S::A, T::B,alpha::Float64) where
        {A<:NetalignMeasure,B<:NetalignMeasure}
        if dim(S)!=dim(T) error("Sizes of measures A and B do not match: $(dim(S)) vs $(dim(T))") end
        new(S,T,alpha)
    end
end
ConvexCombMeasure(S::A,T::B,alpha::Float64) where {A<:NetalignMeasure,B<:NetalignMeasure} =
    ConvexCombMeasure{A,B}(S,T,alpha)

immutable ConvexCombScore{A<:NetalignScore,B<:NetalignScore} <: NetalignScore
    s :: A
    t :: B
    score :: Float64
end
ConvexCombScore(s::A,t::B) where {A<:NetalignScore,B<:NetalignScore} =
    ConvexCombScore{A,B}(s,t)

function measure(m::ConvexCombMeasure,f::Vector{Int})
    a = measure(m.S,f)
    b = measure(m.T,f)
    s = m.alpha * score(a) + (1.0-m.alpha) * score(b)
    ConvexCombScore(a,b,s)
end

dim(m::ConvexCombMeasure,d::Int) = dim(m.S,d)
dim(m::ConvexCombMeasure) = dim(m,2)

"""
    S3Measure(G1::SparseMatrixCSC,G2::SparseMatrixCSC)

S3 edge conservation measure.
"""
immutable S3Measure <: NetalignMeasure
    G1 :: SparseMatrixCSC{Int,Int}
    G2 :: SparseMatrixCSC{Int,Int}
    function S3Measure(G1::SparseMatrixCSC,G2::SparseMatrixCSC)
        if size(G1,1) > size(G2,1) error("Network size order") end
        new(G1,G2)
    end
end

immutable S3Score <: NetalignScore
    Nc :: Int
    Nn :: Int
    score :: Float64 # Nc/(Nc+Nn)
end

function measure(meas::S3Measure,f::Vector{Int})
    # Assumes that all edge weights are 1
    h = view(f,1:size(meas.G1,1))
    w = nonzeros(meas.G1 + meas.G2[h,h])
    Nc = count(x->x==2,w)
    Nn = length(w) - Nc
    Nc,cr = divrem(Nc,2)
    Nn,nr = divrem(Nn,2) # need this & above so it works well with sim. ann. code
    if cr!=0 || nr!=0; error("G1 and G2 need to be symmetric"); end
    score = Nc/(Nc + Nn)
    S3Score(Nc,Nn,score)
end

dim(m::S3Measure,d::Int) = d==1 ? size(m.G1,1) : d==2 ? size(m.G2,1) : error("d")
dim(m::S3Measure) = dim(m,2)

"""
    DS3Measure(G1::SparseMatrixCSC,G2::SparseMatrixCSC)
DS3 event conservation measure.
"""
immutable DS3Measure <: NetalignMeasure
    G1 :: SparseMatrixCSC{Events,Int}
    G2 :: SparseMatrixCSC{Events,Int}
    function DS3Measure(G1::SparseMatrixCSC,G2::SparseMatrixCSC)
        if size(G1,1) > size(G2,1) error("Bad args") end
        new(G1,G2)
    end
end

immutable DS3Score <: NetalignScore
    Tc :: Float64
    Tn :: Float64
    score :: Float64 # Tc/(Tc+Tn)
end

function measure(m::DS3Measure,f::Vector{Int})
    G1 = m.G1
    G2 = m.G2
    h = view(f,1:size(G1,1))
    # Important: assumes edges (of G1 and G2) are symmetric,
    I1,J1,V1 = findnz(G1)
    I2,J2,V2 = findnz(G2[h,h])
    S = sparse(vcat(I1,I2),vcat(J1,J2),vcat(V1,V2),
               size(G1,1),size(G1,1),
               mergesorted)
    Tc,Tn = 0.0,0.0
    for tis in nonzeros(S)
        Tcp,Tnp = cet_ncet(tis)
        Tc += Tcp
        Tn += Tnp
    end
    Tc /= 2
    Tn /= 2 # this & above required for sim. ann. code
    score = Tc / (Tc + Tn)
    DS3Score(Tc,Tn,score)
end

dim(m::DS3Measure,d::Int) =
    d==1 ? size(m.G1,1) : d==2 ? size(m.G2,1) : error("d")
dim(m::DS3Measure) = dim(m,2)

"""
    WECMeasure(G1::SparseMatrixCSC,G2::SparseMatrixCSC,S::AbstractMatrix)

Weighted edge conservation from WAVE.
"""
immutable WECMeasure{T<:AbstractMatrix} <: NetalignMeasure
    G1 :: SparseMatrixCSC{Int,Int}
    G2 :: SparseMatrixCSC{Int,Int}
    S :: T # matrix type
    function WECMeasure{T}(G1::SparseMatrixCSC,
                           G2::SparseMatrixCSC,S::T) where {T<:AbstractMatrix}
        if size(G1,1) > size(G2,1) || size(S,1)!=size(G1,1) ||
             size(S,2)!=size(G2,1) error("Network/matrix dims") end
        new(G1,G2,S)
    end
end
WECMeasure(G1::SparseMatrixCSC,G2::SparseMatrixCSC,S::AbstractMatrix) =
    WECMeasure{typeof(S)}(G1,G2,S)

immutable WECScore <: NetalignScore
    score :: Float64
end

function measure(m::WECMeasure,f::Vector{Int})
    # Assumes that all edge weights are 1
    h = view(f,1:size(m.G1,1))
    I,J,w = findnz(m.G1 + m.G2[h,h])
    wix = w.==2  # conserved edges
    I = I[wix]
    KI = sub2ind(size(m.S), I, h[I])
    score = sum(m.S[KI]) / min(nnz(m.G1),nnz(m.G2))
    WECScore(score)
end

dim(m::WECMeasure,d::Int) = size(m.S,d)
dim(m::WECMeasure) = dim(m,2)

"""
    DWECMeasure(G1::SparseMatrixCSC{Events},G2::SparseMatrixCSC{Events},S::AbstractMatrix)

Dynamic weighted edge conservation from DynaWAVE.
"""
immutable DWECMeasure{T} <: NetalignMeasure
    # IMP: symmetric "adjacency" matrices
    G1 :: SparseMatrixCSC{Events,Int}
    G2 :: SparseMatrixCSC{Events,Int}
    S :: T # matrix type
    activitysum1 :: Float64
    activitysum2 :: Float64
    function DWECMeasure{T}(G1::SparseMatrixCSC{Events},
                            G2::SparseMatrixCSC{Events},S::T) where {T}
        if size(G1,1) > size(G2,1) || size(S,1)!=size(G1,1) ||
             size(S,2)!=size(G2,1) error("Bad args") end
        new(G1, G2, S, networkactivity(G1), networkactivity(G2))
    end
end
DWECMeasure(G1::SparseMatrixCSC{Events}, G2::SparseMatrixCSC{Events},S::AbstractMatrix) =
    DWECMeasure{typeof(S)}(G1,G2,S)

immutable DWECScore <: NetalignScore
    score :: Float64
end

function measure(m::DWECMeasure,f::Vector{Int})
    G1 = m.G1
    G2 = m.G2
    h = view(f,1:size(G1,1))
    I1,J1,V1 = findnz(G1)
    I2,J2,V2 = findnz(G2[h,h])
    I,J,W = findnz(sparse(vcat(I1,I2),vcat(J1,J2),vcat(V1,V2),
                          size(G1,1),size(G1,1), mergesorted))
    w = map(es -> cet_ncet(es)[1], W)
    KI = sub2ind(size(m.S), I, h[I])
    score = dot(w,m.S[KI]) / min(m.activitysum1,m.activitysum2)
    DWECScore(score)
end

dim(m::DWECMeasure,d::Int) = size(m.S,d)
dim(m::DWECMeasure) = size(m.S,2)

"""
    NetalMeasure(G1::SparseMatrixCSC,G2::SparseMatrixCSC)

Measure from NETAL paper.
"""
immutable NetalMeasure <: NetalignMeasure
    G1 :: SparseMatrixCSC{Int,Int}
    G2 :: SparseMatrixCSC{Int,Int}
    dep1 :: Vector{Float64}
    dep2 :: Vector{Float64}
    expci1 :: Vector{Float64}
    expci2 :: Vector{Float64}
    maxdeg1 :: Int
    maxdeg2 :: Int
    function NetalMeasure(G1::SparseMatrixCSC,G2::SparseMatrixCSC)
        if size(G1,1) > size(G2,1) error("Bad args") end
        deg1 = degree(G1)
        deg2 = degree(G2)
        dep1 = dependency(G1,deg1)
        dep2 = dependency(G2,deg2)
        expci1 = expci(G1,dep1)
        expci2 = expci(G2,dep2)
        new(G1,G2, dep1,dep2, expci1,expci2, maximum(deg1),maximum(deg2))
    end
end

immutable NetalScore <: NetalignScore
    score :: Float64
end

function measure(meas::NetalMeasure,f::AbstractVector{Int})
    # Assumes that all edge weights are 1
    h = view(f,1:size(meas.G1,1))
    I,J,w = findnz(meas.G1 + meas.G2[h,h])
    score = count(x -> x==2, w) / min(nnz(meas.G1),nnz(meas.G2))
    NetalScore(score)
end

dim(m::NetalMeasure,d::Int) =
    d==1 ? size(m.G1,1) : d==2 ? size(m.G2,1) : error("d")
dim(m::NetalMeasure) = dim(m,2)

"""
   LCCSMeasure(G::SparseMatrixCSC...)

Largest common conserved sub-graph, from multiMAGNA++ paper. Works for >2 networks.
"""
immutable LCCSMeasure <: NetalignMeasure
    G
    function LCCSMeasure(G::SparseMatrixCSC...)
        if length(G) < 2 error("Need 2 of more networks") end
        new(G)
    end
end

immutable LCCSScore <: NetalignScore
    edgefrac :: Float64
    nodefrac :: Float64
    edgecov :: Int
    nodecov :: Int
    score :: Float64
end

function compositegraph(G::Tuple,
                        fs::Tuple{Vector{<:Integer}})
    k = length(G)
    ncs = size(G[end],1)
    I,J,_ = findnz(G[end])
    M = Vector(k)
    M[k] = sparse(I,J,1,ncs,ncs)
    for l = 1:k-1
        f = fs[l]
        Ic,Jc,Vc = findnz(G[l])
        fi = f[Ic]
        fj = f[Jc]
        M[l] = sparse(fi,fj,1,ncs,ncs)
        append!(I, fi)
        append!(J, fj)
    end
    M,sparse(I,J,1,ncs,ncs)
end

# Only works for 1-1 MNA
# fs[i] is permutation from G[i] to G[end]
function measure(meas::LCCSMeasure,fs::Tuple{Vector{<:Integer}})
    for i = 2:length(meas.G)
        if size(meas.G[i-1],1) > size(meas.G[i],1) error("Network size order") end
    end
    k = length(meas.G)
    length(fs) == k-1 || error("number of permutations")
    ncs = size(meas.G[end],1)
    M,R = compositegraph(meas.G, fs) # wrt G[end]

    I,J,V = findnz(R)
    ix = V.==k
    I = I[ix]
    J = J[ix]
    A = sparse(I,J,1,size(R,1),size(R,2))

    if nnz(A)==0
        efrac,nfrac,nc,ec = (0.0,0.0,0.0,0,0)
    end

    Acc,p = largest_component(A)
    nc = sum(p)
    ec = nnz(Acc)

    # find the maximum possible number of edges for each graph
    # want the inverse map from cluster to each graph
    # see magna 2014 bioinf paper
    emin = minimum([nnz(M[l][p,p]) for l in 1:k])
    nmin = minimum(map(G->size(G,1), meas.G))
    if nmin!=0 && emin!=0
        efrac = ec/emin
        nfrac = nc/nmin
    else
        nfrac = 0.0
        efrac = 0.0
    end
    return LCCSScore(efrac,nfrac,ec,nc,sqrt(efrac*nfrac))
end

immutable CIQMeasure <: NetalignMeasure
    G
    function CIQMeasure(G::SparseMatrixCSC...)
        if length(G) < 2 error("Need 2 of more networks") end
        new(G)
    end
end
immutable CIQScore <: NetalignScore
    ciqn :: Float64
    ciqd :: Float64
    score :: Float64
end
# Only works for 1-1 MNA
# fs[i] is permutation from G[i] to G[end]
function measure(meas::CIQMeasure,fs::Tuple{Vector{<:Integer}})
    for i = 2:length(meas.G)
        if size(meas.G[i-1],1) > size(meas.G[i],1) error("Network size order") end
    end
    
    k = length(meas.G)
    length(fs) == k-1 || error("number of permutations")
    ncs = size(meas.G[end],1)
    M,R = compositegraph(meas.G, fs)

    if nnz(R)==0
        ciqn,ciqd = (0.0,0.0)
    end

    scl = ones(Int,ncs,k)
    for l = 1:k-1
        f = fs[l]
        for i = 1:size(meas.G[l],1)
            scl[f[i],l] += 1
        end
    end

    u,v,r = findnz(triu(R,1))
    scl = scl .> 0
    # >0 if there is a node in each cluster that both come from the same network
    s = vec(sum(scl[u,:] .& scl[v,:],2))

    cs = r./s
    cs[(s.==0) .| (r.<=1)] = 0.0

    ciqd = sum(r)
    ciqn = dot(r,cs)

    ciq = ciqd!=0 ? ciqn/ciqd : 0.0
    CIQScore(ciqn,ciqd,ciq)
end

function gofilter(known)
    counts = Dict{String,Int}()
    for node in known, goterm in node
        counts[goterm] = get(counts,goterm,0) + 1
    end
    goodterms = keys(filter((k,v) -> v>1, counts))
    map(goterms -> intersect(goterms,goodterms), known)
end

immutable GOCMeasure <: NetalignMeasure
    known1 :: Vector # known[i] contains annotated go for node i
    known2 :: Vector
    k :: Int
    function GOCMeasure(known1::Vector,known2::Vector,k::Int=1)
        if length(known1) > length(known2) error("Network size order") end
        new(gofilter(known1),gofilter(known2),k)
    end
end
immutable GOCScore <: NetalignScore
    num :: Int
    den :: Int
    score :: Float64
end
function measure(meas::GOCMeasure,f::Vector{<:Integer})
    s = 0
    c = 0
    for (x,y) in zip(meas.known1,meas.known2[view(f,1:length(meas.known1))])
        ixy = length(intersect(x,y))
        mxy = min(length(x),length(y))
        if mxy >= meas.k
            s += Int(ixy > 0) #/mxy
            c += 1
        end
    end
    GOCScore(s, c, c==0 ? 0.0 : s/c)
end

include("ghost.jl")

"""
GHOST signature and node similarity
"""
immutable GhostMeasure <: NetalignMeasure
    G1 :: SparseMatrixCSC{Int,Int}
    G2 :: SparseMatrixCSC{Int,Int}
    evs1::Vector{Vector{Vector{Float64}}}
    evs2::Vector{Vector{Vector{Float64}}}
    k :: Int
    function GhostMeasure(G1::SparseMatrixCSC,G2::SparseMatrixCSC,k=4;verbose=true)
        if size(G1,1) > size(G2,1) error("Network size order") end
        verbose && println("Calculating GHOST signatures")
        m = size(G1,1)
        n = size(G2,1)
        evs1 = ghost_signature(G1,k)
        evs2 = ghost_signature(G2,k)
        new(G1,G2,evs1,evs2,k)
    end
end

immutable GhostScore <: NetalignScore
    score :: Float64
end

function measure(meas::GhostMeasure,f::Vector{Int})
    s = 0.0
    n = length(f)
    for i = 1:n
        s += score(meas, i, f[i])
    end
    GhostScore(s/n)
end

function score(meas::GhostMeasure, u::Int, v::Int)
    ev1 = meas.evs1[u]
    ev2 = meas.evs2[v]
    d = 0.0
    for r = 1:meas.k
        d += sqrt(js_divergence(ev1[r], ev2[r]))
    end
    1.0 - d/meas.k
end

dim(m::GhostMeasure,d::Int) = d==1 ? size(m.G1,1) : d==2 ? size(m.G2,1) : error("d")
dim(m::GhostMeasure) = dim(m,2)
