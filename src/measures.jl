export
ConvexCombMeasure, ConvexCombScore,
S3Measure, S3Score, DS3Measure, DS3Score,
WECMeasure, WECScore, DWECMeasure, DWECScore,
NetalMeasure, NetalScore, NullMeasure, NullScore

"""
Null measure
Always returns 0.0 score
"""
immutable NullMeasure <: NetalignMeasure
    G
    function NullMeasure(G::SparseMatrixCSC...)
        for i = 2:length(G)
            if size(G[i-1],1) > size(G[i],1) error("Network size order") end
        end
        new(G)
    end
end

immutable NullScore <: NetalignScore
    score :: Float64
end

measure(m::NullMeasure,f::Vector{<:Integer}) = NullScore(0.0)
dim(m::NullMeasure,d::Int) = size(m.G[d],1)
dim(m::NullMeasure) = size(m.G[d],length(m.G))

"""
Convex combination of two NAMs
"""
immutable ConvexCombMeasure{A<:NetalignMeasure,B<:NetalignMeasure} <: NetalignMeasure
    S :: A
    T :: B
    alpha :: Float64 # alpha S + (1-alpha) T
    function ConvexCombMeasure{A,B}(S::A, T::B,alpha::Float64) where
        {A<:NetalignMeasure,B<:NetalignMeasure}
        if dim(S)!=dim(T) error("Bad args") end
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
S3 edge conservation measure
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
DS3 event conservation
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
Weighted edge conservation from WAVE
"""
immutable WECMeasure{T} <: NetalignMeasure
    G1 :: SparseMatrixCSC{Int,Int}
    G2 :: SparseMatrixCSC{Int,Int}
    S :: T # matrix type
    function WECMeasure{T}(G1::SparseMatrixCSC,
                           G2::SparseMatrixCSC,S::T) where {T}
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
Weighted dynamic edge conservation from DynaWAVE
"""
immutable DWECMeasure{T} <: NetalignMeasure
    # IMP: symmetric "adjacency" matrices
    G1 :: SparseMatrixCSC{Events,Int}
    G2 :: SparseMatrixCSC{Events,Int}
    S :: T # matrix type
    activitysum1 :: Float64
    activitysum2 :: Float64
    function DWECMeasure{T}(G1::SparseMatrixCSC,
                            G2::SparseMatrixCSC,S::T) where {T}
        if size(G1,1) > size(G2,1) || size(S,1)!=size(G1,1) ||
             size(S,2)!=size(G2,1) error("Bad args") end
        new(G1, G2, S, totalactivity(G1), totalactivity(G2))
    end
end

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

"""
Measure from NETAL paper
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

function measure(meas::NetalMeasure,f::Vector{Int})
    # Assumes that all edge weights are 1
    h = view(f,1:size(meas.G1,1))
    I,J,w = findnz(meas.G1 + meas.G2[h,h])
    score = count(x -> x==2, w) / min(nnz(meas.G1),nnz(meas.G2))
    NetalScore(score)
end

dim(m::NetalMeasure,d::Int) =
    d==1 ? size(m.G1,1) : d==2 ? size(m.G2,1) : error("d")
dim(m::NetalMeasure) = dim(m,2)
