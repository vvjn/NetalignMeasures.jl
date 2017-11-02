export alnfillrandom, alnfillrandom!, aln2perm
# Given an alignment f from set V1 to V2, that does not
# map all points in V1, map the rest of the points randomly
alnfillrandom(f::AbstractVector{<:Integer}, n::Integer) = alnfillrandom!(copy(f),n)
function alnfillrandom!(f::AbstractVector{<:Integer}, n::Integer)
    zpts = find(iszero, f)
    nzpts = find(!iszero, f)
    remf = setdiff(1:n, f[nzpts])
    shuffle!(remf)
    f[zpts] = remf[1:length(zpts)]
    f
end

function aln2perm(f::AbstractVector{Int},n::Integer)
    m = length(f)
    p = zeros(Int, n)
    p[1:m] = f
    mapped_to = zeros(Int,n)
    for i = 1:n
        if p[i] > 0
            mapped_to[p[i]] = 1
        end
    end
    j = 1
    for i = 1:n
        if mapped_to[i] == 0
            mapped_to[j] = i
            j += 1
        end
    end
    rp = randperm(n-m)
    for i = 1:(n-m)
        p[m+i] = mapped_to[rp[i]]
    end
    p
end


