export ghost

# GHOST's node similarity
function lapevals(A::DenseMatrix)
    D = Diagonal(vec(sum(A,1)))
    S = Diagonal(sqrt.(1./D.diag))
    L = S * (D - A) * S
    LAPACK.syevr!('N', 'V', 'U', L, -0.01, 2.01, 0, 0, 1e-3)[1]
end

# Calculates topological signature of a node
# based on its local neiborhood's spectral signature
# It induces a local neighbood graph and calculates
# its spectra. It does so for local neighborhoods graphs created
# by inducing on nodes near v from hops r=1..k
function multiscale_spectra(G::SparseMatrixCSC,k::Integer,v,debug=false)
    evals = Vector{Vector{Float64}}(k)
    dists = bfs(G,v)[1]
    pf = find(x -> (x <= k) & (x >= 0), dists)
    G = G[pf,pf]
    dists = dists[pf]
    for r = 1:k
        p = find(x -> x <= r, dists) # calculate neighbors
        if length(p) > 1
            evals[r] = lapevals(full(G[p,p]))
        else
            evals[r] = Float64[]
        end
    end
    evals
end

function ghost_signature(G::SparseMatrixCSC,k::Integer,u,n=50,mode=:bins)
    evs = multiscale_spectra(G,k,u)
    bins = -1/n:(1/n):(2+1/n)
    if mode == :struct
        s = 0.01
        p = pdf.(Normal(0,s),linspace(-7.5s,7.5s,round(7.5s*n)))
    end
    map(evs) do x
        hx = Float64.(fit(Histogram, x, bins, closed=:right).weights .> 0)

        if mode == :struct
            resize!(hx, length(hx)+length(p)-1)
            hx[length(bins):end] = 0.0
            filt!(hx, p, 1, hx) #conv(p,cx)
        end

        nx = norm(hx,1)
        nx==0 ? hx : hx./nx
    end
end

function ghost_signature(G::SparseMatrixCSC,k::Integer)
    n = size(G,1)
    evs = Vector{Vector{Vector{Float64}}}(n)
    for v = 1:n
        print("\rG: $v/$n")
        evs[u] = ghost_signature(G,k,v)
    end
    println()
    evs
end
