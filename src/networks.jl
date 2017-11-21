export adjnodes, adjnzval, degree, dependency, expci

""" degree vector of an adjacency graph """
degree(G::SparseMatrixCSC) = vec(full(sum(G,1)))

"""
 Given an adj. matrix of a graph, return index of nodes
 adjacent to node i
"""
adjnodes(A::SparseMatrixCSC,i) = A.rowval[nzrange(A,i)]

"""
 Given an adj. matrix of a graph, return values of nodes
 adjacent to node i
"""
adjnzval(A::SparseMatrixCSC,i) = A.nzval[nzrange(A,i)]

"""
probability that any of i's neighbors will be conserved (netal)
"""
function dependency(G::SparseMatrixCSC, deg=degree(G))
    dep = 1./deg
    map!(x -> ifelse(isinf(x), 0, x), dep, dep)
end

"""
 approx expected number of conserved interactions (netal)
"""
function expci(G::SparseMatrixCSC, dep=dependency(G))
    map(i -> sum(dep[adjnodes(G,i)]), 1:size(G,1))
end

"""
expci(..)[i,j] = I[i,j] = expected # of conserved edges if i is aligned to j
"""
function expci(G1::SparseMatrixCSC,G2::SparseMatrixCSC,
               deg1=degree(G1),
               deg2=degree(G2),
               expci1=expci(G1,dependency(G1,deg1)),
               expci2=expci(G2,dependency(G2,deg2));
               normalize=true)
    z = normalize ? max(maximum(deg1),maximum(deg2)) : 1
    min.(expci1./z, (expci2./z)')
end
