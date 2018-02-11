export Events, networkactivity, nodeactivity, cet_ncet, fixevents,
fixevents!, widenevents!, snapshot, snapshots, mintime, maxtime,
meannodeactivity, numevents, mergeevents, mergeevents!, flatten

# Functions for working with dynamic networks
import Base: zero, ==

"""
An dynamic edge in a dynamic network is defined as a set of events. An
event is an interaction between two nodes from start time t0 to stop
time t1, t0 <= t1, and is represent as a tuple (t0,t1). The set of
events is represented as a Vector, sorted wrt start time t0, and if
two start times are equal, then stop time t1.  Note: Unless specified,
it is assumed that two events in an event set will not overlap,
i.e. there exists no two tuples (t0,t1) and (s0,s1)
s.t. intersect([t0,t1],[s0,s1]) != null

A dynamic network can be represented using a
SparseMatrixCSC{Events} structure where an element in the
matrix (i.e. an edge) is represented using the following structure
"""
immutable Events
    timestamps:: Vector{Tuple{Float64,Float64}}
end
Events() = Events(Vector{Tuple{Float64,Float64}}())
zero(::Type{Events}) = Events()
(==)(x::Events, y::Events) = x.timestamps == y.timestamps
mergesorted(x::Events, y::Events) =
    Events(mergesorted(x.timestamps,y.timestamps))

numevents(G::SparseMatrixCSC{Events}) =
    sum(map(x->length(x.timestamps), G.nzval))

flatten(G::SparseMatrixCSC{<:Number,Int}) = G
flatten(G::SparseMatrixCSC{Events,Int}) =
    SparseMatrixCSC(G.m,G.n,G.colptr,G.rowval,ones(Int,length(G.nzval)))

"""
    Total time during which events in event set are active
    i.e. active duration of an edge
"""
nodeactivity(evs::Events) = sum(x -> x[2]-x[1], evs.timestamps)

meannodeactivity(G::SparseMatrixCSC{Events}) = mean(map(nodeactivity, G.nzval))
"""
Input: adj. matrix of a dynamic network
Output: adj. matrix of static network containing corresponding
 duration of each edge in the dynamic network
"""
nodeactivity(G::SparseMatrixCSC{Events}) =
    SparseMatrixCSC(G.m,G.n,copy(G.colptr),copy(G.rowval),
                    map(nodeactivity, G.nzval))

"""
Sum of active duration over all events in G
"""
networkactivity(G::SparseMatrixCSC{Events}) = sum(nodeactivity, G.nzval)

function mintime(G::SparseMatrixCSC{Events})
    t_min = Inf
    for events in G.nzval, x in events
        t_min = min(t_min,x[1])
    end
    t_min
end
function maxtime(G::SparseMatrixCSC{Events})
    t_max = -Inf
    for events in G.nzval, x in events
        t_max = max(t_max,x[2])
    end
    t_max
end

"""
    If event intersects with [t_s,t_e] then add edges to snapshot
"""
function snapshot(G::SparseMatrixCSC{Events},t_s::Real,t_e::Real)
    m,n = size(G)
    I = Int[]
    J = Int[]
    for col = 1:n
        for j = nzrange(G,col)
            row = G.rowval[j]
            events = G.nzval[j]
            Tc,Tn = cet_ncet(events,[(t_s,t_e)])
            if Tc > 0
                push!(I,row)
                push!(J,col)
            end
        end
    end
    sparse(I,J,1,m,n)
end

""" Convert to snapshots with window size t_w
"""
function snapshots(G::SparseMatrixCSC{Events},t_w::Real,
                   t_min::Real=-Inf,t_max::Real=Inf)
    if isinf(t_min)
        t_min = mintime(G)
    end
    if isinf(t_max)
        t_max = maxtime(G)
    end
    map(t -> snapshot(G,t,t+t_w), t_min:t_w:t_max)
end

"""
Make each event in network wider by adding making each event
occur `pre` time ahead and `post` time later
"""
function widenevents!(G::SparseMatrixCSC{Events},pre=0.5,post=0.5)
    m,n = size(G)
    for i = 1:n
        for j = nzrange(G,i)
            row = G.rowval[j]
            events = G.nzval[j]
            for k = 1:length(events)
                x = events[k]
                events[k] = (x[1]-pre,x[2]+post)
            end
            fixevents!(G.nzval[j])
        end
    end
    G
end

"""
Given a vector of timestamps (sorted wrt the start time),
with possibly overlapping durations,
return a vector of timestamps with no overlapping durations
"""
fixevents(x::Events) = Events(fixevents(x.timestamps))
function fixevents(tis::Vector)
    n = length(tis)
    tout = similar(tis,0)
    a,b = tis[1]
    i = 2
    while i <= n
        c,d = tis[i]
        i += 1
        if b < c
            push!(tout,(a,b))
            a,b = c,d
        else
            b = d
        end
    end
    push!(tout,(a,b))
    tout
end
mergeevents(x) = fixevents(x)

"""
Given a vector of timestamps (sorted wrt the start time),
with possibly overlapping durations,
modify vector of timestamps s.t. it has no overlapping durations
"""
function fixevents!(x::Events)
    fixevents!(x.timestamps)
    x
end
function fixevents!(tis::Vector)
    n = length(tis)
    a,b = tis[1]
    i = 2
    j = 1
    while i <= n
        c,d = tis[i]
        i += 1
        if b < c
            tis[j] = (a,b)
            j += 1
            a,b = c,d
        else
            b = d
        end
    end
    tis[j] = (a,b)
    resize!(tis,j)
end
mergeevents!(x) = fixevents!(x)

"""
Given two vectors of timestamps (sorted wrt the start time),
with possibly overlapping durations,
return a vector of timestamps with no overlapping durations
"""
fixevents(x::Events, y::Events) =
    Events(fixevents(x.timestamps, y.timestamps))
fixevents(tis1::Vector, tis2::Vector) =
    fixevents!(mergesorted(tis1,tis2))
mergeevents(x,y) = fixevents(x,y)

"""
Calculates CET and NCET of the first `n` elements in a sorted array of events
As in, the events from two node pairs are combined together
and `tis` is the combined events
Sortedness is not checked
"""
cet_ncet(x::Events, n::Int=length(x.timestamps)) =
    cet_ncet(x.timestamps, n)
function cet_ncet(tis::Array,n::Int=length(tis))
    T = eltype(eltype(tis))
    if n > length(tis)
        error("n is larger than array size")
    end
    n==0 && return (zero(T),zero(T))
    (T_c,T_n) = (zero(T),zero(T))
    a,b = tis[1]
    i = 2
    while i <= n
        c,d = tis[i]
        i += 1
        if b <= c
            T_n += b-a
            a,b = c,d
        elseif b > d
            T_n += c-a
            T_c += d-c
            a,b = d,b
        else
            T_n += c-a
            T_c += b-c
            a,b = b,d
        end
    end
    T_n += b-a
    (T_c,T_n)
end

"""
Calculates CET and NCET of the events in tsl and tsr
Assumes the events in each array are non-overlapping
Sortedness is not checked
"""
cet_ncet(x::Events,y::Events) =
    cet_ncet(x.timestamps, y.timestamps)
function cet_ncet(tsl::Array,tsr::Array)
    T = eltype(eltype(tsl))
    nl = length(tsl)
    nr = length(tsr)
    (T_c,T_n) = (zero(T),zero(T))

    il = 1
    ir = 1

    a,b = (zero(T),zero(T))

    if il<nl && ir<nr
        if tsl[il][1] <= tsr[ir][1]
            a,b = tsl[il]
            il += 1
        else
            a,b = tsr[ir]
            ir += 1
        end
    end

    while il <= nl && ir <= nr
        if tsl[il][1] <= tsr[ir][1]
            c,d = tsl[il]
            il += 1
        else
            c,d = tsr[ir]
            ir += 1
        end

        if b <= c
            T_n += b-a
            a,b = c,d
        elseif b > d
            T_n += c-a
            T_c += d-c
            a,b = d,b
        else
            T_n += c-a
            T_c += b-c
            a,b = b,d
        end
    end

    while il <= nl
        c,d = tsl[il]
        il += 1

        if b <= c
            T_n += b-a
            a,b = c,d
        elseif b > d
            T_n += c-a
            T_c += d-c
            a,b = d,b
        else
            T_n += c-a
            T_c += b-c
            a,b = b,d
        end
    end

    while ir <= nr
        c,d = tsr[ir]
        ir += 1

        if b <= c
            T_n += b-a
            a,b = c,d
        elseif b > d
            T_n += c-a
            T_c += d-c
            a,b = d,b
        else
            T_n += c-a
            T_c += b-c
            a,b = b,d
        end
    end

    T_n += b-a
    (T_c,T_n)
end
