import Base.Order: Ordering, Forward, ord, lt

# merge sorted vectors vl and vr into v
# from indices lo to hi in v
function mergesorted!(v::AbstractVector,
                      lo::Int, hi::Int,
                      vl::AbstractVector,
                      lol::Int, hil::Int,
                      vr::AbstractVector,
                      lor::Int, hir::Int,                      
                      order::Ordering)
    c = lol
    p = lor
    nl = hil
    nr = hir
    i = lo
    @inbounds while c <= nl && p <= nr && i <= hi
        if lt(order, vr[p], vl[c])
            v[i] = vr[p]
            p = p+1
            i = i+1
        else
            v[i] = vl[c]
            c = c+1
            i = i+1
        end
    end
    @inbounds while p <= nr && i <= hi
        v[i] = vr[p]
        i = i+1
        p = p+1
    end
    @inbounds while c <= nl && i <= hi
        v[i] = vl[c]
        i = i+1
        c = c+1
    end
    v
end

function mergesorted!(v::AbstractVector, vl::AbstractVector,
                      vr::AbstractVector, order::Ordering)
    inds = indices(v,1)
    indsl = indices(vl,1)
    indsr = indices(vr,1)
    mergesorted!(v,first(inds),last(inds),vl,first(indsl),last(indsl),
                 vr,first(indsr),last(indsr),order)
end

function mergesorted!(v::AbstractVector,
                      vl::AbstractVector,
                      vr::AbstractVector,
                      lt=isless,
                      by=identity,
                      rev::Bool=false,
                      order::Ordering=Forward)
    ordr = ord(lt,by,rev,order)
    mergesorted!(v, vl, vr, ordr)
end

function mergesorted(vl::AbstractVector,
                      vr::AbstractVector,
                      lt=isless,
                      by=identity,
                      rev::Bool=false,
                      order::Ordering=Forward)
    v = similar(promote_type(typeof(vl),typeof(vr)), length(vl)+length(vr))
    ordr = ord(lt,by,rev,order)
    mergesorted!(v, vl, vr, ordr)
end
