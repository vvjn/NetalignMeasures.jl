# NetalignMeasures documentation

Network alignment measures including S3, DS3, WEC, DWEC, as well as
node similarities for network alignment including GDV similarity,
degree similarity, etc.

## Installation

NetalignMeasures can be installed as follows.

```julia
Pkg.clone("https://github.com/vvjn/NetalignMeasures.jl")
```

## Overview

```@meta
CurrentModule = NetalignMeasures
```

``` @docs
NetalignMeasure
```

## Measures

``` @docs
S3Measure
DS3Measure
WECMeasure
DWECMeasure
LCCSMeasure
CIQMeasure
GOCMeasure
GhostMeasure
ConvexCombMeasure
NodeSimMeasure
NodeSimMeasure(::Val{:gdvs},gdv1::AbstractMatrix{T},gdv2::AbstractMatrix{T}) where {T}
NodeSimMeasure(::Val{:lgraalgdvs}, gdv1::AbstractMatrix{T},gdv2::AbstractMatrix{T}) where {T}
NodeSimMeasure(::Val{:hgraalgdvs}, G1::SparseMatrixCSC,G2::SparseMatrixCSC,
                        gdv1::AbstractMatrix{T},gdv2::AbstractMatrix{T},alpha::T=0.5) where {T}
NodeSimMeasure(::Val{:pcagdvs}, gdv1::AbstractMatrix,gdv2::AbstractMatrix)
NodeSimMeasure(::Val{:relative}, v1::AbstractVector,v2::AbstractVector)
NodeSimMeasure(::Val{:degree}, G1::SparseMatrixCSC,G2::SparseMatrixCSC)
NodeSimMeasure(::Val{:clustering}, G1::SparseMatrixCSC,G2::SparseMatrixCSC)
NodeSimMeasure(::Val{:eccentricity}, G1::SparseMatrixCSC,
                        G2::SparseMatrixCSC)
NodeSimMeasure(::Val{:migraal}, Xs::AbstractMatrix...)
NodeSimMeasure(::Val{:evalues}, E::AbstractMatrix, clamp=1e-101)
NodeSimMeasure(::Val{:ghost}, G1::SparseMatrixCSC,
                        G2::SparseMatrixCSC,k::Integer=4;verbose=true)


```
