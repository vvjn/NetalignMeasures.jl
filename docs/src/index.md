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
GhostMeasure
ConvexCombMeasure
NodeSimMeasure
NullMeasure
NetalMeasure

```

## Dynamic Networks

### Types

``` @docs
Events
```

### Functions

``` @docs
cet_ncet
networkactivity
snapshots
nodeactivity
fixevents
fixevents!
snapshot
widenevents!
```

## Static Networks

### Functions

``` @docs
degree
dependency
expci
adjnodes
adjnzval
```

# Other functions

``` @docs
pca
```
