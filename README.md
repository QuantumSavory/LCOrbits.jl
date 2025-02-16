A Julia wrapper for the resources at

- https://www.ii.uib.no/~larsed/entanglement/ (related to https://arxiv.org/abs/1011.5464)
- https://zenodo.org/records/3757948 (related to https://arxiv.org/abs/1910.03969)


All of the wrapped resources are under the original authors's copyright.
This repository is just for convenient hosting and wrapping in a Julia library.

## `lcorbits_summary`

Use to get the database from `arxiv:1011.5464`

Keys are:

- `EC` - the equivalency class index (in the standard order suggested by the authors)
- `orbitsize` - the size of the orbit
- `vertices` - the number of vertices of the graph
- `contains2color` - whether there is a two colorable graph in the orbit
- `minedges` - the minimal edges (and the corresponding chromatic index, and how many non-isomorphic graphs have that property)
- `minchromidx` - the minimal chromatic index (and the corresponding number of edges, and how many non-isomorphic graphs have that property)
- `same_min_edgeschromidx` - whether the previous two are the same
- `schmidtmin` - Schmidt measure lower bound
- `schmidtmax` - Schmidt measure upper bound
- `schmidtsame` - whether the bounds match
- `minedge_repr` - the minimal edge representation
- `minedge_repr_colors` - coloring for the minimal edge representation
- `minchromidx_repr` - the minimal chromatic index representation
- `minchromidx_repr_colors` - coloring for the minimal chromatic index representation
- `same_min_edgeschromidx_repr` - whether the minimal edge representation and the minimal chromatic index representation are the same
- `bipartrankidx6` - rank index for bipartite splits with 6 vertices in the smaller partition
- `bipartrankidx5` - rank index for bipartite splits with 5 vertices in the smaller partition
- `bipartrankidx4` - rank index for bipartite splits with 4 vertices in the smaller partition
- `bipartrankidx3` - rank index for bipartite splits with 3 vertices in the smaller partition
- `bipartrankidx2` - rank index for bipartite splits with 2 vertices in the smaller partition
- `cardmult` - the cardinality-multiplicities


```julia-repr
julia> using LCOrbits, GraphMakie, CairoMakie

julia> orbits = LCOrbits.lcorbits_summary(10); # only graphs of size 10

julia> orbit42_10 = orbits[42,:] # the 42nd orbit
DataFrameRow
 Row │ EC     orbitsize  vertices  contains2color  minedges     ⋯
     │ Int64  Int64      Int64     Bool            NamedTuple   ⋯
─────┼────────────────────────────────────────────────────────────────────────────────
  42 │   628        120        10            true  (minedges = 9, chromidx = 5, size ⋯

julia> graphplot(orbit42_10.minedge_repr; edge_color=orbit42_10.minedge_repr_colors)
```

![](./docs/src/graphplot.png)