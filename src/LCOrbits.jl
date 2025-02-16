module LCOrbits

using Artifacts: @artifact_str
using CSV: CSV
using DataFrames: DataFrame, select!, Not
using Graphs: Graph, add_edge!, ne, edges, Edge

# from larsed https://arxiv.org/abs/1011.5464
lcf02 = joinpath(artifact"lcorbit02", "entanglement2")
lcf03 = joinpath(artifact"lcorbit03", "entanglement3")
lcf04 = joinpath(artifact"lcorbit04", "entanglement4")
lcf05 = joinpath(artifact"lcorbit05", "entanglement5")
lcf06 = joinpath(artifact"lcorbit06", "entanglement6")
lcf07 = joinpath(artifact"lcorbit07", "entanglement7")
lcf08 = joinpath(artifact"lcorbit08", "entanglement8")
lcf09 = joinpath(artifact"lcorbit09", "entanglement9")
lcf10 = joinpath(artifact"lcorbit10", "entanglement10")
lcf11 = joinpath(artifact"lcorbit11", "entanglement11")
lcf12 = joinpath(artifact"lcorbit12", "entanglement12")

lcfs = [lcf02, lcf03, lcf04, lcf05, lcf06, lcf07, lcf08, lcf09, lcf10, lcf11, lcf12]

# from sammorley https://arxiv.org/abs/1910.03969

function parsetuple(str)
    tuple((parse(Int,i) for i in split(str[2:end-1],','))...)
end

function parsevec(str)
    Int[parse(Int,i) for i in split(str[2:end-1],',')]
end

function parsegraphcol(str,vertices)
    groups = split(str[2:end-1], ")(")
    ncolors = length(groups)
    graph = Graph(vertices)
    colors_dict = Dict{Edge{Int},Int}()
    for (color, group) in enumerate(groups)
        edges = split(group,',')
        for edge in edges
            e1, e2 = minmax(parse.(Int, split(edge,'-'))...)
            p = Edge(e1+1,e2+1) # switch to Julia indexing
            add_edge!(graph, p)
            colors_dict[p] = color
            #println(p)
        end
    end
    colors = Int[colors_dict[e] for e in edges(graph)]
    return graph, colors
end

bipartranknames(n) = (n>=4 ? [Symbol(:_bipartrankidx,i) for i in n÷2:-1:2] : [:bipartrank0])

function get_raw_df_larsed(n)
    CSV.read(
        LCOrbits.lcfs[n-1],
        DataFrame;
        header=[:EC, :orbitsize, :vertices, :_minedges, :_minchromidx, :_schmidt, bipartranknames(n)..., :_cardmult, :contains2color, :_minedgerepr, :_minchromidxrepr],
        delim='\t',
        truestrings=["yes"],
        falsestrings=["no"],
        missingstring=nothing,
        ignorerepeated=true,
        types=Dict(:_schmidt=>CSV.String15)
    )
end

function clean_up_df_larsed(n)
    df = get_raw_df_larsed(n)
    l = size(df, 1)

    minedges = Vector{@NamedTuple{minedges::Int, chromidx::Int, size::Int}}(undef, l)
    minchromidx = Vector{@NamedTuple{minchromidx::Int, edges::Int, size::Int}}(undef, l)
    same_min_edgeschromidx = Vector{Bool}(undef, l)

    schmidtmin = Vector{Int}(undef, l)
    schmidtmax = Vector{Int}(undef, l)
    schmidtsame = Vector{Bool}(undef, l)

    minedgerepr = Vector{Graph{Int}}(undef, l)
    minedgerepr_colors = Vector{Vector{Int}}(undef, l)
    minchromidxrepr = Vector{Graph{Int}}(undef, l)
    minchromidxrepr_colors = Vector{Vector{Int}}(undef, l)
    same_min_edgeschromidx_repr = Vector{Bool}(undef, l)

    bipartrankidx_ = Dict([bridx=>Vector{Union{NTuple{bridx,Int},Missing}}(missing, l) for bridx in 6:-1:2]...)

    cardmult = Vector{Union{Missing,Vector{Int}}}(missing, l)

    for (i,row) in enumerate(eachrow(df))
        e = parsetuple(row._minedges)
        minedges[i] = e
        minchromidx[i], same_min_edgeschromidx[i] = (row._minchromidx == "-") ? (e,true) : (parsetuple(row._minchromidx),false)

        if '<' in row._schmidt
            schmidtsame[i] = false
            schmidtmin[i], schmidtmax[i] = parse.(Int, split(row._schmidt,'<'))
        else
            schmidtsame[i] = true
            schmidtmin[i] = schmidtmax[i] = parse.(Int, row._schmidt)
        end

        mer, merc = parsegraphcol(row._minedgerepr, row.vertices)
        minedgerepr[i], minedgerepr_colors[i] = mer, merc
        if row._minchromidxrepr == "-"
            minchromidxrepr[i], minchromidxrepr_colors[i] = mer, merc
            same_min_edgeschromidx_repr[i] = true
        else
            minchromidxrepr[i], minchromidxrepr_colors[i] = parsegraphcol(row._minchromidx, row.vertices)
            same_min_edgeschromidx_repr[i] = false
        end

        for _br in n÷2:-1:2
            bridxsym = Symbol(:_bipartrankidx,_br)
            bipartrankidx_[_br][i] = parsetuple(row[bridxsym])
        end

        if row._cardmult != "-"
            cardmult[i] = parsevec(row._cardmult)
        end
    end

    df[!,:minedges] = minedges
    df[!,:minchromidx] = minchromidx
    df[!,:same_min_edgeschromidx] = same_min_edgeschromidx

    df[!,:schmidtmin] = schmidtmin
    df[!,:schmidtmax] = schmidtmax
    df[!,:schmidtsame] = schmidtsame

    df[!,:minedge_repr] = minedgerepr
    df[!,:minedge_repr_colors] = minedgerepr_colors
    df[!,:same_min_edgeschromidx_repr] = same_min_edgeschromidx_repr
    df[!,:minchromidx_repr] = minchromidxrepr
    df[!,:minchromidx_repr_colors] = minchromidxrepr_colors

    for _br in 12÷2:-1:2
        bridxsym = Symbol(:bipartrankidx,_br)
        df[!,bridxsym] = bipartrankidx_[_br]
    end

    df[!,:cardmult] = cardmult

    return select!(df, Not(:_minedges,:_minchromidx,:_schmidt,bipartranknames(n)...,:_cardmult,:_minedgerepr,:_minchromidxrepr))
end

"""Get the database from `arxiv:1011.5464`

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
"""
function lcorbits_summary end

function lcorbits_summary(n)
    (2≤n≤12) || error("The arxiv:1011.5464 database contains only graph states for a number of qubits 2≤n≤12")
    return clean_up_df_larsed(n)
end

function lcorbits_summary()
    vcat([lcorbits_summary(n) for n in 2:12]...)
end

end
