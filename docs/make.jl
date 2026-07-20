using Documenter
using LCOrbits

makedocs(
    sitename = "LCOrbits.jl",
    modules = [LCOrbits],
    doctest = false,
    warnonly = :missing_docs,
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/QuantumSavory/LCOrbits.jl.git",
    devbranch = "master",
    push_preview = true,
)
