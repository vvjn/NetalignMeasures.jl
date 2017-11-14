using Documenter, NetalignMeasures

makedocs(
           format = :html,
           sitename = "NetalignMeasures",
           pages = [
                    "index.md"
           ],
    modules = [NetalignMeasures]
       )

deploydocs(
           repo = "github.com/vvjn/NetalignMeasures.jl.git",
           target = "build",
           deps   = nothing,
           make   = nothing,
           julia = "0.6"
)
