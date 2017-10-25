using Documenter, NetalignMeasures

makedocs()

deploydocs(
           deps   = Deps.pip("mkdocs", "python-markdown-math"),           
           repo = "github.com/vvjn/NetalignMeasures.jl.git",
           julia = "0.6"
)

# makedocs(
#     format = :html,
#     sitename = "NetalignUtils",
#     pages = [
#         "index.md"
#     ]
# )

