using PhaseBases
using Documenter

makedocs(;
    modules=[PhaseBases],
    authors="Oleg Soloviev",
    repo="https://github.com/olejorik/PhaseBases.jl/blob/{commit}{path}#L{line}",
    sitename="PhaseBases.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://olejorik.github.io/PhaseBases.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/olejorik/PhaseBases.jl",
)
