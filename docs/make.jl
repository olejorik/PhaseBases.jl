using PhaseBases
using Documenter, Literate

DocMeta.setdocmeta!(PhaseBases, :DocTestSetup, :(using PhaseBases); recursive=true)

# ---- Convert example scripts to Markdown via Literate.jl ----
@info "Current dir = $(@__DIR__)"
tutorials_src = joinpath(@__DIR__, "..", "examples")
tutorials_dst = joinpath(@__DIR__, "src", "examples")
mkpath(tutorials_dst)
for f in readdir(tutorials_src; join=true)
    endswith(f, ".jl") || continue
    Literate.markdown(f, tutorials_dst)
end

makedocs(;
    sitename = "PhaseBases.jl",
    modules  = [PhaseBases],
    authors  = "Oleg Soloviev",
    repo     = "https://github.com/olejorik/PhaseBases.jl/blob/{commit}{path}#L{line}",
    checkdocs = :exports,
    format   = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://olejorik.github.io/PhaseBases.jl/stable/",
        assets     = String[],
    ),
    clean = false,
    pages = [
        "Home"  => "index.md",
        "About" => "about.md",
        "Examples" => [
            "Symbolic Zernike Polynomials" => "examples/SymbolicZernike.md",
            "Zernike Pyramids (grid layout)"  => "examples/ZernikePyramids.md",
            "Zernike Pyramids (BBox layout)"  => "examples/ZernikePyramids2.md",
        ],
        "Reference" => ["API" => "api.md"],
    ],
)

deploydocs(;
    repo       = "github.com/olejorik/PhaseBases.jl",
    target     = "build",
    devbranch  = "main",
)

