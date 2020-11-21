using SignalTemps
using Documenter

makedocs(;
    modules=[SignalTemps],
    authors="Matt Karikomi <mattkarikomi@gmail.com> and contributors",
    repo="https://github.com/mkarikom/SignalTemps.jl/blob/{commit}{path}#L{line}",
    sitename="SignalTemps.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mkarikom.github.io/SignalTemps.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mkarikom/SignalTemps.jl",
)
