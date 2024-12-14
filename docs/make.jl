using AJD
using Documenter

DocMeta.setdocmeta!(AJD, :DocTestSetup, :(using AJD); recursive=true)

makedocs(;
    modules=[AJD],
    authors="Nicholas Gericke <nicholas.gericke@campus.tu-berlin.de>",
    sitename="AJD.jl",
    format=Documenter.HTML(;
        canonical="https://gericke-n.github.io/AJD.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
    ],
)

deploydocs(;
    repo="github.com/gericke-n/AJD.jl",
    devbranch="master",
)
