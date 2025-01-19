using AJD
using Documenter

DocMeta.setdocmeta!(AJD, :DocTestSetup, :(using AJD); recursive=true)

makedocs(;
    modules=[AJD],
    authors="Nicholas Gericke <nicholas.gericke@campus.tu-berlin.de>",
    sitename="AJD.jl",
    format=Documenter.HTML(;
        canonical="https://muehlefeldt.github.io/AJD.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting-started.md",
        "Detailed Docs" => "detailed-docs.md",
        "Theory" => "theoretical-background.md",
        "References" => "references.md",
    ],
)

deploydocs(;
    repo="github.com/muehlefeldt/AJD.jl",
    devbranch="master",
)
