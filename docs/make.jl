using Documenter

makedocs(;
    sitename="BioMASS",
    pages = [
        "Home" => "index.md",
        "Getting started with BioMASS" => [
            "Prerequisites and Installation" => "introduction/set_up.md",
            "Overview" => "introduction/overview.md",
            "Visualization of simulation results" => "introduction/visualization.md",
            "Usage" => "introduction/usage.md",
        ],
        "Citation" => "citation.md"
    ]
)
