using Documenter

makedocs(;
    sitename="BioMASS",
    pages = [
        "Home" => "index.md",
        "Getting started with BioMASS" => [
            "Installation" => "introduction/set_up.md",
            "Using BioMASS" => "introduction/usage.md",
        ],
        "Citation" => "citation.md"
    ]
)