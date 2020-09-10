using Documenter

makedocs(;
    sitename="BioMASS",
    pages = [
        "Home" => "index.md",
        "Getting started with BioMASS" => [
            "Installation" => "introduction/set_up.md",
            "Import model" => "introduction/import_model.md",
            "Parameter estimation" => "introduction/parameter_estimation.md",
            "Visualization of simulation results" => "introduction/visualization.md",
            "Sensitivity analysis" => "introduction/sensitivity_analysis.md",
        ]
    ]
)