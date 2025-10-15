push!(LOAD_PATH, "src")
push!(LOAD_PATH, "lib/XSPECModels/src")

ENV["GKSwstype"] = "100"

using Documenter
using SpectralFitting, XSPECModels

SpectralFitting.download_all_model_data()

makedocs(
    warnonly = [:cross_references, :autodocs_block, :missing_docs],
    modules = [SpectralFitting, XSPECModels],
    clean = true,
    sitename = "SpectralFitting.jl",
    pages = [
        "Home" => "index.md",
        "Walkthrough" => "walkthrough.md",
        "Examples" => [
            "Diverse examples" => "examples/examples.md",
            "A quick guide" => "examples/sherpa-example.md",
        ],
        "Models" => [
            "Using models" => "models/using-models.md",
            "Model index" => "models/models.md",
            "Composite models" => "models/composite-models.md",
            "Surrogate models" => "models/surrogate-models.md",
            "XSPEC models" => "models/xspec-models.md",
        ],
        "Fitting" => [
            "parameters.md",
            # "fitting.md",
        ],
        "Datasets" => [
            "Using datasets" => "datasets/datasets.md",
            "Grouping & Rebinning" => "datasets/grouping-rebinning.md",
            "Mission support" => "datasets/mission-support.md",
        ],
        "Why & How" => "why-and-how.md",
        "Reference" => "reference.md",
    ],
)

deploydocs(repo = "github.com/fjebaker/SpectralFitting.jl.git")
