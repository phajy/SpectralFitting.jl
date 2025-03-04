using Test, SpectralFitting, XSPECModels

testdir = get(
    ENV,
    "SF_TEST_SUITE_DATA",
    @__DIR__() * "/../../../spectral-fitting-test-suite/sample-data",
)
@show testdir

has_test_dir = isdir(testdir)
if !has_test_dir
    @warn "No test data found. Skipping some tests."
end
if !has_test_dir && get(ENV, "CI", false)
    error("Missing test dir on CI")
end

include("utils.jl")

@testset "reflection" verbose = true begin
    include("reflection/test-reflection.jl")
end

@testset "api" verbose = true begin
    include("models/test-model-api.jl")
end

@testset "macro" verbose = true begin
    @testset "xspecmodel" begin
        include("macros/test-xspec.jl")
    end
end

@testset "composite-algebra" verbose = true begin
    @testset "model-algebra" begin
        include("composite/test-algebra.jl")
    end
    @testset "model-invocation" begin
        include("composite/test-invocation.jl")
    end
end

@testset "parameters" verbose = true begin
    @testset "model-parameters" begin
        include("parameters/test-model-parameters.jl")
    end
    @testset "fit-params" begin
        include("parameters/test-free-frozen.jl")
    end
end

@testset "model-library" verbose = true begin
    include("models/test-julia-models.jl")
    include("models/test-general-models.jl")
    include("models/test-model-consistency.jl")
    include("models/test-table-models.jl")
    include("models/test-surrogate-models.jl")
    include("models/test-auto-cache.jl")
    include("models/test-as-convolution.jl")
    include("models/test-copy.jl")

    # only test XSPEC models when not using CI
    # since model data access is annoying
    @ciskip @testset "xspec-models" begin
        include("models/test-xspec-models.jl")
        include("models/test-general-xspec-models.jl")
    end
end

@testset "io" verbose = true begin
    @testset "printing" begin
        include("io/test-printing.jl")
    end
    @testset "datasets" begin
        if has_test_dir
            include("datasets/test-ogip.jl")
        else
            @warn "Skipping OGIP dataset tests."
        end
        include("datasets/test-grouping.jl")
        include("datasets/test-units.jl")
        include("datasets/test-binning.jl")
        include("datasets/test-datasets.jl")
    end
    include("io/test-remote-pathname-compression.jl")
end

@testset "fitting" verbose = true begin
    @time include("fitting/test-fit-simple-dataset.jl")
    @time include("fitting/test-cache.jl")
    @time include("fitting/test-binding.jl")
    @time include("fitting/test-results.jl")

    @time include("fitting/test-fit-powerlaw.jl")
    @time include("fitting/test-models.jl")

    @testset "multifits" begin
        @time include("fitting/test-fit-multi.jl")
        @time include("fitting/test-fit-optim.jl")
    end
    if has_test_dir
        @testset "sample-data" begin
            @time include("fitting/test-sample-data.jl")
        end
    else
        @warn "Skipping dataset tests."
    end
end

@testset "simulation" verbose = true begin
    include("simulation/test-simulation.jl")
    if has_test_dir
        @testset "sample-data" begin
            include("simulation/test-sample-data-sim.jl")
        end
    else
        @warn "Skipping simulating real observatory tests."
    end
end

using Aqua

# little bit of aqua
Aqua.test_undefined_exports(SpectralFitting)
Aqua.test_unbound_args(SpectralFitting)
