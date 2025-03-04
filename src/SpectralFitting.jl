module SpectralFitting

import Base
import Printf
import Downloads
import Statistics
import Pkg.MiniProgressBars: MiniProgressBar, start_progress, end_progress, show_progress
import Distributions
import ConstructionBase

import Random

using FITSIO
using SparseArrays
using Surrogates
using LinearAlgebra
using FileIO
using Interpolations
import DataInterpolations
using SpecialFunctions
using Polynomials
using PreallocationTools
using EnumX

# fitting backends
import ForwardDiff
import LsqFit
import Optimization

using DocStringExtensions

abstract type AbstractInstrument end

abstract type AbstractStatistic end
struct ChiSquared <: AbstractStatistic end
struct Cash <: AbstractStatistic end

# unitful units
include("units.jl")
SpectralUnits.@reexport using .SpectralUnits

include("utils.jl")
include("print-utilities.jl")
include("support.jl")

include("fitparam.jl")
include("param-cache.jl")
include("abstract-models.jl")

include("composite-models.jl")

include("reflection.jl")

include("meta-models/wrappers.jl")
include("meta-models/table-models.jl")
include("meta-models/surrogate-models.jl")
include("meta-models/caching.jl")
include("meta-models/functions.jl")

include("datasets/ogip.jl")
include("datasets/datasets.jl")
include("datasets/binning.jl")
include("datasets/grouping.jl")
include("datasets/injectivedata.jl")
include("datasets/binneddata.jl")

include("model-data-io.jl")

# include fitting api
include("fitting/problem.jl")
include("fitting/cache.jl")
include("fitting/config.jl")
include("fitting/result.jl")
include("fitting/statistics.jl")
include("fitting/binding.jl")
include("fitting/methods.jl")

include("simulate.jl")
include("fitting/goodness.jl")

# include julia models
include("julia-models/model-utilities.jl")
include("julia-models/additive.jl")
include("julia-models/multiplicative.jl")
include("julia-models/convolutional.jl")

include("plots-recipes.jl")

function __init__()
    # check if we have the minimum model data already
    _check_model_directory_present()
end

end # module
