@enumx ErrorStatistics begin
    Numeric
    Poisson
    Gaussian
    Unknown
end

abstract type AbstractLayout end
struct OneToOne <: AbstractLayout end
struct ContiguouslyBinned <: AbstractLayout end

export OneToOne, ContiguouslyBinned

# used to check what to default to
# todo: check if can be compile time eval'd else expand for loop or @generated
function common_support(x, y)
    # order of preference is important
    # can normally trivially fallback from one-to-one to contiguous bins to regular bins
    for layout in (ContiguouslyBinned(), OneToOne())
        if supports(layout, x) && supports(layout, y)
            return layout
        end
    end
    error("No common support between $(typeof(x)) and $(typeof(y)).")
end

function _support_reducer(x::OneToOne, y)
    if supports(x, y)
        return x
    else
        error("No common support!!")
    end
end
function _support_reducer(x::ContiguouslyBinned, y)
    if supports(x, y)
        return x
    else
        _support_reducer(OneToOne(), y)
    end
end
function _support_reducer(x, y)
    common_support(x, y)
end

common_support(args::Vararg) = reduce(_support_reducer, args)

supports_contiguosly_binned(::Type) = false
supports_one_to_one(::Type) = false

supports_contiguosly_binned(x) = supports_contiguosly_binned(typeof(x))
supports_one_to_one(x) = supports_one_to_one(typeof(x))

supports(layout::AbstractLayout, x) = supports(layout, typeof(x))
supports(::ContiguouslyBinned, T::Type) = supports_contiguosly_binned(T)
supports(::OneToOne, T::Type) = supports_one_to_one(T)

"""
we make special exceptions if only a single support is given

above is used for both models and
datasets to address compatability. if the models and data have different support,
they are incompatible

Cases to address:

x to y

contiguously binned x to y

if not contiguously binned can split into multiple seperate datasets

- things like background subtraction should happend _before_ it is wrapped
simplifies everything conceptually
- might need to have some kind of `prepare!` api that does all preparation before a fit is done
but an issue with that is that it would modify the underlying data representation without
the user's knowledge. can avoid that via duplication, but that seems excessive
- only thing is though is if that is called internally, that then forces the user to have the
data in a specific format during fit. maybe they don't want background subtraction
- think it's good to keep a `prepare!` API that does all the sensible default things,
but that it isn't called by the fitting routines
- this is then hopefully minimal boilerplate. we just add a whole load of warnings

how do we handle things like response folding or background subtraction in fit?
- i think we need a `apply_objective_transformation!(y, x, dataset)` function for this, that
does any folding or whatever else. the difficulty here would be communicating back
what format the data is in. ideally need some kind of bi-directional transformations
- 

multiple distinct data sets?
- for this we just have many AbstractDatasets
i think it's probably easiest to just have one "underlying" dataset being wrapped in this
wrapping structure, and to handle multiplexing seperately
"""


"""
    abstract type AbstractDataset
    
Abstract type for use in fitting routines. High level representation of some underlying 
data structures. 

Fitting data is considered to have an *objective* and a *domain*. As
the domain may be, for example, energy bins (high and low), or
fourier frequencies (single value), the purpose of this abstraction
is to provide some facility for translating between these representations
for the models to fit with. This is done by checking that the [`AbstractLayout`](@ref)
of the model and data are compatible, or at least have compatible translations.

Must implement a minimal set of accessor methods. These are paired with
`objective` and `domain` parlance. Note that these functions are prefixed with
`make_*` and not `get_*` to represent that there may be allocations or work going
into the translation. Usage of these functions should be sparse in the interest of
performance.

The arrays returned by the `make_*` functions must correspond to the [`AbstractLayout`](@ref)
specified by the caller.

[`ContiguouslyBinned`](@ref) (domain should be `length(objective) + 1`, where the limits of the
``n^\text{th}`` bin are `domain[n]` and `domain[n+1]` respectively.

- [`make_objective_variance`](@ref)
- [`make_objective`](@ref)
- [`make_domain_variance`](@ref)
- [`make_domain`](@ref)

"""
abstract type AbstractDataset end

function Base.show(io::IO, ::MIME"text/plain", data::AbstractDataset)
    buff = IOBuffer()
    _printinfo(buff, data)
    s = String(take!(buff))
    print(io, encapsulate(s))
end


"""
    make_objective

Returns the array used as the target for model fitting. The array must correspond to the data
[`AbstractLayout`](@ref) specified by the `layout` parameter.

In as far as it can be guarunteed, the memory in the returned array will not be mutated by any fitting procedures.

Domain for this objective should be returned by [`make_domain`](@ref).
"""
make_objective(layout::AbstractLayout, dataset::AbstractDataset) =
    error("Layout $(layout) is not implemented for $(typeof(dataset))")
make_objective_variance(layout::AbstractLayout, dataset::AbstractDataset) =
    error("Layout $(layout) is not implemented for $(typeof(dataset))")

"""
    make_domain

Returns the array used as the domain for the modelling
"""
make_domain(layout::AbstractLayout, dataset::AbstractDataset) =
    error("Layout $(layout) is not implemented for $(typeof(dataset))")
make_domain_variance(layout::AbstractLayout, dataset::AbstractDataset) =
    error("Layout $(layout) is not implemented for $(typeof(dataset))")

objective_transformer(layout::AbstractLayout, dataset::AbstractDataset) =
    error("Not implemented for $(layout) and $(typeof(dataset))")

"""
Must support the same API, but may also have some query methods for specific internals.
"""
abstract type AbstractMultiDataset <: AbstractDataset end

export make_domain,
    make_domain_variance,
    make_objective,
    make_objective_variance,
    normalize!,
    objective_transformer

include("spectrum.jl")
include("response.jl")
include("spectral-dataset.jl")

# mission specifics
include("xmm-newton.jl")
