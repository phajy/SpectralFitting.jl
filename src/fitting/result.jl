export FittingResult,
    MultiFittingResult,
    AbstractFittingResult,
    FittingResultSlice,
    invoke_result,
    update_model!

function _pretty_print_result(model, u, σ, chi2)
    ppx2 = prettyfloat(chi2)
    ppu = join((prettyfloat(i) for i in u), ", ")
    ppσ = isnothing(σ) ? nothing : join((prettyfloat(i) for i in σ), ", ")
    """
      Model: $(model)
      . u     : [$(ppu)]
      . σᵤ    : [$(ppσ)]
      . χ²    : $(ppx2)
    """
end

abstract type AbstractFittingResult end

invoke_result(res::AbstractFittingResult) = invoke_result(res, res.u)

struct FittingResultSlice{P<:AbstractFittingResult,V,U,T} <: AbstractFittingResult
    index::Int
    parent::P
    domain::V
    objective::V
    variance::V
    u::U
    σu::Union{Nothing,U}
    χ2::T
end

get_cache(f::FittingResultSlice) = f.parent.config.cache
get_model(f::FittingResultSlice) = f.parent.config.prob.model.m[f.index]
get_dataset(f::FittingResultSlice) = f.parent.config.prob.data.d[f.index]
fit_statistic(f::FittingResultSlice) = fit_statistic(f.parent.config)

estimated_error(r::FittingResultSlice) = r.σu
estimated_params(r::FittingResultSlice) = r.u

function invoke_result(slice::FittingResultSlice{P}, u) where {P}
    @assert length(u) == length(slice.u)
    cache = if P <: MultiFittingResult
        get_cache(slice).caches[slice.index]
    else
        get_cache(slice)
    end
    _invoke_and_transform!(cache, slice.domain, u)
end

function _pretty_print(slice::FittingResultSlice)
    "FittingResultSlice:\n" *
    _pretty_print_result(get_model(slice), slice.u, slice.σu, slice.χ2)
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(slice::FittingResultSlice))
    print(io, encapsulate(_pretty_print(slice)))
end

struct FittingResult{T,U,C} <: AbstractFittingResult
    χ2::T
    u::U
    σu::Union{Nothing,U}
    config::C
end

estimated_error(r::FittingResult) = r.σu
estimated_params(r::FittingResult) = r.u

function invoke_result(result::FittingResult, u)
    @assert length(u) == length(result.u)
    _invoke_and_transform!(result.config.cache, result.config.model_domain, u)
end

function Base.getindex(result::FittingResult, i)
    if i == 1
        @views FittingResultSlice(
            1,
            result,
            result.config.model_domain[:],
            result.config.objective[:],
            result.config.variance[:],
            result.u[:],
            isnothing(result.σu) ? nothing : result.σu[:],
            result.χ2,
        )
    else
        throw(BoundsError())
    end
end

function _pretty_print(res::FittingResult)
    "FittingResult:\n" * _pretty_print_result(res.config.cache.model, res.u, res.σu, res.χ2)
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(res::FittingResult))
    print(io, encapsulate(_pretty_print(res)))
end

struct MultiFittingResult{T,U,C} <: AbstractFittingResult
    χ2s::Vector{T}
    us::U
    σus::Union{Nothing,U}
    config::C
end

estimated_error(r::MultiFittingResult) = r.σus
estimated_params(r::MultiFittingResult) = r.us

function Base.getindex(result::MultiFittingResult, i::Int)
    u = result.us[i]
    σu = isnothing(result.σus) ? nothing : result.σus[i]
    chi2 = result.χ2s[i]
    d_start, d_end = _get_range(result.config.cache.domain_mapping, i)
    o_start, o_end = _get_range(result.config.cache.objective_mapping, i)

    @views FittingResultSlice(
        i,
        result,
        result.config.model_domain[d_start:d_end],
        result.config.objective[o_start:o_end],
        result.config.variance[o_start:o_end],
        u,
        σu,
        chi2,
    )
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(res::MultiFittingResult))
    total_χ2 = prettyfloat(sum(res.χ2s))

    buff = IOBuffer()
    buff_c = IOContext(buff, io)

    println(buff_c, "MultiFittingResult:")
    print(buff_c, " ")
    for i = 1:length(res.us)
        slice = res[i]
        b = _pretty_print_result(get_model(slice), slice.u, slice.σu, slice.χ2)
        r = indent(b, 1)
        print(buff_c, r)
    end
    text = String(take!(buff))
    print(io, encapsulate(text) * "Σχ² = $(total_χ2)")
end

function update_model!(
    model::AbstractSpectralModel,
    result::Union{<:FittingResult,<:FittingResultSlice},
)
    free_params = filter(isfree, parameter_tuple(model))
    for (i, f) in enumerate(free_params)
        set_value!(f, result.u[i])
    end
    model
end

function update_model!(
    model::FittableMultiModel,
    result::Union{<:FittingResult,<:FittingResultSlice},
)
    @assert length(model.m) == 1
    update_model!(model.m[1], result)
end

function update_model!(multimodel::FittableMultiModel, result::MultiFittingResult)
    error("TODO")
end

function finalize(config::FittingConfig, params, final_stat, ; σparams = nothing)
    FittingResult(final_stat, params, σparams, config)
end

function finalize(
    config::FittingConfig{Impl,<:MultiModelCache},
    params,
    final_stat,
    ;
    σparams = nothing,
) where {Impl}
    domain = config.model_domain
    cache = config.cache
    statistic = fit_statistic(config)
    results = map(enumerate(cache.caches)) do (i, ch)
        p = @views params[cache.parameter_mapping[i]]
        σp = @views isnothing(σparams) ? nothing : σparams[cache.parameter_mapping[i]]

        domain_start, domain_end = _get_range(cache.domain_mapping, i)
        objective_start, objective_end = _get_range(cache.objective_mapping, i)

        d = @views domain[domain_start:domain_end]

        output = _invoke_and_transform!(ch, d, p)

        chi2 = measure(
            statistic,
            config.objective[objective_start:objective_end],
            output,
            config.variance[objective_start:objective_end],
        )
        (; chi2, p, σp)
    end

    unc = getindex.(results, :σp)
    unc_or_nothing = if any(isnothing, unc)
        nothing
    else
        unc
    end
    MultiFittingResult(
        getindex.(results, :chi2),
        getindex.(results, :p),
        unc_or_nothing,
        config,
    )
end

function determine_layout(result::FittingResultSlice)
    dataset = get_dataset(result)
    with_units(
        common_support(get_model(result), dataset),
        preferred_units(dataset, fit_statistic(result)),
    )
end

function residuals(result::FittingResultSlice)
    y = invoke_result(result, result.u)
    y_residual = @. (result.objective - y) / sqrt(result.variance)
    y_residual
end
residuals(result::FittingResult; kwargs...) = residuals(result[1]; kwargs...)
