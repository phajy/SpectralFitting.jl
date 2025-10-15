mutable struct Spectrum{T} <: AbstractDataset
    channels::Vector{Int}
    quality::Vector{Int}
    grouping::Vector{Int}

    # this could be counts or flux
    data::Vector{T}
    units::SpectralUnits.RateOrCount

    exposure_time::T
    background_scale::T
    area_scale::T

    error_statistics::SpectralFitting.ErrorStatistics.T
    errors::Union{Nothing,Vector{T}}
    systematic_error::T

    telescope_name::String
    instrument::String
end

function normalize!(spectrum::Spectrum)
    if spectrum.units == u"counts"
        @. spectrum.data /= spectrum.exposure_time
        @. spectrum.errors /= spectrum.exposure_time
        spectrum.units = u"counts/ s"
    end
    spectrum
end

supports(::Type{<:Spectrum}) = (ContiguouslyBinned(),)

function make_objective(::ContiguouslyBinned, dataset::Spectrum)
    dataset.data
end

function make_objective_variance(::ContiguouslyBinned, dataset::Spectrum)
    dataset.errors .^ 2
end

function make_model_domain(::ContiguouslyBinned, dataset::Spectrum)
    @warn "Spectrum doesn't know the energy values by default. Domain is channels. Proceed only if you know what you are doing."
    dataset.channels
end

function make_output_domain(::ContiguouslyBinned, dataset::Spectrum)
    @warn "Spectrum doesn't know the energy values by default. Domain is channels. Proceed only if you know what you are doing."
    dataset.channels
end

isgrouped(spectrum::Spectrum) = all(==(1), spectrum.grouping)

regroup!(spectrum::Spectrum) = regroup!(spectrum, spectrum.grouping)
function regroup!(spectrum::Spectrum{T}, grouping) where {T}
    if all(==(1), grouping)
        # no op
        return spectrum
    end

    itt = GroupingIterator(grouping)
    for grp in itt
        spectrum.channels[grp[1]] = grp[1]
        regroup_vector!(spectrum.data, grp)
        regroup_quality_vector!(spectrum.quality, grp)
        if !isnothing(spectrum.errors)
            vs = spectrum.data[grp[1]]
            if spectrum.units == u"counts"
                spectrum.errors[grp[1]] = count_error(vs, 1.0)
            elseif spectrum.units == u"counts / s"
                vs = vs * spectrum.exposure_time
                es = count_error(vs, 1.0)
                spectrum.errors[grp[1]] = es / spectrum.exposure_time
            else
                error(
                    "No method for grouping errors with given spectral units ($(spectrum.units)).",
                )
            end
        end
    end

    resize!(spectrum, length(itt))
    spectrum.grouping .= 1
    spectrum
end

function group_min_counts!(spectrum::Spectrum, min_counts::Int)
    NEW_GRP = 1
    CONTINUE_GRP = 0

    function _counts(x)
        if spectrum.units == u"counts"
            convert(Int, x)
        elseif spectrum.units == u"counts / s"
            convert(Int, x * spectrum.exposure_time)
        end
    end

    sum::Int = 0
    for (i, f) in enumerate(spectrum.data)
        c = _counts(f)
        sum += c
        if sum >= min_counts
            spectrum.grouping[i] = NEW_GRP
            sum = 0
        else
            spectrum.grouping[i] = CONTINUE_GRP
        end
    end
end

"""
    group_min_snr!(spectrum::Spectrum, min_snr::Real; background::Union{Nothing,Spectrum} = nothing)

Group spectrum channels to achieve a minimum signal-to-noise ratio (S/N) per bin.

The signal-to-noise ratio is calculated using the formula from XMM SAS specgroup:

    S/N = (source - background × areanorm) / √(source + background × areanorm²)

where `areanorm = (backscal_source / backscal_background) × (exposure_source / exposure_background)`.

If no background is provided, the S/N is calculated as:

    S/N = source / √source

# Arguments
- `spectrum::Spectrum`: The spectrum to be grouped (modified in-place)
- `min_snr::Real`: Minimum signal-to-noise ratio required per bin
- `background::Union{Nothing,Spectrum}`: Optional background spectrum. If provided, 
  background subtraction is accounted for in the S/N calculation.

# Notes
- Both spectrum and background must have units of counts (not count rate)
- The grouping is applied sequentially from low to high channels
- Bad quality channels are included in the grouping but may affect the S/N calculation

# Example
```julia
# Group to minimum S/N of 5 without background
group_min_snr!(spec, 5.0)

# Group to minimum S/N of 3 with background
group_min_snr!(spec, 3.0; background=bgd_spec)
```
"""
function group_min_snr!(
    spectrum::Spectrum,
    min_snr::Real;
    background::Union{Nothing,Spectrum} = nothing,
)
    NEW_GRP = 1
    CONTINUE_GRP = 0

    # Helper function to convert to counts
    function _counts(x, exposure, units)
        if units == u"counts"
            x
        elseif units == u"counts / s"
            x * exposure
        else
            error("Unsupported units: $(units)")
        end
    end

    # Calculate area normalization if background is provided
    areanorm = if !isnothing(background)
        (spectrum.background_scale / background.background_scale) *
        (spectrum.exposure_time / background.exposure_time)
    else
        0.0
    end

    # Initialize accumulators
    source_sum = 0.0
    background_sum = 0.0

    for (i, source_data) in enumerate(spectrum.data)
        # Convert to counts for the calculation
        source_counts = _counts(source_data, spectrum.exposure_time, spectrum.units)
        source_sum += source_counts
        
        if !isnothing(background)
            background_counts = _counts(background.data[i], background.exposure_time, background.units)
            background_sum += background_counts
        end
        
        # Calculate S/N for current accumulated bin
        snr = if !isnothing(background)
            # With background: S/N = (S - B*areanorm) / sqrt(S + B*areanorm^2)
            signal = source_sum - background_sum * areanorm
            noise_sq = source_sum + background_sum * (areanorm^2)
            noise = sqrt(max(noise_sq, 0.0))  # Ensure non-negative
            noise > 0.0 ? signal / noise : 0.0
        else
            # Without background: S/N = S / sqrt(S)
            source_sum > 0.0 ? source_sum / sqrt(source_sum) : 0.0
        end
        
        # If we've met the S/N threshold, close this group and prepare to start a new one
        if snr >= min_snr
            spectrum.grouping[i] = NEW_GRP
            source_sum = 0.0
            background_sum = 0.0
        else
            spectrum.grouping[i] = CONTINUE_GRP
        end
    end

    spectrum
end

function Base.resize!(spectrum::Spectrum, n::Int)
    resize!(spectrum.channels, n)
    resize!(spectrum.data, n)
    resize!(spectrum.grouping, n)
    resize!(spectrum.quality, n)
    if !isnothing(spectrum.errors)
        resize!(spectrum.errors, n)
    end
end

function drop_channels!(spectrum::Spectrum, indices)
    deleteat!(spectrum.channels, indices)
    deleteat!(spectrum.data, indices)
    deleteat!(spectrum.quality, indices)
    deleteat!(spectrum.grouping, indices)
    if !isnothing(spectrum.errors)
        deleteat!(spectrum.errors, indices)
    end
    length(indices)
end

_readable_boolean(b) = b ? "yes" : "no"

function _printinfo(io::IO, spectrum::Spectrum)
    dmin, dmax = prettyfloat.(extrema(spectrum.data))
    is_grouped = isgrouped(spectrum) |> _readable_boolean
    num_bad = count(!=(GOOD_QUALITY), spectrum.quality)
    has_bad = num_bad > 0 ? "yes ($num_bad)" : "no"
    descr = """Spectrum: $(spectrum.telescope_name)[$(spectrum.instrument)]
      Units                 : $(spectrum.units)
      . Exposure time       : $(spectrum.exposure_time)
      . Channels            : $(length(spectrum.channels))
      . Data (min/max)      : ($dmin, $dmax)
      . Grouped             : $is_grouped
      . Bad channels        : $has_bad
    """
    print(io, descr)
end

error_statistic(spec::Spectrum) = spec.error_statistics

function subtract_background!(spectrum::Spectrum, background::Spectrum)
    @assert spectrum.units == u"counts"
    # errors added in quadrature
    # TODO: this needs fixing to propagate errors properly
    data_variance = spectrum.errors .^ 2
    background_variance = background.errors .^ 2
    _subtract_background!(
        spectrum.errors,
        data_variance,
        background_variance,
        spectrum.area_scale,
        background.area_scale,
        spectrum.background_scale,
        background.background_scale,
        spectrum.exposure_time,
        background.exposure_time,
    )
    @. spectrum.errors = √abs(spectrum.errors)
    _subtract_background!(
        spectrum.data,
        spectrum.data,
        background.data,
        spectrum.area_scale,
        background.area_scale,
        spectrum.background_scale,
        background.background_scale,
        spectrum.exposure_time,
        background.exposure_time,
    )
    spectrum
end

"""
Does the background subtraction and returns units of counts. That means we have
multiplied through by a factor ``t_D`` relative to the reference equation (2.3)
in the XSPEC manual.
"""
_subtract_background!(output, spec, back, aD, aB, bD, bB, tD, tB) =
    @. output = (spec / (aD)) - (tD / tB) * _scaled_background(back, aB, bD, bB)

_scaled_background(back, aB, bD, bB) = (bD / bB) * (back / aB)


export Spectrum
