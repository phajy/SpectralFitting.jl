export AbstractCompositeOperator,
    AdditionOperator,
    MultiplicationOperator,
    ConvolutionOperator,
    operation_symbol,
    CompositeModel

"""
    abstract type AbstractCompositeOperator

Superype of all composition operators. Used to implement the model algebra of
[`AbstractSpectralModelKind`](@ref) through a trait system.

The Julia symbol corresponding to a given `AbstractCompositeOperator` may be
obtained through [`operation_symbol`](@ref).
"""
abstract type AbstractCompositeOperator end
"""
    AdditionOperator <: AbstractCompositeOperator
    AdditionOperator()

Corresponds to the `:(+)` symbol.
"""
struct AdditionOperator <: AbstractCompositeOperator end
"""
    MultiplicationOperator <: AbstractCompositeOperator
    MultiplicationOperator()

Corresponds to the `:(*)` symbol.
"""
struct MultiplicationOperator <: AbstractCompositeOperator end
"""
    ConvolutionOperator <: AbstractCompositeOperator
    ConvolutionOperator()

Has no corresponding symbol, since it invokes a function call `C(A)`.
"""
struct ConvolutionOperator <: AbstractCompositeOperator end

# define these on type, since mostly used in the generation context
# which only operates on un-instantiated types
"""
    operation_symbol(::AbstractCompositeOperator)
    operation_symbol(::Type{<:AbstractCompositeOperator})

Obtain the model symbol from a given [`AbstractCompositeOperator`](@ref).
"""
operation_symbol(::Type{<:AdditionOperator}) = :(+)
operation_symbol(::Type{<:MultiplicationOperator}) = :(*)
operation_symbol(::Type{<:ConvolutionOperator}) = nothing
operation_symbol(::O) where {O<:AbstractCompositeOperator} = operation_symbol(O)

# algebra grammar
add_models(_, _, ::M1, ::M2) where {M1,M2} =
    throw("Left and right models must be Additive.")
add_models(m1, m2, ::Additive, ::Additive) = CompositeModel(m1, m2, AdditionOperator())
add_models(m1::M1, m2::M2) where {M1,M2} = add_models(m1, m2, modelkind(M1), modelkind(M2))
Base.:+(m1::M1, m2::M2) where {M1<:AbstractSpectralModel,M2<:AbstractSpectralModel} =
    add_models(m1, m2, modelkind(M1), modelkind(M2))

mult_models(_, _, ::M1, ::M2) where {M1,M2} = throw("Left model must be Multiplicative.")
mult_models(m1, m2, ::Multiplicative, ::AbstractSpectralModelKind) =
    CompositeModel(m1, m2, MultiplicationOperator())
mult_models(m1::M1, m2::M2) where {M1,M2} =
    mult_models(m1, m2, modelkind(M1), modelkind(M2))
Base.:*(m1::M1, m2::M2) where {M1<:AbstractSpectralModel,M2<:AbstractSpectralModel} =
    mult_models(m1, m2, modelkind(M1), modelkind(M2))

conv_models(_, _, ::M1, ::M2) where {M1,M2} =
    throw("Left model must be Convolutional and right model must be Additive.")
conv_models(m1, m2, ::Convolutional, ::Additive) =
    CompositeModel(m1, m2, ConvolutionOperator())
conv_models(m1::M1, m2::M2) where {M1,M2} =
    conv_models(m1, m2, modelkind(M1), modelkind(M2))
(m1::AbstractSpectralModel)(m2::M2) where {M2<:AbstractSpectralModel} = conv_models(m1, m2)

"""
    CompositeModel{M1,M2,O} <: AbstractSpectralModel
    CompositeModel(left_model, right_model, op::AbstractCompositeOperator)

Type resulting from operations combining any number of [`AbstractSpectralModel`](@ref) via the
model algebra defined from [`AbstractSpectralModelKind`](@ref).

Each operation binary operation in the model algebra is encoded in the parametric types of the
`CompositeModel`, where the operation is given by an [`AbstractCompositeOperator`](@ref).
Composite models adopt the model kind of the `right` model, i.e. `M2`, and obey the model algebra
accordingly.

Composite models very rarely need to be constructed directly, and are instead obtained by regular
model operations.

# Example

```julia
model = PhotoelectricAbsorption() * (PowerLaw() + BlackBody())
typeof(model) <: CompositeModel # true
```
"""
struct CompositeModel{M1,M2,O,T,K} <: AbstractSpectralModel{T,K}
    left::M1
    right::M2

    CompositeModel{M1,M2,O}(
        m1::M1,
        m2::M2,
    ) where {M1<:AbstractSpectralModel{T},M2<:AbstractSpectralModel{T,K},O} where {T,K} =
        new{M1,M2,O,T,K}(copy(m1), copy(m2))
end

CompositeModel(
    m1::M1,
    m2::M2,
    ::O,
) where {M1<:AbstractSpectralModel,M2<:AbstractSpectralModel,O} =
    CompositeModel{M1,M2,O}(copy(m1), copy(m2))

function Base.copy(m::CompositeModel{M1,M2,O}) where {M1,M2,O}
    CompositeModel{M1,M2,O}(getfield(m, :left), getfield(m, :right))
end

function implementation(::Type{<:CompositeModel{M1,M2}}) where {M1,M2}
    if (implementation(M1) isa XSPECImplementation) ||
       (implementation(M2) isa XSPECImplementation)
        XSPECImplementation()
    else
        JuliaImplementation()
    end
end

# invocation wrappers
function invokemodel(e, m::CompositeModel)
    fluxes = allocate_model_output(m, e)
    invokemodel!(fluxes, e, m)
    view(fluxes, :, 1)
end

function invokemodel!(f, e, model::CompositeModel)
    @assert size(f, 2) == objective_cache_count(model) "Too few flux arrays allocated for this model."
    (model, parameter_tuple(model))
    composite_model_call!(f, e, model, parameter_tuple(model))
end
function invokemodel!(f, e, model::CompositeModel, parameters::AbstractVector)
    @assert size(f, 2) == objective_cache_count(model) "Too few flux arrays allocated for this model."
    composite_model_call!(f, e, model, parameters)
end

function Base.show(io::IO, @nospecialize(model::CompositeModel))
    destructure = destructure_model(model)
    expr = "$(destructure.expression)"

    for (sym, model) in destructure.model_map
        s = "$sym"
        name = "$(Base.typename(typeof(model)).name)"
        expr = replace(expr, s => name)
    end
    print(io, "CompositeModel[")
    printstyled(io, "$(expr)", color = :cyan)
    print(io, "]")
end

function _print_param(io, free, name, val, q0, q1, q2, q3, q4; binding = nothing)
    print(io, lpad("$name", q0), " ->")
    if val isa FitParam
        info = get_info_tuple(val)
        print(io, lpad(info[1], q1 + 1))
        if free
            print(io, " ± ", rpad(info[2], q2))
            print(io, " ∈ [", lpad(info[3], q3), ", ", rpad(info[4], q4), "]")
        end

        if !isnothing(binding)
            print(io, "   ")
            printstyled(io, lpad(binding, 7), color = :magenta)
        elseif free
            printstyled(io, lpad("FREE", 7), color = :green)
        else
            print(io, "  ")
            printstyled(io, lpad("FROZEN", 15 + q1 + q2 + q3 + q4), color = :cyan)
        end
    end
    println(io)
end

function _printinfo(io::IO, @nospecialize(model::CompositeModel); bindings = nothing)
    expr_buffer = 5
    sym_buffer = 5

    destructed = destructure_model(model)

    info_tuples = [get_info_tuple(p) for (_, p) in destructed.parameter_map]
    q1, q2, q3, q4 = map(1:4) do i
        maximum(j -> length("$(j[i])"), info_tuples) + 1
    end

    param_offset = sym_buffer + maximum(destructed.parameter_map) do (s, _)
        length("$s")
    end

    println(io, "CompositeModel with $(length(destructed.model_map)) model components:")
    print(io, " "^expr_buffer)
    printstyled(io, destructed.expression, color = :cyan)
    println(io)
    println(io, "Model key and parameters:")

    param_index = 1
    for (sym, m) in destructed.model_map
        param_syms = destructed.parameter_symbols[sym]
        basename = Base.typename(typeof(m)).name

        printstyled(io, lpad("$sym", sym_buffer), color = :cyan)
        print(io, " => ")
        printstyled(io, "$basename", color = :cyan)
        println(io)

        for ps in param_syms
            param = destructed.parameter_map[ps]
            free = param isa FitParam ? !isfrozen(param) : true
            val, binding = if !isnothing(bindings) && !isempty(bindings)
                get(bindings, param_index, param => nothing)
            else
                param, nothing
            end
            _print_param(io, free, ps, val, param_offset, q1, q2, q3, q4; binding)
            param_index += 1
        end
    end
end

# explicitly disable interface for ConstructionBase.jl
ConstructionBase.setproperties(::CompositeModel, ::NamedTuple) =
    throw("Cannot be used with `CompositeModel`.")
ConstructionBase.constructorof(::Type{<:CompositeModel}) =
    throw("Cannot be used with `CompositeModel`.")

@inline @generated function _composite_parameter_symbols(model::CompositeModel)
    syms = [:($(Meta.quot(s))) for s in Reflection.get_parameter_symbols(model)]
    :(($(syms...),))
end

@inline @generated function _composite_model_symbols(model::CompositeModel)
    info = Reflection.get_info(model, :model; T = Float64)
    syms = [:($(Meta.quot(s))) for (s, _) in info.models]
    :(($(syms...),))
end

@inline @generated function _composite_model_map(model::CompositeModel)
    info = Reflection.get_info(model, :model; T = Float64)
    syms = ([:($s) for (s, _) in info.models]...,)
    values = [i.lens for (_, i) in info.models]
    :(NamedTuple{$(syms)}(($(values...),)))
end

function Base.propertynames(model::CompositeModel)
    (_composite_parameter_symbols(model)..., _composite_model_symbols(model)...)
end

# TODO: really ensure this is type stable as it could be a performance killer
Base.@constprop aggressive function _get_property(model::CompositeModel, symb::Symbol)
    if symb in _composite_model_symbols(model)
        return _composite_model_map(model)[symb]
    else
        return parameter_named_tuple(model)[symb]
    end
end

function Base.getproperty(model::CompositeModel, symb::Symbol)
    _get_property(model, symb)
end

function Base.setproperty!(model::CompositeModel, symb::Symbol, value::FitParam)
    set!(getproperty(model, symb), value)
end

function Base.setproperty!(model::CompositeModel, symb::Symbol, x)
    error(
        "Only `FitParam` may be directly set with another `FitParam`. Use `set_value!` and related API otherwise.",
    )
end
