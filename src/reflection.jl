
module Reflection

import ..SpectralFitting:
    AbstractSpectralModel,
    CompositeModel,
    FitParam,
    Convolutional,
    Multiplicative,
    Additive,
    operation_symbol,
    AbstractSpectralModelKind,
    modelkind

"""
    const Lens = Union{Symbol,Expr}

Used to create a lens for accessing a field of a given struct. For example,
a lens down to a field `:a` of a nested struct might be something like

```julia
getfield(getfield(instance, :b), :a)
```
"""
const Lens = Union{Symbol,Expr}
const DEFAULT_LENS::Lens = :(model)
const DOMAIN_SYMBOL::Symbol = :(domain)
const OBJECTIVE_SYMBOL::Symbol = :(objective)

const MODEL_SYMBOL_LOOKUP = Dict{AbstractSpectralModelKind,Symbol}(
    Convolutional() => :c,
    Multiplicative() => :m,
    Additive() => :a,
)

"""
    struct ModelInfo
        "All parameter symbols for the model."
        symbols::Vector{Symbol}
        "Unique symbols generated for the parameter assignment when buildin the function call."
        generated_symbols::Vector{Symbol}
        "Additional closure parameters that need to be handled when invoking the model."
        closure_symbols::Vector{Symbol}
        "Unique closure generated symbols."
        generated_closure_symbols::Vector{Symbol}
        "The lens that takes you to this model from some parent."
        lens::Lens
        "The model type itself."
        model::Type
    end

All models are parsed into a `ModelInfo` struct relative to their parent (in the
case of composite models).

The symbols field contains all of the model parameter symbols _as they are in
the structure_, not as they have been generated. Recall when the invocation
expressions are generated, we create anonymous paramter names to avoid conflicts.
These are the `generated_symbols` instead.
"""
struct ModelInfo
    "All parameter symbols for the model."
    symbols::Vector{Symbol}
    "Unique symbols generated for the parameter assignment when buildin the function call."
    generated_symbols::Vector{Symbol}
    "Additional closure parameters that need to be handled when invoking the model."
    closure_symbols::Vector{Symbol}
    "Unique closure generated symbols."
    generated_closure_symbols::Vector{Symbol}
    "The lens that takes you to this model from some parent."
    lens::Lens
    "The model type itself."
    model::Type
end

function make_invoke_call(info::ModelInfo, objective::Symbol, domain::Symbol, T::Type)
    constructor = make_constructor(
        info.model,
        info.generated_closure_symbols,
        info.generated_symbols,
        T,
    )
    :(invokemodel!($objective, $domain, $(constructor)))
end

function reassemble_model(model::Type{<:AbstractSpectralModel}, parameters)
    T = eltype(parameters)
    info = get_info(model, DEFAULT_LENS)
    closure_assigns = Expr[]
    for c_lens in closure_parameter_lenses(info)
        push!(closure_assigns, c_lens)
    end

    params = map(i -> :(parameters[$i]), eachindex(info.symbols))
    make_constructor(info.model, closure_assigns, params, T)
end

"""
    struct CompositeModelInfo
        "The parameter symbols of the model with the respective lens to the actual parameter."
        parameter_symbols::Vector{Pair{Symbol,Lens}}
        "Each model assigned to a unique symbol."
        models::Vector{Pair{Symbol,ModelInfo}}
        "The expression representing the folding operations of this composite model."
        model_expression::Expr
        "Constructor and objective folding expressions, used in generating the invocation call."
        expressions::Vector{Expr}
        "The maximum number of objective caches this model will need."
        maximum_objective_cache_count::Int
        "How many objective caches are currently active."
        objective_cache_count::Int
    end

The composite equivalent of [`ModelInfo`](@ref), augmented to track the model
symbol (`a1`, `m3`, etc.), and the model parameters (`K_1`, `a_3`, etc.)
"""
struct CompositeModelInfo
    "The parameter symbols of the model with the respective lens to the actual parameter."
    parameter_symbols::Vector{Pair{Symbol,Lens}}
    "The parameter symbols of the model with the respective lens to the actual parameter."
    model_symbols::Dict{Symbol,Vector{Symbol}}
    "Each model assigned to a unique symbol."
    models::Vector{Pair{Symbol,ModelInfo}}
    "The expression representing the folding operations of this composite model."
    model_expression::Expr
    "Constructor and objective folding expressions, used in generating the invocation call."
    expressions::Vector{Expr}
    "The maximum number of objective caches this model will need."
    maximum_objective_cache_count::Int
end

function CompositeModelInfo(
    info_map::Vector{Pair{Symbol,ModelInfo}},
    model_expr::Expr,
    expressions::Vector{Expr},
    max_obj_count::Int,
)
    parameter_map = Pair{Symbol,Lens}[]
    model_symbol_map = Dict{Symbol,Vector{Symbol}}()

    for (ms, info) in info_map
        unique_symbs = Symbol[]
        for (s, l) in zip(info.symbols, parameter_lenses(info))
            unique_symbol = make_unique_symbol(parameter_map, s; delim = "_")
            push!(parameter_map, unique_symbol => l)
            push!(unique_symbs, unique_symbol)
        end
        model_symbol_map[ms] = unique_symbs
    end

    CompositeModelInfo(
        parameter_map,
        model_symbol_map,
        info_map,
        model_expr,
        expressions,
        max_obj_count,
    )
end

function make_unique_symbol(symbol_pairs, new; delim = "")
    i = 1
    symb = Symbol(new, delim, i)
    while any(i -> i[1] == symb, symbol_pairs)
        i += 1
        symb = Symbol(new, delim, i)
    end
    symb
end

mutable struct CompositeAggregation
    info_map::Vector{Pair{Symbol,ModelInfo}}
    expressions::Vector{Expr}
    max_obj_count::Int
    current_obj::Int
    CompositeAggregation() = new(Pair{Symbol,ModelInfo}[], Expr[], 0, 0)
end

"""
    set_objective_count!(a::CompositeAggregation, o::Int)

Set the objective count to a specific value. Will bump the maximum objective
count if the new value exceeds the current maximum.
"""
function set_objective_count!(a::CompositeAggregation, o::Int)
    @assert o >= 0
    a.current_obj = o
    if a.current_obj > a.max_obj_count
        a.max_obj_count = a.current_obj
    end
    a.current_obj
end
inc_ob_count!(a::CompositeAggregation) = set_objective_count!(a, a.current_obj + 1)
dec_ob_count!(a::CompositeAggregation) = set_objective_count!(a, a.current_obj - 1)
make_objective_symbol(a::CompositeAggregation) = Symbol(OBJECTIVE_SYMBOL, a.current_obj)

"""
    get_info(model::Type{<:AbstractSpectralModel}, lens::Lens)

Returns a [`ModelInfo`](@ref) struct for a given model.
"""
function get_info(@nospecialize(model::Type{<:AbstractSpectralModel}), lens::Lens)
    symbs = [get_parameter_symbols(model)...]
    gen_symbs = map(Base.gensym, symbs)

    closure_symbs = [get_closure_symbols(model)...]
    gen_closure_symbs = map(Base.gensym, closure_symbs)

    ModelInfo(symbs, gen_symbs, closure_symbs, gen_closure_symbs, lens, model)
end
function get_info(
    @nospecialize(model::Type{<:CompositeModel}),
    lens::Lens;
    T::Type = Float64,
)
    agg = CompositeAggregation()
    expr = _add_composite_info!(agg, model, lens, T)
    CompositeModelInfo(agg.info_map, expr, agg.expressions, agg.max_obj_count)
end

"""
Used exclusively to do recursive [`CompositeModel`](@ref) parsing.
"""
function _add_composite_info!(
    agg::CompositeAggregation,
    model::Type{<:AbstractSpectralModel},
    lens::Lens,
    T::Type,
)
    model_symbol = MODEL_SYMBOL_LOOKUP[modelkind(model)]
    unique_model_symbol = make_unique_symbol(agg.info_map, model_symbol)

    info = get_info(model, lens)
    push!(agg.info_map, unique_model_symbol => info)

    if modelkind(model) !== Convolutional()
        inc_ob_count!(agg)
    end

    obj = make_objective_symbol(agg)
    NumType = T <: FitParam ? T.parameters[1] : T

    invoke_call = make_invoke_call(info, obj, DOMAIN_SYMBOL, NumType)
    push!(agg.expressions, invoke_call)

    unique_model_symbol
end
function _add_composite_info!(
    agg::CompositeAggregation,
    model::Type{<:CompositeModel},
    lens::Lens,
    T::Type,
)
    right_run = :(getfield($lens, :right))
    left_run = :(getfield($lens, :left))
    # order is very important
    r = _add_composite_info!(agg, model.parameters[2], right_run, T)
    l = _add_composite_info!(agg, model.parameters[1], left_run, T)
    op = operation_symbol(model.parameters[3])
    if isnothing(op)
        :(($l)($r))
    else
        add_objective_reduction!(agg, op)
        Expr(:call, op, l, r)
    end
end

"""
    add_objective_reduction!(ra::CompositeAggregation, op::Symbol)

Reduces the current objective count and applies the reduction operation `op` to
them. For example, if `op` is `:+`, and the objective count is 3, then after
this function has been called the objective count will be 2 and the reduction
expression

```julia
@. flux2 = flux2 + flux3
```

will have been added to the [`CompositeAggregation`](@ref).
"""
function add_objective_reduction!(a::CompositeAggregation, op::Symbol)
    r = make_objective_symbol(a)
    dec_ob_count!(a)
    l = make_objective_symbol(a)
    expr = Expr(:call, op, l, r)
    push!(a.expressions, :(@.($l = $expr)))
end

"""
    assemble_objective_unpack(N)

Assembles the statements for unpacking the objective cache into a number of
views. Assembles the part of the model call that looks like:
```julia
objective1 = view(objectives, :, 1), objective2 = ...
```
for `N` objectives slices.
"""
function assemble_objective_unpack(N)
    symbols = (Symbol(OBJECTIVE_SYMBOL, i) for i = 1:N)
    unpacks = [:($s = view(objectives, :, $i)) for (i, s) in enumerate(symbols)]
    unpacks
end

"""
    assemble_composite_model_call(model::Type{<:CompositeModel}, parameters::Type{<:AbstractVector})

Assemble the full composite model call, with objective unpacking via
[`assemble_objective_unpack`](@ref), closure and parameter assignments, model
invocation, and objective reduction. Uses [`assemble_fast_call`](@ref) to put
the final function body together.
"""
function assemble_composite_model_call(
    @nospecialize(model::Type{<:CompositeModel}),
    parameters,
)
    # propagate information about free parameters to allow for AD
    info = get_info(model, DEFAULT_LENS; T = eltype(parameters))
    unpack = assemble_objective_unpack(info.maximum_objective_cache_count)
    i = 0
    parameter_assigns = Expr[]
    closure_assigns = Expr[]
    for (_, info) in info.models
        for s in info.generated_symbols
            assignment = :($(s) = parameters[$(i += 1)])
            push!(parameter_assigns, assignment)
        end

        for (c_lens, c) in
            zip(closure_parameter_lenses(info), info.generated_closure_symbols)
            c_assingment = :($c = $c_lens)
            push!(closure_assigns, c_assingment)
        end
    end
    assemble_fast_call(unpack, closure_assigns, parameter_assigns, info.expressions)
end

function assemble_fast_call(
    unpack::Vector{Expr},
    closure::Vector{Expr},
    params::Vector{Expr},
    body::Vector{Expr},
)
    quote
        @fastmath begin
            @inbounds let ($(unpack...))
                $(closure...)
                $(params...)
                $(body...)
                return $(Symbol(OBJECTIVE_SYMBOL, 1))
            end
        end
    end
end

function assemble_parameter_named_tuple(model::Type{<:CompositeModel})
    info = get_info(model, DEFAULT_LENS; T = Float64)
    keys = (first.(info.parameter_symbols)...,)
    values = (last.(info.parameter_symbols)...,)
    :(NamedTuple{$(keys)}(($(values...),)))
end
function assemble_parameter_named_tuple(model::Type{<:AbstractSpectralModel})
    info = get_info(model, DEFAULT_LENS)
    keys = (info.symbols...,)
    values = (parameter_lenses(info)...,)
    :(NamedTuple{$(keys)}(($(values...),)))
end

# API for models

function get_parameter_symbols(model::Type{<:AbstractSpectralModel})
    fieldnames(model)
end
function get_parameter_symbols(model::Type{<:CompositeModel})
    info = get_info(model, DEFAULT_LENS)
    first.(info.parameter_symbols)
end

function get_closure_symbols(::Type{<:AbstractSpectralModel})
    ()
end

"""
    make_constructor(model::Type{<:AbstractSpectralModel}, closures::Vector, params::Vector, T::Type)

Create a constructor expression for the model.  Should return something similar
to
```julia
:(PowerLaw{T}(arg1, arg2, arg3)))
```
unpacking the `closures` and `params` vectors in the appropriate places.
"""
function make_constructor(
    model::Type{<:AbstractSpectralModel},
    closures::Vector,
    params::Vector,
    T::Type,
)
    :($(model.name.wrapper){$(model.parameters[1:(end-1)]...),$(T)}(
        $(closures...),
        $(params...),
    ))
end

"""
    parameter_lenses(::Type{<:AbstractSpectralModel}, info::ModelInfo)

Return a vector of lenses [`Lens`](@ref) that refer to each of the models
parameters. The lenses should be relative to `model.lens`.
"""
function parameter_lenses(::Type{<:AbstractSpectralModel}, info::ModelInfo)
    map(info.symbols) do symb
        :(getfield($(info.lens), $(Meta.quot(symb))))
    end
end
parameter_lenses(info::ModelInfo) = parameter_lenses(info.model, info)

function closure_parameter_lenses(::Type{<:AbstractSpectralModel}, info::ModelInfo)
    map(info.closure_symbols) do symb
        :(getfield($(info.lens), $(Meta.quot(symb))))
    end
end
closure_parameter_lenses(info::ModelInfo) = closure_parameter_lenses(info.model, info)

end # module Reflection

# public API wrappers

struct DestructuredModel{T}
    "Maps the model symbol to the specific model."
    model_map::Vector{Pair{Symbol,AbstractSpectralModel}}
    "Maps model symbols to the parameter symbols."
    parameter_symbols::Dict{Symbol,Vector{Symbol}}
    "Maps each parameter symbol to its parameter."
    parameter_map::Dict{Symbol,T}
    "Model expression"
    expression::Expr
end

@inline function destructure_model(model::AbstractSpectralModel)
    error("Can only destructure composite models.")
end

@inline @generated function destructure_model(model::CompositeModel)
    info = Reflection.get_info(model, :model; T = Float64)
    lenses = [:($(Meta.quot(s)) => $(i.lens)) for (s, i) in info.models]
    params = [:($(Meta.quot(s)) => $(lens)) for (s, lens) in info.parameter_symbols]
    quote
        DestructuredModel(
            Pair{Symbol,AbstractSpectralModel}[$(lenses...)],
            Dict{Symbol,Vector{Symbol}}($(info.model_symbols...)),
            Dict($(params...)),
            $(Meta.quot(info.model_expression)),
        )
    end
end

@inline @generated function parameter_named_tuple(model::AbstractSpectralModel)
    Reflection.assemble_parameter_named_tuple(model)
end

@inline function parameter_tuple(model::AbstractSpectralModel)
    nt = parameter_named_tuple(model)
    Tuple(nt)
end

@inline @generated function remake_model_with_parameters(
    model::AbstractSpectralModel,
    parameters::Base.AbstractVecOrTuple,
)
    Reflection.reassemble_model(model, parameters)
end

@inline @generated function objective_cache_count(model::CompositeModel)
    info = Reflection.get_info(model, :model; T = Float64)
    :($(info.maximum_objective_cache_count))
end
@inline function objective_cache_count(model::AbstractSpectralModel)
    1
end

@inline @generated function composite_model_call!(
    objectives,
    domain,
    model::CompositeModel,
    parameters,
)
    Reflection.assemble_composite_model_call(model, parameters)
end
