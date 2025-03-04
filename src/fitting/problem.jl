struct FittableMultiDataset{D}
    d::D
    FittableMultiDataset(data::Vararg{<:AbstractDataset}) = new{typeof(data)}(data)
end

function Base.getindex(multidata::FittableMultiDataset, i)
    multidata.d[i]
end

function Base.merge(m1::FittableMultiDataset, m2::FittableMultiDataset)
    FittableMultiDataset(m1.d..., m2.d...)
end

struct FittableMultiModel{M}
    m::M
    FittableMultiModel(model::Vararg{<:AbstractSpectralModel}) = new{typeof(model)}(model)
end

function translate_bindings(
    model_index::Int,
    m::FittableMultiModel,
    bindings::Vector{Vector{Pair{Int,Int}}},
)
    # map the parameter index to a string to display
    translation = Dict{Int,Pair{paramtype(m.m[1]),String}}()
    for b in bindings
        # we skip the first item in the list since that is the root of the binding
        root = b[1]

        params = parameter_named_tuple(m.m[first(root)])
        symbol = propertynames(params)[last(root)]
        value = params[last(root)]

        for pair in @views b[2:end]
            # check if this binding applies to the current model
            if first(pair) == model_index
                translation[last(pair)] = value => "~Model $(first(root)) $(symbol)"
            end
        end
    end

    translation
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(model::FittableMultiModel))
    buff = IOBuffer()
    buff_c = IOContext(buff, io)
    println(buff_c, "Models:")
    for (i, m) in enumerate(model.m)
        buf = IOBuffer()

        println(buf)
        printstyled(buf, "Model $i", color = :yellow)
        print(buf, ": ")

        _printinfo(buf, m)
        print(buff_c, indent(String(take!(buf)), 2))
    end
    print(io, encapsulate(String(take!(buff))))
end

function Base.getindex(multimodel::FittableMultiModel, i)
    multimodel.m[i]
end

function _multi_constructor_wrapper(
    T::Union{<:Type{<:FittableMultiDataset},<:Type{<:FittableMultiModel}},
    args,
)
    if args isa T
        args
    else
        if args isa Tuple
            @show eltype(args)
            if eltype(args) <: T
                reduce(merge, args)
            else
                T(args...)
            end
        else
            T(args)
        end
    end
end

"""
    _accumulated_indices(items)

`items` is a tuple or vector of lengths `n1, n2, ...`

Returns a tuple or array with same length as items, which gives the index boundaries of
an array with size `n1 + n2 + ...`.
"""
function _accumulated_indices(items)
    total::Int = 0
    map(items) do item
        total = item + total
        total
    end
end

struct FittingProblem{M<:FittableMultiModel,D<:FittableMultiDataset,B}
    model::M
    data::D
    bindings::B
    function FittingProblem(m::FittableMultiModel, d::FittableMultiDataset)
        bindings = Vector{Pair{Int,Int}}[]
        new{typeof(m),typeof(d),typeof(bindings)}(m, d, bindings)
    end
end

FittingProblem(pairs::Vararg{<:Pair}) = FittingProblem(first.(pairs), last.(pairs))

function FittingProblem(m, d)
    _m = _multi_constructor_wrapper(FittableMultiModel, m)
    _d = _multi_constructor_wrapper(FittableMultiDataset, d)
    FittingProblem(_m, _d)
end

function model_count(prob::FittingProblem)
    return length(prob.model.m)
end

function data_count(prob::FittingProblem)
    return length(prob.data.d)
end

function Base.show(io::IO, ::MIME"text/plain", @nospecialize(prob::FittingProblem))
    buff = IOBuffer()
    buff_c = IOContext(buff, io)

    println(buff_c, "FittingProblem:")
    println(buff_c, "  . Models     : $(length(prob.model.m))")
    println(buff_c, "  . Datasets   : $(length(prob.data.d))")
    println(buff_c, "  Parameter Summary:")

    factor = div(length(prob.data.d), length(prob.model.m))

    total = 0
    free = 0
    frozen = 0
    for m in prob.model.m
        for p in parameter_tuple(m)
            total += factor
            if isfrozen(p)
                frozen += factor
            else
                free += factor
            end
        end
    end

    bound = 0
    for binding in prob.bindings
        bound += length(Set(first.(binding))) - 1
    end
    free = free - bound

    println(buff_c, "  . Total      : $(total)")

    print(buff_c, "  . ")
    printstyled(buff_c, "Frozen", color = :cyan)
    println(buff_c, "     : $(frozen)")

    print(buff_c, "  . ")
    printstyled(buff_c, "Bound", color = :cyan)
    println(buff_c, "      : $(bound)")

    print(buff_c, "  . ")
    printstyled(buff_c, "Free", color = :cyan)
    println(buff_c, "       : $(free)")

    print(io, encapsulate(String(take!(buff))))
end


"""
    details(prob::FittingProblem)

Show details about the fitting problem, including the specific model parameters that are bound together.
"""
function details(prob::FittingProblem; color = true)
    buff = IOBuffer()
    buff_ctx = IOContext(buff, :color => color)

    # TODO: this is a horribly implemented function
    print(buff_ctx, "Models:")
    for (i, m) in enumerate(prob.model.m)
        buf = IOBuffer()
        buf_ctx = IOContext(buf, buff_ctx)

        println(buf_ctx)
        printstyled(buf_ctx, "Model $i", color = :yellow)
        print(buf_ctx, ": ")

        _printinfo(buf_ctx, m; bindings = translate_bindings(i, prob.model, prob.bindings))
        print(buff_ctx, indent(String(take!(buf)), 2))
    end

    print(encapsulate(String(take!(buff))))
end

export FittingProblem, FittableMultiModel, FittableMultiDataset, details
