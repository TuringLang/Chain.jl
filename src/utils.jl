@generated function initnamemap(namemap::NamedTuple{names}) where names
    if :parameters in names
        :((; $((:($name = Symbol.(namemap.$name)) for name in names)...)))
    else
        :((; parameters = Symbol[],
          $((:($name = Symbol.(namemap.$name)) for name in names)...)))
    end
end

# fallback
function initnamemap(namemap)
    (; parameters = Symbol[], (Symbol(key) => Symbol.(values) for (key, values) in namemap)...)
end

"""
    namemap_intersect(namemap::NamedTuple, names)

This is an internal function used to remove values from a name map
and return a new name map.
"""
@generated function namemap_intersect(namemap::NamedTuple{sections}, names) where {sections}
    :((; $((:($section = intersect(namemap.$section, names)) for section in sections)...)))
end

"""
    string2symbol(x)

Convert strings to symbols.

If `x isa String`, the corresponding `Symbol` is returned. Likewise, if
`x isa AbstractVector{String}`, the corresponding vector of `Symbol`s is returned. In all
other cases, input `x` is returned.
"""
string2symbol(x) = x
string2symbol(x::String) = Symbol(x)
string2symbol(x::AbstractVector{String}) = Symbol.(x)

#################### Mathematical Operators ####################
function cummean(x::AbstractArray)
    return mapslices(cummean, x, dims = 1)
end

function cummean(x::AbstractVector)
    y = similar(x, Float64)
    xs = 0.0
    fill!(y, xs)

    c = 0
    for i in 1:length(x)
        if !ismissing(x[i])
            c += 1
            xs += x[i]
        end
        y[i] = xs / c
    end
    return y
end

## Csorgo S and Faraway JJ. The exact and asymptotic distributions of the
## Cramer-von Mises statistic. Journal of the Royal Statistical Society,
## Series B, 58: 221-234, 1996.
function pcramer(q::Real)
    p = 0.0
    for k in 0:3
        c1 = 4.0 * k + 1.0
        c2 = c1^2 / (16.0 * q)
        p += gamma(k + 0.5) / factorial(k) * sqrt(c1) * exp(-c2) * besselk(0.25, c2)
    end
    return p / (pi^1.5 * sqrt(q))
end

function _dict2namedtuple(d::Dict)
    t_keys = ntuple(x -> Symbol(collect(keys(d))[x]), length(keys(d)))
    t_vals = ntuple(x -> collect(values(d))[x], length(values(d)))
    return NamedTuple{t_keys}(t_vals)
end

function _namedtuple2dict(n::NamedTuple)
    return Dict(pairs(n))
end

# adapted from https://github.com/JuliaLang/julia/blob/3fdfb6734bb6ed7fc805a484183077dd7924b7c0/base/namedtuple.jl#L222
function merge_union(a::NamedTuple{an}, b::NamedTuple{bn}) where {an,bn}
    if @generated
        names = Base.merge_names(an, bn)
        types = merge_union_types(names, a, b)

        values = map(names) do n
            if Base.sym_in(n, an)
                if Base.sym_in(n, bn)
                    :(union(getfield(a, $(QuoteNode(n))), getfield(b, $(QuoteNode(n)))))
                else
                    :(getfield(a, $(QuoteNode(n))))
                end
            else
                :(getfield(b, $(QuoteNode(n))))
            end
        end
        
        return :(NamedTuple{$names,$types}(($(values...),)))
    else
        names = Base.merge_names(an, bn)
        types = merge_union_types(names, typeof(a), typeof(b))

        values = map(names) do n
            if Base.sym_in(n, an)
                if Base.sym_in(n, bn)
                    union(getfield(a, n), getfield(b, n))
                else
                    getfield(a, n)
                end
            else
                getfield(b, n)
            end
        end
        
        return NamedTuple{names,types}(values)
    end
end

# adapted from https://github.com/JuliaLang/julia/blob/3fdfb6734bb6ed7fc805a484183077dd7924b7c0/base/namedtuple.jl#L192
Base.@pure function merge_union_types(names::Tuple{Vararg{Symbol}}, a::Type{<:NamedTuple}, b::Type{<:NamedTuple})
    an = Base._nt_names(a)
    bn = Base._nt_names(b)

    types = map(names) do n
        if Base.sym_in(n, an)
            if Base.sym_in(n, bn)
                promote_type(fieldtype(a, n), fieldtype(b, n))
            else
                fieldtype(a, n)
            end
        else
            fieldtype(b, n)
        end
    end

    return Tuple{types...}
end

# promote element types of the tail of a NamedTuple
function promote_eltype_namedtuple_tail(::NamedTuple{k,v}) where {k,v}
    return promote_eltype_tuple_type(Base.tuple_type_tail(v))
end

# promote element types of a tuple
promote_eltype_tuple_type(::Type{Tuple{}}) = Any
promote_eltype_tuple_type(::Type{Tuple{T}}) where T = T
function promote_eltype_tuple_type(t::Type{<:Tuple})
    Base.promote_eltype(Base.tuple_type_head(t), promote_eltype_tuple_type(Base.tuple_type_tail(t)))
end

"""
    cskip(x)

Wrapper for `collect(skipmissing(x))`.
"""
cskip(x) = collect(skipmissing(x))

# check if `T` and recursive `eltype` of `T` are concrete types
function isconcretetype_recursive(T)
    return isconcretetype(T) && (eltype(T) === T || isconcretetype_recursive(eltype(T)))
end

"""
    concretize(x)

Return a container whose elements are as concrete as possible by, e.g., replacing
element types of `Real` or `Union{Missing,Float64}` with `Float64`, if appropriate.

Note that this function is only type stable if the containers are already concretely
typed. In this case the input `x` is returned, and hence input and output share the
same underlying data. In general, however, the input and the output do not share the same
data.
"""
concretize(x) = x
function concretize(x::AbstractArray)
    if isconcretetype_recursive(typeof(x))
        return x
    else
        xnew = map(concretize, x)
        T = mapreduce(typeof, promote_type, xnew)
        if T <: eltype(xnew)
            return convert(AbstractArray{T}, xnew)
        else
            return xnew
        end
    end
end

function concretize(x::Chains)
    value = x.value
    if isconcretetype_recursive(value)
        return x
    else
        return Chains(concretize(value), x.logevidence, x.name_map, x.info)
    end
end

function isstrictlyincreasing(x::AbstractVector{Int})
    return isempty(x) || _isstrictlyincreasing_nonempty(x)
end

_isstrictlyincreasing_nonempty(x::AbstractRange{Int}) = step(x) > 0
function _isstrictlyincreasing_nonempty(x::AbstractVector{Int})
    i = first(x)
    for j in Iterators.drop(x, 1)
        j > i || return false
        i = j
    end
    return true
end
