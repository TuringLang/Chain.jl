#################### Chains ####################

## Constructors ##

# Default name map.
const DEFAULT_MAP = Dict{Symbol, Vector{Symbol}}(:parameters => Symbol[])

# Constructor to handle a vector of vectors.
Chains(val::AbstractVector{<:AbstractVector{<:Union{Missing, Real}}}, args...; kwargs...) =
	Chains(copy(reduce(hcat, val)'), args...; kwargs...)

# Constructor to handle a 1D array.
Chains(val::AbstractVector{<:Union{Missing, Real}}, args...; kwargs...) =
	Chains(reshape(val, :, 1, 1), args...; kwargs...)

# Constructor to handle a 2D array
Chains(val::AbstractMatrix{<:Union{Missing, Real}}, args...; kwargs...) =
    Chains(reshape(val, size(val, 1), size(val, 2), 1), args...; kwargs...)

# Constructor to handle parameter names that are not Symbols.
function Chains(
    val::AbstractArray{<:Union{Missing,Real},3},
    parameter_names::AbstractVector,
    args...;
    kwargs...
)
    return Chains(val, Symbol.(parameter_names), args...; kwargs...)
end

# Generic chain constructor.
function Chains(
    val::AbstractArray{<:Union{Missing, Real},3},
    parameter_names::AbstractVector{Symbol} = Symbol.(:param, 1:size(val, 2)),
    name_map = (parameters = parameter_names,);
    start::Int = 1,
    thin::Int = 1,
    evidence = missing,
    info::NamedTuple = NamedTuple(),
    sorted::Bool = false
)
    # Make sure that we have a `:parameters` index and # Copying can avoid state mutation.
    _name_map = initnamemap(name_map)

    # Preclean the name_map of names that aren't in the
    # parameter_names vector.
    for names in _name_map
        filter!(x -> x ∈ parameter_names, names)
    end

    # Store unassigned variables.
    unassigned = Set(Symbol[])

    # Check that all parameters are assigned.
    for param in parameter_names
        if all(param ∉ names for names in _name_map)
            push!(unassigned, param)
        end
    end

    # Assign all unassigned parameter names.
    append!(_name_map[:parameters], unassigned)

    # Construct the AxisArray.
    arr = AxisArray(val;
                    iter = range(start, step=thin, length=size(val, 1)),
                    var = parameter_names,
                    chain = 1:size(val, 3))

    # Create the new chain.
    chains = Chains(arr, evidence, _name_map, info)

    if sorted
        return sort(chains)
    else
        return chains
    end
end

# Retrieve a new chain with only a specific section pulled out.
Chains(::Chains, ::Tuple{}; kwargs...) = throw(ArgumentError("no section specified"))
Chains(c::Chains, section; kwargs...) = Chains(c, (section,); kwargs...)
function Chains(c::Chains, section::AbstractArray; kwargs...)
	Chains(c, ntuple(i -> section[i], length(section)); kwargs...)
end
function Chains(c::Chains, section::NTuple{<:Any,Symbol}; sorted::Bool=false)
    # Make sure the section exists first.
	section ⊆ keys(c.name_map) ||
		throw(ArgumentError("$section is not a subset of the chain's name map"))

	# Create the new section map.
    name_map = NamedTuple{section}(getindex.(Ref(c.name_map), section))

    # Extract wanted values.
    value = c.value[:, mapreduce(collect, vcat, name_map), :]

    # Create the new chain.
    chain = Chains(value, c.logevidence, name_map, c.info)

    if sorted
        return sort(chain)
    else
        return chain
    end
end

#################### Indexing ####################

function _sym2index(chains::Chains, syms::AbstractVector{Symbol}; sorted::Bool = false)
    # Start by looking up the symbols in the list of parameter names.
    names_of_params = names(chains)
    indices = indexin(syms, names_of_params)

    allsyms = Symbol[]
    sizehint!(allsyms, length(syms))

    for (index, sym) in zip(indices, syms)
        if index === nothing
            prefix = string(sym) * "["
            for name in names_of_params
                startswith(string(name), prefix) && push!(allsyms, name)
            end
        else
            push!(allsyms, sym)
        end
    end

    return sorted ? sort!(allsyms, lt=natural) : allsyms
end

Base.getindex(c::Chains, i1::T) where T<:Union{AbstractUnitRange, StepRange} = c[i1, :, :]
Base.getindex(c::Chains, i1::Integer) = c[i1:i1, :, :]

Base.getindex(c::Chains, v::String) = c[:, Symbol(v), :]
Base.getindex(c::Chains, v::AbstractVector{String}) = c[:, Symbol.(v), :]

Base.getindex(c::Chains, v::Symbol) = c[[v]]
function Base.getindex(c::Chains, v::AbstractVector{Symbol})
    syms = _sym2index(c, v)
    return c[:, syms, :]
end

function Base.getindex(c::Chains, i...)
    # Make sure things are in array form to preserve the axes.
    ind = [typeof(i[1]) <: Integer ? (i[1]:i[1]) : i[1],
           typeof(i[2]) <: Union{AbstractArray, Colon} ?  i[2] : [i[2]],
           typeof(i[3]) <: Union{AbstractArray, Colon} ?  i[3] : [i[3]]
    ]

    # Check to see if we received a symbol or a string in i[2].
    if !isa(ind[2], Colon)
        if eltype(ind[2]) <: Symbol
            ind[2] = _sym2index(c, ind[2])
        else
            ind[2] = _sym2index(c, Symbol.(ind[2]))
        end
    end

    newval = getindex(c.value, ind...)
    names = newval.axes[2].val
    new_name_map = _trim_name_map(names, c.name_map)
    return Chains(newval, c.logevidence, new_name_map, c.info)
end

Base.setindex!(c::Chains, v, i...) = setindex!(c.value, v, i...)
Base.lastindex(c::Chains) = lastindex(c.value, 1)
Base.lastindex(c::Chains, d::Integer) = lastindex(c.value, d)

"""
    Base.get(c::Chains, v::Symbol; flatten=false)
    Base.get(c::Chains, vs::Vector{Symbol}; flatten=false)

Returns a `NamedTuple` with `v` as the key, and matching paramter
names as the values.

Passing `flatten=true` will return a `NamedTuple` with keys ungrouped.

Example:

```julia
x = get(c, :param1)
x = get(c, [:param1, :param2])
```
"""
Base.get(c::Chains, v::Symbol; flatten = false) = get(c, [v], flatten=flatten)
function Base.get(c::Chains, vs::Vector{Symbol}; flatten = false)
    pairs = Dict()
    for v in vs
        syms = _sym2index(c, [v])
        len = length(syms)
        val = ()
        if len > 1
            val = ntuple(i -> c.value[:,syms[i],:], length(syms))
        elseif len == 1
            val = c.value[:,syms[1],:]
        else
            continue
        end

        if flatten
            for i in eachindex(syms)
                pairs[syms[i]] = val[i]
            end
        else
            pairs[v] = val
        end
    end
    return _dict2namedtuple(pairs)
end

"""
    get(c::Chains; section::Union{Vector{Symbol}, Symbol; flatten=false}

Returns all parameters in a given section(s) as a `NamedTuple`.

Passing `flatten=true` will return a `NamedTuple` with keys ungrouped.

Example:

```julia
x = get(chn, section = :parameters)
x = get(chn, section = [:internals, :parameters])
```
"""
function Base.get(
    c::Chains;
    section::Union{Symbol,AbstractVector{Symbol}},
    flatten = false
)
    names = Set(Symbol[])
    regex = r"[^\[]*"
    _section = section isa Symbol ? (section,) : section
    for v in _section
        v in keys(c.name_map) || error("section $v does not exist")

        # If the name contains a bracket,
        # split it so get can group them correctly.
        if flatten
            append!(names, c.name_map[v])
        else
            for name in c.name_map[v]
                m = match(regex, string(name))
                push!(names, Symbol(m.match))
            end
        end
    end

    return get(c, collect(names); flatten = flatten)
end

"""
    get_params(c::Chains; flatten=false)

Returns all parameters packaged as a `NamedTuple`. Variables with a bracket
in their name (as in "P[1]") will be grouped into the returned value as P.

Passing `flatten=true` will return a `NamedTuple` with keys ungrouped.

Example:

```julia
x = get_params(chn)
x.P
```
"""
get_params(c::Chains; flatten = false) = get(c, section = sections(c), flatten=flatten)

#################### Base Methods ####################

function Base.show(io::IO, c::Chains)
    print(io, "Object of type Chains, with data of type $(summary(c.value.data))\n\n")
    println(io, header(c))

    # Show summary stats.
    show(io, describe(c))
end

Base.keys(c::Chains) = names(c)
Base.size(c::Chains) = size(c.value)
Base.size(c::Chains, ind) = size(c)[ind]
Base.length(c::Chains) = length(range(c))
Base.first(c::Chains) = first(c.value[Axis{:iter}].val)
Base.step(c::Chains) = step(c.value[Axis{:iter}].val)
Base.last(c::Chains) = last(c.value[Axis{:iter}].val)

Base.convert(::Type{Array}, chn::Chains) = convert(Array, chn.value)

#################### Auxilliary Functions ####################

"""
    range(c::Chains)

Returns the range used in a `Chains` object.
"""
function Base.range(c::Chains)
    return c.value[Axis{:iter}].val
end

"""
    chains(c::Chains)

Returns the names or symbols of each chain in a `Chains` object.
"""
function chains(c::Chains)
    return c.value[Axis{:chain}].val
end

"""
    names(chains::Chains)

Return the parameter names in the `chains`.
"""
Base.names(chains::Chains) = chains.value[Axis{:var}].val

"""
    names(chains::Chains, section::Symbol)

Return the parameter names of a `section` in the `chains`.
"""
Base.names(chains::Chains, section::Symbol) = convert(Vector{Symbol}, chains.name_map[section])

"""
    names(chains::Chains, sections)

Return the parameter names of the `sections` in the `chains`.
"""
function Base.names(c::Chains, sections)
    names = Symbol[]
    for section in sections
        append!(names, c.name_map[section])
    end
    return names
end

"""
    get_sections(c::Chains, sections::Vector = [])

Returns multiple `Chains` objects, each containing only a single section.
"""
function get_sections(c::Chains, sections::Vector = [])
    sections = length(sections) == 0 ? collect(keys(c.name_map)) : sections
    return [Chains(c, section) for section in sections]
end

# Return a new chain for each section.
function get_sections(c::Chains, section::Union{Symbol, String})
    return get_sections(c, [section])
end

"""
    sections(c::Chains)

Retrieve a list of the sections in a chain.
"""
sections(c::Chains) = collect(keys(c.name_map))

"""
    header(c::Chains; section=missing)

Returns a string containing summary information for a `Chains` object.
If the `section` keyword is used, this function prints only the relevant section
header.

Example:

```julia
# Printing the whole header.
header(chn)

# Print only one section's header.
header(chn, section = :parameter)
```
"""
function header(c::Chains; section=missing)
    rng = range(c)

    # Function to make section strings.
    section_str(sec, arr) = string(
        "$sec",
        repeat(" ", 18 - length(string(sec))),
        "= $(join(map(string, arr), ", "))\n"
    )

    # Set up string array.
    section_strings = String[]

    # Get section lines.
    if section isa Missing
        for (sec, nms) in pairs(c.name_map)
            section_string = section_str(sec, nms)
            push!(section_strings, section_string)
        end
    else
        section in keys(c.name_map) ||
            throw(ArgumentError("$section not found in name map."))
        section_string = section_str(section, c.name_map[section])
        push!(section_strings, section_string)
    end

    #

    # Return header.
    return string(
        ismissing(c.logevidence) ? "" : "Log evidence      = $(c.logevidence)\n",
        "Iterations        = $(first(c)):$(last(c))\n",
        "Thinning interval = $(step(c))\n",
        "Chains            = $(join(map(string, chains(c)), ", "))\n",
        "Samples per chain = $(length(range(c)))\n",
        section_strings...
    )
end

function indiscretesupport(c::Chains,
                           bounds::Tuple{Real, Real}=(0, Inf))
  nrows, nvars, nchains = size(c.value)
  result = Array{Bool}(undef, nvars * (nrows > 0))
  for i in 1:nvars
    result[i] = true
    for j in 1:nrows, k in 1:nchains
      x = c.value[j, i, k]
      if !isinteger(x) || x < bounds[1] || x > bounds[2]
        result[i] = false
        break
      end
    end
  end
  result
end

function link(c::Chains)
  cc = copy(c.value.data)
  for j in axes(cc, 2)
    x = cc[:, j, :]
    if minimum(x) > 0.0
      cc[:, j, :] = maximum(x) < 1.0 ? logit.(x) : log.(x)
    end
  end
  cc
end

"""
    _trim_name_map(names::Vector, name_map::NamedTuple)

This is an internal function used to remove values from a name map
and return a new name_map.
"""
function _trim_name_map(names::Vector, name_map::NamedTuple)
    n = _namedtuple2dict(name_map)
    for (key, values) in n
        intersection = values ∩ names
        if length(intersection) > 0
            n[key] = intersection
        else
            delete!(n, key)
        end
    end
    return _dict2namedtuple(n)
end

### Chains specific functions ###
"""
    sort(c::Chains)

Returns a new column-sorted version of `c`, using natural sort order.
"""
function Base.sort(c::Chains)
    v = c.value
    x, y, z = size(v)
    unsorted = collect(zip(1:y, v.axes[2].val))
    sorted = sort(unsorted, by = x -> string(x[2]), lt=natural)
    new_axes = (v.axes[1], Axis{:var}([n for (_, n) in sorted]), v.axes[3])
    new_v = copy(v.data)
    for i in eachindex(sorted)
        new_v[:, i, :] = v[:, sorted[i][1], :]
    end

    aa = AxisArray(new_v, new_axes...)

    # Sort the name map too:
    namemap = deepcopy(c.name_map)
    for names in namemap
        sort!(names, by=string, lt=natural)
    end

    return Chains(aa, c.logevidence, namemap, c.info)
end

"""
    setinfo(c::Chains, n::NamedTuple)

Returns a new `Chains` object with a `NamedTuple` type `n` placed in the `info` field.

Example:

```julia
new_chn = setinfo(chn, NamedTuple{(:a, :b)}((1, 2)))
```
"""
function setinfo(c::Chains, n::NamedTuple)
    return Chains(c.value, c.logevidence, c.name_map, n)
end

"""
    set_section(chains::Chains, namemap)

Create a new `Chains` object from `chains` with the provided `namemap` mapping of parameter
names.

Both chains share the same underlying data. Any parameters in the chain that are unassigned
will be placed into the `:parameters` section.
"""
function set_section(chains::Chains, namemap)
    # Initialize the name map.
    _namemap = initnamemap(namemap)

    # Make sure all the names are in the new name map.
    newnames = Set(Symbol[])
    names_of_params = names(chains)
    for names in _namemap
        filter!(x -> x ∈ names_of_params, names)
        for name in names
            push!(newnames, name)
        end
    end
    missingnames = setdiff(names_of_params, newnames)

    # Assign everything that is missing to :parameters.
    if !isempty(missingnames)
        @warn "Section mapping does not contain all parameter names, " *
            "$missingnames assigned to :parameters."
        for name in missingnames
            push!(_namemap.parameters, name)
        end
    end

    return Chains(chains.value, chains.logevidence, _namemap, chains.info)
end

function _clean_sections(c::Chains, sections::Union{Vector{Symbol}, Symbol})
    sections = sections isa AbstractArray ? sections : [sections]
    ks = collect(keys(c.name_map))
    return ks ∩ sections
end

#################### Concatenation ####################

Base.cat(c::Chains, cs::Chains...; dims = Val(1)) = _cat(dims, c, cs...)
Base.cat(c::T, cs::T...; dims = Val(1)) where T<:Chains = _cat(dims, c, cs...)

Base.vcat(c::Chains, cs::Chains...) = _cat(Val(1), c, cs...)
Base.vcat(c::T, cs::T...) where T<:Chains = _cat(Val(1), c, cs...)

Base.hcat(c::Chains, cs::Chains...) = _cat(Val(2), c, cs...)
Base.hcat(c::T, cs::T...) where T<:Chains = _cat(Val(2), c, cs...)

AbstractMCMC.chainscat(c::Chains, cs::Chains...) = _cat(Val(3), c, cs...)

_cat(dim::Int, cs::Chains...) = _cat(Val(dim), cs...)

function _cat(::Val{1}, c1::Chains, args::Chains...)
    # check inputs
    thin = step(c1)
    all(c -> step(c) == thin, args) || throw(ArgumentError("chain thinning differs"))
    nms = names(c1)
    all(c -> names(c) == nms, args) || throw(ArgumentError("chain names differ"))
    chns = chains(c1)
    all(c -> chains(c) == chns, args) || throw(ArgumentError("sets of chains differ"))

    # concatenate all chains
    data = mapreduce(c -> c.value.data, vcat, args; init = c1.value.data)
    value = AxisArray(data;
                      iter = range(first(c1); length = size(data, 1), step = thin),
                      var = nms,
                      chain = chns)

    return Chains(value, missing, c1.name_map, c1.info)
end

function _cat(::Val{2}, c1::Chains, args::Chains...)
    # check inputs
    rng = range(c1)
    all(c -> range(c) == rng, args) || throw(ArgumentError("chain ranges differ"))
    chns = chains(c1)
    all(c -> chains(c) == chns, args) || throw(ArgumentError("sets of chains differ"))

    # combine names and sections of parameters
    nms = names(c1)
    n = length(nms)
    for c in args
        nms = union(nms, names(c))
        n += length(names(c))
        n == length(nms) || throw(ArgumentError("non-unique parameter names"))
    end

    name_map = mapreduce(c -> c.name_map, merge_union, args; init = c1.name_map)

    # concatenate all chains
    data = mapreduce(c -> c.value.data, hcat, args; init = c1.value.data)
    value = AxisArray(data; iter = rng, var = nms, chain = chns)

    return Chains(value, missing, name_map, c1.info)
end

function _cat(::Val{3}, c1::Chains, args::Chains...)
    # check inputs
    rng = range(c1)
    all(c -> range(c) == rng, args) || throw(ArgumentError("chain ranges differ"))
    nms = names(c1)
    all(c -> names(c) == nms, args) || throw(ArgumentError("chain names differ"))

    # concatenate all chains
    data = mapreduce(c -> c.value.data, (x, y) -> cat(x, y; dims = 3), args;
                     init = c1.value.data)
    value = AxisArray(data; iter = rng, var = nms, chain = 1:size(data, 3))

    return Chains(value, missing, c1.name_map, c1.info)
end

function pool_chain(c::Chains)
    data = c.value.data
    pool_data = reshape(permutedims(data, [1, 3, 2]), :, size(data, 2), 1)
    return Chains(pool_data, names(c), c.name_map; info=c.info)
end

function set_names(c::Chains, d::Dict; sorted::Bool=true)
    # Set new parameter names.
    params = names(c)
    new_params = replace(params, d...)

    # Iterate through each pair in the name map to
    # create a new one.
    new_map = Dict()
    for (section, parameters) in pairs(c.name_map)
        new_map[section] = replace(parameters, d...)
    end

    # Return a new chains object.
    return Chains(
        c.value.data,
        new_params,
        new_map;
        info=c.info,
        evidence=c.logevidence,
        sorted=sorted,
    )
end
