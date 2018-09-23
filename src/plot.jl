
@shorthands meanplot
@shorthands autocorplot
@shorthands mixeddensity
@shorthands traceplot

struct _TracePlot; c; val; end
struct _MeanPlot; c; val;  end
struct _DensityPlot; c; val;  end
struct _HistogramPlot; c; val;  end
struct _AutocorPlot; lags; val;  end

# define alias functions for old syntax
const translationdict = Dict(
                        :traceplot => _TracePlot,
                        :meanplot => _MeanPlot,
                        :density => _DensityPlot,
                        :histogram => _HistogramPlot,
                        :autocorplot => _AutocorPlot
                      )

const supportedplots = push!(collect(keys(translationdict)), :mixeddensity)

@recipe function f(c::AbstractChains, i::Int; colordim = :chain, barbounds = (0, Inf), maxlag = nothing)
    st = get(plotattributes, :seriestype, :traceplot)

    if colordim == :parameter
        title := "Chain $(c.chains[i])"
        labels := c.names
        val = c.value[:, :, i]
    elseif colordim == :chain
        title := c.names[i]
        labels := map(k -> "Chain $(c.chains[k])", 1:size(c)[3])
        val = c.value[:, i, :]
    else
        throw(ArgumentError("`colordim` must be one of `:chain` or `:parameter`"))
    end

    if st == :mixeddensity
        discrete = MCMCChain.indiscretesupport(c, barbounds)
        st = (discrete[i] && colordim == :chain) ? :histogram : :density
        seriestype := st
    end

    if st == :autocorplot
        lags = 0:(maxlag === nothing ? round(Int, 10 * log10(length(c.range))) : maxlag)
        ac = MCMCChain.autocor(c, lags=collect(lags))
        val = colordim == :parameter ? ac.value[:, :, i]' : ac.value[i, :, :]
        _AutocorPlot(lags, val)
    elseif st ∈ supportedplots
        translationdict[st](c, val)
    else
        c.range, val
    end
end

@recipe function f(p::_DensityPlot)
    xaxis --> "Sample value"
    yaxis --> "Density"
    [collect(skipmissing(p.val[:,k])) for k in axes(p.val, 2)]
end

@recipe function f(p::_HistogramPlot)
    xaxis --> "Sample value"
    yaxis --> "Frequency"
    fillalpha --> 0.7
    [collect(skipmissing(p.val[:,k])) for k in axes(p.val, 2)]
end

@recipe function f(p::_MeanPlot)
    seriestype := :path
    xaxis --> "Iteration"
    yaxis --> "Mean"
    p.c.range, MCMCChain.cummean(p.val)
end

@recipe function f(p::_AutocorPlot)
    seriestype := :path
    xaxis --> "Lag"
    yaxis --> "Autocorrelation"
    p.lags, p.val
end

@recipe function f(p::_TracePlot)
    seriestype := :path
    xaxis --> "Iteration"
    yaxis --> "Sample value"
    p.c.range, p.val
end

@recipe function f(c::MCMCChain.AbstractChains;
                   width = 500,
                   height = 250,
                   colordim = :chain
                  )

    ptypes = get(plotattributes, :seriestype, (:traceplot, :mixeddensity))

    # sanity check
    ptypes = ptypes isa AbstractVector || ptypes isa Tuple ? ptypes : (ptypes,)
    @assert all(map(ptype -> ptype ∈ supportedplots, ptypes))

    nrows, nvars, nchains = size(c.value)
    ntypes = length(ptypes)
    N = colordim == :chain ? nvars : nchains
    layout := (N, ntypes)
    size --> (ntypes*width, N*height)
    indices = reshape(1:N*ntypes, ntypes, N)'

    legend --> false

    for (j, ptype) in enumerate(ptypes)
        for i in 1:N
            @series begin
                subplot := indices[i, j]
                colordim := colordim
                seriestype := ptype
                c, i
            end
        end
    end
end

function plot(c::AbstractChains, psyms::Vector{Symbol}; args...)
    @assert all(map(psym -> haskey(translationdict, psym), psyms))
    @warn "This syntax is deprecated, please use plot(c; seriestype = $psyms) instead."
    return plot(c; seriestype = map(psym -> translationdict[psym], psyms))
end

function plot(c::AbstractChains, psym::Symbol; args...)
    @assert haskey(translationdict, psym)
    @warn "This syntax is deprecated, please use plot(c, seriestype = $(translationdict[psym])) instead"
    return plot(c; seriestype = psym, args...)
end
