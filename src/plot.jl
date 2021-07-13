@shorthands meanplot
@shorthands autocorplot
@shorthands mixeddensity
@shorthands pooleddensity
@shorthands traceplot
@shorthands corner
@shorthands violinplot

struct _TracePlot; c; val; end
struct _MeanPlot; c; val;  end
struct _DensityPlot; c; val;  end
struct _HistogramPlot; c; val;  end
struct _AutocorPlot; lags; val;  end
struct _ViolinPlot; par; val; end

# define alias functions for old syntax
const translationdict = Dict(
                        :traceplot => _TracePlot,
                        :meanplot => _MeanPlot,
                        :density => _DensityPlot,
                        :histogram => _HistogramPlot,
                        :autocorplot => _AutocorPlot,
                        :pooleddensity => _DensityPlot,
                        :violinplot => _ViolinPlot
                      )

const supportedplots = push!(collect(keys(translationdict)), :mixeddensity, :corner)

@recipe f(c::Chains, s::Symbol) = c, [s]

@recipe function f(
    chains::Chains, i::Int;
    colordim = :chain,
    barbounds = (-Inf, Inf),
    maxlag = nothing,
    append_chains = false,
    sections = chains.name_map[:parameters],
    combined = true
)
    st = get(plotattributes, :seriestype, :traceplot)
    c = append_chains || st == :pooleddensity ? pool_chain(chains) : chains

    if colordim == :parameter
        title --> "Chain $(MCMCChains.chains(c)[i])"
        label --> string.(names(c))
        val = c.value[:, :, i]
    elseif colordim == :chain
        title --> string(names(c)[i])
        label --> map(x -> "Chain $x", MCMCChains.chains(c))
        val = c.value[:, i, :]
    else
        throw(ArgumentError("`colordim` must be one of `:chain` or `:parameter`"))
    end

    if st == :mixeddensity || st == :pooleddensity
        discrete = indiscretesupport(c, barbounds)
        st = if colordim == :chain
            discrete[i] ? :histogram : :density
        else
            # NOTE: It might make sense to overlay histograms and density plots here.
            :density
        end
        seriestype := st
    end

    if st == :autocorplot
        lags = 0:(maxlag === nothing ? round(Int, 10 * log10(length(range(c)))) : maxlag)
        ac = autocor(c; sections = nothing, lags = lags)
        ac_mat = convert(Array, ac)
        val = colordim == :parameter ? ac_mat[:, :, i]' : ac_mat[i, :, :]
        _AutocorPlot(lags, val)
    elseif st ∈ supportedplots
        translationdict[st](c, val)
    else
        range(c), val
    end

    total_chains = i
    if st == :violinplot
        n_iter, n_par, n_chains = size(chains)
        if combined
            colordim := :chain
            par = string.(reshape(repeat(sections, inner = n_iter), n_iter, n_par))[:,i]
            val = Array(chains)[:,i]
            _ViolinPlot(par, val)
        elseif combined == false
            if colordim == :chain
                par_names = ["$(sections[i]).Chain $j" for i in 1:n_par, j in 1:n_chains]
                pars = string.(reshape(repeat(vec(par_names), inner = n_iter), (n_iter, n_par, n_chains)))
                val = chains.value[:,i,:]
                par = pars[:,i,:]
            elseif colordim == :parameter
                par_vec = repeat(sections, inner = n_iter)
                pars = string.(reshape(repeat(par_vec, n_chains, 1), (n_iter, n_par, n_chains)))
                val = chains.value[:,:,i]
                par = pars[:,:,i]
                label --> string.(names(c))
            else
                throw(ArgumentError("`colordim` must be one of `:chain` or `:parameter`"))
            end
            _ViolinPlot(par, val)
        else
            throw(ArgumentError("In `ViolinPlots` `Chains` can be combined or separated "))
        end
    elseif st ∈ supportedplots
        translationdict[st](c, val)
    else
        range(c), val
    end
end

@recipe function f(p::_DensityPlot)
    xaxis --> "Sample value"
    yaxis --> "Density"
    trim --> true
    [collect(skipmissing(p.val[:,k])) for k in 1:size(p.val, 2)]
end

@recipe function f(p::_HistogramPlot)
    xaxis --> "Sample value"
    yaxis --> "Frequency"
    fillalpha --> 0.7
    bins --> 25
    trim --> true
    [collect(skipmissing(p.val[:,k])) for k in 1:size(p.val, 2)]
end

@recipe function f(p::_MeanPlot)
    seriestype := :path
    xaxis --> "Iteration"
    yaxis --> "Mean"
    range(p.c), cummean(p.val)
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
    range(p.c), p.val
end

@recipe function f(
    chains::Chains,
    parameters::AbstractVector{Symbol};
    colordim = :chain
)
    colordim != :chain &&
        error("Symbol names are interpreted as parameter names, only compatible with ",
              "`colordim = :chain`")

    ret = indexin(parameters, names(chains))
    any(y === nothing for y in ret) && error("Parameter not found")

    return chains, Int.(ret)
end

@recipe function f(
    chains::Chains,
    parameters::AbstractVector{<:Integer} = Int[];
    sections = _default_sections(chains),
    width = 500,
    height = 250,
    colordim = :chain,
    append_chains = false
)
    _chains = isempty(parameters) ? Chains(chains, _clean_sections(chains, sections)) : chains
    c = append_chains ? pool_chain(_chains) : _chains
    ptypes = get(plotattributes, :seriestype, (:traceplot, :mixeddensity))
    ptypes = ptypes isa Symbol ? (ptypes,) : ptypes
    @assert all(ptype -> ptype ∈ supportedplots, ptypes)
    ntypes = length(ptypes)
    nrows, nvars, nchains = size(c)
    isempty(parameters) && (parameters = colordim == :chain ? (1:nvars) : (1:nchains))
    N = length(parameters)

    if :corner ∉ ptypes
        size --> (ntypes*width, N*height)
        legend --> false

        multiple_plots = N * ntypes > 1
        if multiple_plots
            layout := (N, ntypes)
        end

        i = 0
        for par in parameters
            for ptype in ptypes
                i += 1

                @series begin
                    if multiple_plots
                        subplot := i
                    end
                    colordim := colordim
                    seriestype := ptype
                    c, par
                end
            end
        end
    else
        ntypes > 1 && error(":corner is not compatible with multiple seriestypes")
        Corner(c, names(c)[parameters])
    end
end

struct Corner
    c
    parameters
end

@recipe function f(corner::Corner)
    label --> permutedims(corner.parameters)
    compact --> true
    size --> (600, 600)
    ar = collect(Array(corner.c.value[:, corner.parameters,i]) for i in chains(corner.c))
    RecipesBase.recipetype(:cornerplot, vcat(ar...))
end

@recipe function f(p::_ViolinPlot)
    @series begin
        seriestype := :violin
        p.par, p.val
    end

    @series begin
        seriestype := :boxplot
        bar_width --> 0.1
        linewidth --> 2
        fillalpha --> 0.8
        p.par, p.val
    end
end
