@shorthands meanplot
@shorthands autocorplot
@shorthands mixeddensity
@shorthands pooleddensity
@shorthands traceplot
@shorthands corner
@shorthands ppcplot

struct _TracePlot; c; val; end
struct _MeanPlot; c; val;  end
struct _DensityPlot; c; val;  end
struct _HistogramPlot; c; val;  end
struct _AutocorPlot; lags; val;  end
struct _PPCPlot; y_obs; y_pred; ymean_pred; end

# define alias functions for old syntax
const translationdict = Dict(
                        :traceplot => _TracePlot,
                        :meanplot => _MeanPlot,
                        :density => _DensityPlot,
                        :histogram => _HistogramPlot,
                        :autocorplot => _AutocorPlot,
                        :pooleddensity => _DensityPlot,
                        :ppcplot => _PPCPlot
                      )

const supportedplots = push!(collect(keys(translationdict)), :mixeddensity, :corner)

@recipe f(c::Chains, s::Symbol) = c, [s]

@recipe function f(
    chains::Chains, i::Int;
    colordim = :chain,
    barbounds = (-Inf, Inf),
    maxlag = nothing,
    append_chains = false
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

@recipe function f(
    yobs_data,
    ypred_data::Chains,
    check_type = :posterior,
    n_samples::Int = 50;
    yvar_name::AbstractVector{<:String}
    )

    st = get(plotattributes, :seriestype, :traceplot)

    if st == :ppcplot
        N = n_samples <= size(ypred_data)[1] ? n_samples : size(ypred_data)[1]
        index = sample(1:size(ypred_data)[1], N, replace = false, ordered = true)
        if check_type == :posterior
            title --> "Posterior predictive check"
        elseif check_type == :prior
            title --> "Prior predictive check"
        else
            throw(ArgumentError("`check_type` must be one of `prior` or `posterior`"))
        end

        if ndims(yobs_data) == 1
            #y_obs = yobs_data
            n_yvar = 1
            n_yval = size(yobs_data)[1] #y values per sample
            predictions = ypred_data.value.data[index,:,:]
            ymean_pred = vec(mean(ypred_data.value.data, dims = 1))
            @series begin
                seriestype := :density
                label --> "Y obs"
                y_obs = yobs_data
            end
            @series begin
                seriestype := :density
                label --> "Y mean"
                ymean_pred
            end
            for i in 1:N
                @series begin
                    seriestype := :density
                    seriesalpha --> 0.2
                    linecolor --> "#BBBBBB"
                    label --> nothing
                    y_pred = predictions[i,:,:]
                end
            end
        elseif ndims(yobs_data) >= 1 ## ordered by columns
            n_yval = size(yobs_data)[1]
            n_yvar = size(yobs_data)[2] ## number of dependent variables
            #n_iter = size(ypred_data)[1]
            y_obs = yobs_data ## observed data for y
            ymean_pred = reshape(mean(ypred_data.value.data, dims = 1), (n_yvar, n_yval))
            #predictions = reshape(ypred_data.value.data, (n_iter,n_yval,n_yvar))[index,:,:]
            #sections = [MCMCChains.group(ypred_data, Symbol(yvar_name[i])).value.data[index,:,:]
            #            for i in 1:n_yvar]
            layout := (1, n_yvar)
            k = 0
            for j in 1:n_yvar
                sections = MCMCChains.group(ypred_data, Symbol(yvar_name[j]))
                predictions = sections.value.data[index,:,:]
                #sections = dropdims(group(ypred_data, yvar_name[j]).value.data, dims = 3)
                #predictions = reshape(ypred_data.value.data, (n_iter,n_yval,n_yvar))[index,:,:]
                k += 1
                @series begin
                    subplot := k
                    seriestype := :density
                    label --> "Yobs"
                    y_obs[:,j]
                end
                @series begin
                    subplot := k
                    seriestype := :density
                    label --> "Ymean"
                    ymean_pred[j,:]
                end
                for i in 1:N
                    #y_pred = sections[j][i,:,:]
                    @series begin
                        subplot := k
                        seriestype := :density
                        seriesalpha --> 0.2
                        linecolor --> "#BBBBBB"
                        label --> nothing
                        y_pred = predictions[i,:,:]
                    end
                end
            end
        else
            throw(ArgumentError("Observed data must have `dim > 1`"))
        end
    else

    end
end

@recipe function f(p::_PPCPlot)
    #title --> "Posterior predictive check"
    p.y_obs, p.y_pred
end