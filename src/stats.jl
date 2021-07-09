#################### Posterior Statistics ####################

"""
    autocor(chains[; lags = [1, 5, 10, 50], demean = true, append_chains = true, kwargs...])

Compute the autocorrelation of each parameter for the chain.

Setting `append_chains=false` will return a vector of dataframes containing the
autocorrelations for each chain.
"""
function autocor(
    chains::Chains;
    lags::AbstractVector{<:Integer} = [1, 5, 10, 50],
    demean::Bool = true,
    append_chains = false,
    kwargs...
)
    funs = Function[]
    func_names = @. Symbol("lag ", lags)
    for i in lags
        push!(funs, x -> autocor(x, [i], demean=demean)[1])
    end

    return summarize(
        chains, funs...;
        func_names = func_names,
        append_chains = append_chains,
        name = "Autocorrelation",
        kwargs...
    )
end

"""
    cor(chains[; sections, append_chains = true, kwargs...])

Compute the Pearson correlation matrix for the chain.

Setting `append_chains=false` will return a vector of dataframes containing a correlation
matrix for each chain.
"""
function cor(
    chains::Chains;
    sections = _default_sections(chains),
    append_chains = true,
    kwargs...
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Obstain names of parameters.
    names_of_params = names(_chains)

    if append_chains
        df = chaindataframe_cor("Correlation", names_of_params, to_matrix(_chains))
        return df
    else
        vector_of_df = [
            chaindataframe_cor(
                "Correlation - Chain $i", names_of_params, data
            )
            for (i, data) in enumerate(to_vector_of_matrices(_chains))
        ]
        return vector_of_df
    end
end

function chaindataframe_cor(name, names_of_params, chains::AbstractMatrix; kwargs...)
    # Compute the correlation matrix.
    cormat = cor(chains)

    # Summarize the results in a named tuple.
    nt = (; parameters = names_of_params,
          zip(names_of_params, (cormat[:, i] for i in axes(cormat, 2)))...)

    # Create a ChainDataFrame.
    return ChainDataFrame(name, nt; kwargs...)
end

"""
    changerate(chains[; sections, append_chains = true, kwargs...])

Compute the change rate for the chain.

Setting `append_chains=false` will return a vector of dataframes containing the change
rates for each chain.
"""
function changerate(
    chains::Chains{<:Real};
    sections = _default_sections(chains),
    append_chains = true,
    kwargs...
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Obstain names of parameters.
    names_of_params = names(_chains)

    if append_chains
        df = chaindataframe_changerate("Change Rate", names_of_params, _chains.value.data)
        return df
    else
        vector_of_df = [
            chaindataframe_changerate(
                "Change Rate - Chain $i", names_of_params, data
            )
            for (i, data) in enumerate(to_vector_of_matrices(_chains))
        ]
        return vector_of_df
    end
end

function chaindataframe_changerate(name, names_of_params, chains; kwargs...)
    # Compute the change rates.
    changerates, mvchangerate = changerate(chains)

    # Summarize the results in a named tuple.
    nt = (; zip(names_of_params, changerates)..., multivariate = mvchangerate)

    # Create a ChainDataFrame.
    return ChainDataFrame(name, nt; kwargs...)
end

changerate(chains::AbstractMatrix{<:Real}) = changerate(reshape(chains, Val(3)))
function changerate(chains::AbstractArray{<:Real,3})
    niters, nparams, nchains = size(chains)

    changerates = zeros(nparams)
    mvchangerate = 0.0

    for chain in 1:nchains, iter in 2:niters
        isanychanged = false

        for param in 1:nparams
            # update if the sample is different from the one in the previous iteration
            if chains[iter-1, param, chain] != chains[iter, param, chain]
                changerates[param] += 1
                isanychanged = true
            end
        end

        mvchangerate += isanychanged
    end

    factor = nchains * (niters - 1)
    changerates ./= factor
    mvchangerate /= factor

    changerates, mvchangerate
end

describe(c::Chains; args...) = describe(stdout, c; args...)

"""
    describe(io, chains[;
             q = [0.025, 0.25, 0.5, 0.75, 0.975],
             etype = :bm,
             kwargs...])

Print the summary statistics and quantiles for the chain.
"""
function describe(
    io::IO,
    chains::Chains;
    q = [0.025, 0.25, 0.5, 0.75, 0.975],
    etype = :bm,
    kwargs...
)
    dfs = vcat(summarystats(chains; etype = etype, kwargs...),
               quantile(chains; q = q, kwargs...))
    return dfs
end

function _hpd(x::AbstractVector{<:Real}; alpha::Real=0.05)
    n = length(x)
    m = max(1, ceil(Int, alpha * n))

    y = sort(x)
    a = y[1:m]
    b = y[(n - m + 1):n]
    _, i = findmin(b - a)

    return [a[i], b[i]]
end

function hpd(chn::Chains; alpha::Real=0.05, kwargs...)
    labels = [:lower, :upper]
    l(x) = _hpd(x, alpha=alpha)[1]
    u(x) = _hpd(x, alpha=alpha)[2]
    return summarize(chn, l, u; name = "HPD", func_names = labels, kwargs...)
end

"""
    quantile(chains[; q = [0.025, 0.25, 0.5, 0.75, 0.975], append_chains = true, kwargs...])

Compute the quantiles for each parameter in the chain.

Setting `append_chains=false` will return a vector of dataframes containing the quantiles
for each chain.
"""
function quantile(
    chains::Chains;
    q::AbstractVector = [0.025, 0.25, 0.5, 0.75, 0.975],
    append_chains = true,
    kwargs...
)
    # compute quantiles
    funs = Function[]
    func_names = @. Symbol(100 * q, :%)
    for i in q
        push!(funs, x -> quantile(cskip(x), i))
    end

    return summarize(
        chains, funs...;
        func_names = func_names,
        append_chains = append_chains,
        name = "Quantiles",
        kwargs...
    )
end


"""
    function summarystats(
        chains;
        sections = _default_sections(chains),
        append_chains= true,
        method::AbstractESSMethod = ESSMethod(),
        maxlag = 250,
        etype = :bm,
        kwargs...
    )

Compute the mean, standard deviation, naive standard error, Monte Carlo standard error,
and effective sample size for each parameter in the chain.

Setting `append_chains=false` will return a vector of dataframes containing the summary
statistics for each chain.

When estimating the effective sample size, autocorrelations are computed for at most `maxlag` lags.
"""
function summarystats(
    chains::Chains;
    sections = _default_sections(chains),
    append_chains::Bool = true,
    method::MCMCDiagnosticTools.AbstractESSMethod = ESSMethod(),
    maxlag = 250,
    etype = :bm,
    kwargs...
)
    # Store everything.
    funs = [mean∘cskip, std∘cskip, sem∘cskip, x -> MCMCDiagnosticTools.mcse(cskip(x); method=etype, kwargs...)]
    func_names = [:mean, :std, :naive_se, :mcse]

    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Calculate ESS separately.
    ess_df = MCMCDiagnosticTools.ess_rhat(_chains; sections = nothing, method = method, maxlag = maxlag)

    # Summarize.
    summary_df = summarize(
        _chains, funs...;
        func_names = func_names,
        append_chains = append_chains,
        additional_df = ess_df,
        name = "Summary Statistics",
        sections = nothing
    )

    return summary_df
end

"""
    mean(chains[, params; kwargs...])

Calculate the mean of a chain.
"""
function mean(chains::Chains; kwargs...)
    # Store everything.
    funs = [mean∘cskip]
    func_names = [:mean]

    # Summarize.
    summary_df = summarize(
        chains, funs...;
        func_names = func_names,
        name = "Mean",
        kwargs...
    )

    return summary_df
end

mean(chn::Chains, syms) = mean(chn[:, syms, :])
# resolve method ambiguity with `mean(f, ::AbstractArray)`
mean(chn::Chains, syms::AbstractVector) = mean(chn[:, syms, :])
