"""
    ess_rhat(chains::Chains; duration=compute_duration, kwargs...)

Estimate the effective sample size and the potential scale reduction.

ESS per second options include `duration=MCMCChains.compute_duration` (the default)
and `duration=MCMCChains.wall_duration`.
"""
function MCMCDiagnosticTools.ess_rhat(
    chains::Chains;
    sections = _default_sections(chains), duration = compute_duration, kwargs...
)
    # Subset the chain
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Estimate the effective sample size and rhat
    ess, rhat = MCMCDiagnosticTools.ess_rhat(_chains.value.data; kwargs...)

    # Calculate ESS/minute if available
    dur = duration(chains)

    # Convert to NamedTuple
    nt = if dur === missing
        merge((parameters = names(_chains),), (ess = ess, rhat = rhat))
    else
        merge((parameters = names(_chains),), (ess = ess, rhat = rhat, ess_per_sec=ess/dur))
    end

    return ChainDataFrame("ESS", nt)
end
