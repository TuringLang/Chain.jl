"""
    discretediag(chains::Chains{<:Real}; sections, kwargs...)

Discrete diagnostic where `method` can be
`[:weiss, :hangartner, :DARBOOT, MCBOOT, :billinsgley, :billingsleyBOOT]`.
"""
function MCMCDiagnosticTools.discretediag(
    chains::Chains{<:Real};
    sections = _default_sections(chains),
    kwargs...
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Compute statistics.
    between_chain_vals, within_chain_vals = MCMCDiagnosticTools.discretediag(
        _chains.value.data; kwargs...
    )

    # Create dataframes
    parameters = (parameters = names(_chains),)
    between_chain_df = ChainDataFrame(
        "Chisq diagnostic - Between chains", merge(parameters, between_chain_vals),
    )
    within_chain_dfs = map(1:size(_chains, 3)) do i
        vals = map(val -> val[:, i], within_chain_vals)
        return ChainDataFrame("Chisq diagnostic - Chain $i", merge(parameters, vals))
    end
    dfs = vcat(between_chain_df, within_chain_dfs)

    return dfs
end
