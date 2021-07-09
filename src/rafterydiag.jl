function MCMCDiagnosticTools.rafterydiag(
    chains::Chains; sections = _default_sections(chains), kwargs...
)
    # Subset the chain.
    _chains = Chains(chains, _clean_sections(chains, sections))

    # Compute statistics for each chain.
    chains_array = _chains.value.data
    results = map(axes(chains_array, 3)) do chain
        vec_of_namedtuples = map(axes(chains_array, 2)) do param
            return MCMCDiagnosticTools.rafterydiag(
                cskip(@view(_chains.value.data[:, param, chain])); kwargs...
            )
        end
        namedtuple_of_vecs = Tables.columntable(vec_of_namedtuples)
        return namedtuple_of_vecs
    end

    # Create data frames.
    parameters = (parameters = names(_chains),)
    dfs = [
        ChainDataFrame(
            "Raftery and Lewis diagnostic - Chain $i", merge(parameters, result)
        )
        for (i, result) in enumerate(results)
    ]

    return dfs
end
