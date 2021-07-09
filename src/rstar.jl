"""
    rstar([rng ,]classifier, chains::Chains; kwargs...)

Compute the ``R^*`` convergence diagnostic of MCMC for the `chains`.

The keyword arguments supported here are the same as those in `rstar` for arrays of samples
and chain indices.

# Example

```jldoctest rstar; setup = :(using Random; Random.seed!(100))
julia> using MLJBase, MLJXGBoostInterface, Statistics

julia> chains = Chains(fill(4, 100, 2, 3));
```

One can compute the distribution of the ``R^*`` statistic with the probabilistic classifier.

```jldoctest rstar
julia> distribution = rstar(XGBoostClassifier(), chains);

julia> isapprox(mean(distribution), 1; atol=0.1)
true
```

For deterministic classifiers, a single ``R^*`` statistic is returned.

```jldoctest rstar
julia> @pipeline XGBoostClassifier name = XGBoostDeterministic operation = predict_mode;

julia> value = rstar(XGBoostDeterministic(), chains);

julia> isapprox(value, 1; atol=0.1)
true
```
"""
function MCMCDiagnosticTools.rstar(
    classif::MLJModelInterface.Supervised, chn::Chains; kwargs...
)
    return MCMCDiagnosticTools.rstar(Random.GLOBAL_RNG, classif, chn; kwargs...)
end

function MCMCDiagnosticTools.rstar(
    rng::Random.AbstractRNG, classif::MLJModelInterface.Supervised, chn::Chains; kwargs...
)
    nchains = size(chn, 3)
    nchains <= 1 && throw(DimensionMismatch())

    # collect data
    x = Array(chn)
    y = repeat(chains(chn); inner = size(chn,1))

    return MCMCDiagnosticTools.rstar(rng, classif, x, y; kwargs...)
end
