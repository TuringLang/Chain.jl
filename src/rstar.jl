"""
    rstar([rng ,]classif::Supervised, chains::Chains; kwargs...)
    rstar([rng ,]classif::Supervised, x::AbstractMatrix, y::AbstractVector; kwargs...)

Compute the ``R^*`` convergence diagnostic of MCMC.

This implementation is an adaption of Algorithm 1 & 2, described by [^LambertVehtari2020].
Note that the correctness of the statistic depends on the convergence of the classifier used
internally in the statistic. You can inspect the training of the classifier by adjusting the
verbosity level.

[^LambertVehtari2020]: Lambert & Vehtari (2020). ``R^*``: A robust MCMC convergence diagnostic with uncertainty using gradient-boosted machines. arXiv preprint <https://arxiv.org/abs/2003.07900>.

# Keyword Arguments
* `subset = 0.8` ... Subset used to train the classifier, i.e. 0.8 implies 80% of the samples are used.
* `iterations = 10` ... Number of iterations used to estimate the statistic. If the classifier is not probabilistic, i.e. does not return class probabilities, it is advisable to use a value of one.
* `verbosity = 0` ... Verbosity level used during fitting of the classifier.

# Example
```jldoctest rstar; output = false, filter = r".*"s
using MLJBase, MLJXGBoostInterface

chn = Chains(fill(4, 100, 2, 3))

Rs = rstar(XGBoostClassifier(), chn; iterations=20)
R = round(mean(Rs); digits=0)

# output

1.0
```
"""
function rstar(rng::Random.AbstractRNG, classif::MLJModelInterface.Supervised, x::AbstractMatrix, y::AbstractVector{Int}; iterations = 10, subset = 0.8, verbosity = 0)

    size(x,1) != length(y) && throw(DimensionMismatch())
    iterations >= 1 && ArgumentError("Number of iterations has to be positive!")

    if iterations > 1 && classif isa MLJModelInterface.Deterministic
        @warn("Classifier is not a probabilistic classifier but number of iterations is > 1.")
    elseif iterations == 1 && classif isa MLJModelInterface.Probabilistic
        @warn("Classifier is probabilistic but number of iterations is equal to one.")
    end

    N = length(y)
    K = length(unique(y))

    # randomly sub-select training and testing set
    Ntrain = round(Int, N*subset)
    Ntest = N - Ntrain

    ids = Random.randperm(rng, N)
    train_ids = view(ids, 1:Ntrain)
    test_ids = view(ids, (Ntrain+1):N)

    # train classifier using XGBoost
    fitresult, _ = MLJModelInterface.fit(classif, verbosity, Tables.table(x[train_ids,:]), MLJModelInterface.categorical(y[train_ids]))

    xtest = Tables.table(x[test_ids,:])
    ytest = view(y, test_ids)

    Rstats = map(i -> K*rstar_score(rng, classif, fitresult, xtest, ytest), 1:iterations)
    return Rstats
end

function rstar(classif::MLJModelInterface.Supervised, x::AbstractMatrix, y::AbstractVector{Int}; kwargs...)
    rstar(Random.GLOBAL_RNG, classif, x, y; kwargs...)
end

function rstar(classif::MLJModelInterface.Supervised, chn::Chains; kwargs...)
    return rstar(Random.GLOBAL_RNG, classif, chn; kwargs...)
end

function rstar(rng::Random.AbstractRNG, classif::MLJModelInterface.Supervised, chn::Chains; kwargs...)
    nchains = size(chn, 3)
    nchains <= 1 && throw(DimensionMismatch())

    # collect data
    x = Array(chn)
    y = repeat(chains(chn); inner = size(chn,1))

    return rstar(rng, classif, x, y; kwargs...)
end

function rstar_score(rng::Random.AbstractRNG, classif::MLJModelInterface.Probabilistic, fitresult, xtest, ytest)
    pred = get.(rand.(Ref(rng), MLJModelInterface.predict(classif, fitresult, xtest)))
    return mean(((p,y),) -> p == y, zip(pred, ytest))
end

function rstar_score(rng::Random.AbstractRNG, classif::MLJModelInterface.Deterministic, fitresult, xtest, ytest)
    pred = MLJModelInterface.predict(classif, fitresult, xtest)
    return mean(((p,y),) -> p == y, zip(pred, ytest))
end
