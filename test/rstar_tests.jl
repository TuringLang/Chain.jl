using MCMCChains
using Tables
using MLJ, MLJModels
using Test

val = rand(1000, 8, 4)
colnames = ["a", "b", "c", "d", "e", "f", "g", "h"]
internal_colnames = ["c", "d", "e", "f", "g", "h"]
chn = Chains(val, colnames, Dict(:internals => internal_colnames))

model = @load XGBoostClassifier()

@testset "RStarTable test" begin
    t = MCMCChains.RStarTable(chn)

    @test Tables.istable(typeof(t))
    @test Tables.columnaccess(typeof(t))
    @test Tables.matrix(t) === t.data
    @test t[[1,2]] isa MCMCChains.RStarTable
end

@testset "R star test" begin

    # Compute R* statistic for a mixed chain.
    R = rstar(chn, model; iterations = 10)

    # Resulting R value should be close to one, i.e. the classifier does not perform better than random guessing.
    @test mean(R) ≈ 1 atol=0.1

    # Compute R* statistic for a non-mixed chain.
    niter = 1000
    val = hcat(sin.(1:niter), cos.(1:niter))
    val = cat(val, hcat(cos.(1:niter)*100, sin.(1:niter)*100), dims=3)
    chn_notmixed = Chains(val)

    # Restuling R value should be close to two, i.e. the classifier should be able to learn an almost perfect decision boundary between chains.
    R = rstar(chn_notmixed, model; iterations = 10)
    @test mean(R) ≈ 2 atol=0.1
end