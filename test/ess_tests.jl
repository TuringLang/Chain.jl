using MCMCChains
using FFTW

using Random
using Statistics
using Test

@testset "ESS per second" begin
    c1 = Chains(randn(100,5, 1), info = (start_time=time(), stop_time = time()+1))
    c2 = Chains(randn(100,5, 1), info = (start_time=time()+1, stop_time = time()+2))
    c = chainscat(c1, c2)

    wall = MCMCChains.wall_duration(c)
    compute = MCMCChains.compute_duration(c)

    @test round(wall, digits=1) ≈ round(c2.info.stop_time - c1.info.start_time, digits=1)
    @test compute ≈ (MCMCChains.compute_duration(c1) + MCMCChains.compute_duration(c2))

    s = ess_rhat(c)
    @test length(s[:,:ess_per_sec]) == 5
    @test all(map(!ismissing, s[:,:ess_per_sec]))
end

@testset "ESS and R̂ (chains)" begin
    x = rand(10_000, 40, 10)
    chain = Chains(x)

    for method in (ESSMethod(), FFTESSMethod(), BDAESSMethod())
        # analyze chain
        ess_df = ess_rhat(chain; method = method)

        # analyze array
        ess_array, rhat_array = ess_rhat(x; method = method)

        @test ess_df[:,2] == ess_array
        @test ess_df[:,3] == rhat_array
    end
end

@testset "ESS and R̂ (single sample)" begin # check that issue #137 is fixed
    val = rand(1, 5, 3)
    chain = Chains(val, ["a", "b", "c", "d", "e"])

    for method in (ESSMethod(), FFTESSMethod(), BDAESSMethod())
        # analyze chain
        ess_df = ess_rhat(chain; method = method)

        # analyze array
        ess_array, rhat_array = ess_rhat(val; method = method)

        @test ismissing(ess_df[:,2][1]) # since min(maxlag, niter - 1) = 0
        @test ismissing(ess_df[:,3][1])
    end
end
