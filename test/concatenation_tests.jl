using MCMCChains
using Test

@testset "merge_union" begin
    @test @inferred(MCMCChains.merge_union((a = [1], b = [3.0]), (c = [3], a = [2.5]))) ==
        (a = [1.0, 2.5], b = [3.0], c = [3])
end

@testset "concatenation tests" begin
    v1 = rand(500, 5, 1)
    v2 = rand(500, 5, 1)
    v3 = rand(500, 4, 1)
    v4 = rand(500, 5, 4)

    c1 = Chains(v1, ["a1", "a2", "a3", "a4", "a5"])
    c2 = Chains(v2, ["a1", "a2", "a3", "a4", "a5"], start = 501)
    c3 = Chains(v3, ["z1", "z2", "z3", "z4"])
    c4 = Chains(v4, ["a1", "a2", "a3", "a4", "a5"])

    # Test dim 1
    c1_2 = cat(c1, c2; dims = 1)
    @test c1_2.value.data == cat(v1, v2, dims=1)
    @test range(c1_2) == 1:1:1000
    @test names(c1_2) == names(c1) == names(c2)
    @test chains(c1_2) == chains(c1) == chains(c2)
    @test c1_2.value == vcat(c1, c2).value

    # Test dim 2
    c1_3 = cat(c1, c3; dims = 2)
    @test c1_3.value.data == cat(v1, v3, dims=2)
    @test range(c1_3) == 1:1:500
    @test names(c1_3) == cat(names(c1), names(c3), dims=1)
    @test chains(c1_3) == chains(c1) == chains(c3)
    @test c1_3.value == hcat(c1, c3).value

    # Test dim 3
    c1_4 = cat(c1, c4; dims = 3)
    @test c1_4.value.data == cat(v1, v4, dims=3)
    @test range(c1_4) == 1:1:500
    @test names(c1_4) == names(c1) == names(c4)
    @test length(chains(c1_4)) == length(chains(c1)) + length(chains(c4))
    @test c1_4.value == chainscat(c1, c4).value
end

@testset "vcat" begin
    chn = Chains(rand(10, 5, 2), ["a", "b", "c", "d", "e"], Dict(:internal => ["d", "e"]))
    chn1 = Chains(rand(5, 5, 2), ["a", "b", "c", "d", "e"], Dict(:internal => ["a", "b"]))

    # incorrect iterations
    @test_throws ArgumentError vcat(chn, Chains(rand(2, 5, 2)))

    # incorrect names
    @test_throws ArgumentError vcat(chn, Chains(rand(10, 5, 2), ["a", "b", "c", "d", "f"]; start=11))

    # incorrect number of chains
    @test_throws ArgumentError vcat(chn, Chains(rand(10, 5, 3), ["a", "b", "c", "d", "e"]; start=11))

    # concate the same chain
    chn_shifted = setrange(chn, 11:20)
    chn2 = vcat(chn, chn_shifted)
    @test chn2.value.data == vcat(chn.value.data, chn.value.data)
    @test size(chn2) == (20, 5, 2)
    @test names(chn2) == names(chn)
    @test range(chn2) == 1:20
    @test chn2.name_map == (parameters = [:a, :b, :c], internal = [:d, :e])

    chn2a = cat(chn, chn_shifted)
    @test chn2a.value == chn2.value
    @test chn2a.name_map == chn2.name_map
    @test chn2a.info == chn2.info

    chn2b = cat(chn, chn_shifted; dims = Val(1))
    @test chn2b.value == chn2.value
    @test chn2b.name_map == chn2.name_map
    @test chn2b.info == chn2.info

    chn2c = cat(chn, chn_shifted; dims = 1)
    @test chn2c.value == chn2.value
    @test chn2c.name_map == chn2.name_map
    @test chn2c.info == chn2.info

    # concatenate a different chain
    chn1_shifted = setrange(chn1, 11:15)
    chn3 = vcat(chn, chn1_shifted)
    @test chn3.value.data == vcat(chn.value.data, chn1.value.data)
    @test size(chn3) == (15, 5, 2)
    @test names(chn3) == names(chn)
    @test range(chn3) == 1:15
    # just take the name map of first argument
    @test chn3.name_map == (parameters = [:a, :b, :c], internal = [:d, :e])

    chn3a = cat(chn, chn1_shifted)
    @test chn3a.value == chn3.value
    @test chn3a.name_map == chn3.name_map
    @test chn3a.info == chn3.info

    chn3b = cat(chn, chn1_shifted; dims = Val(1))
    @test chn3b.value == chn3.value
    @test chn3b.name_map == chn3.name_map
    @test chn3b.info == chn3.info

    chn3c = cat(chn, chn1_shifted; dims = 1)
    @test chn3c.value == chn3.value
    @test chn3c.name_map == chn3.name_map
    @test chn3c.info == chn3.info
end

@testset "hcat" begin
    chn = Chains(rand(10, 3, 2), ["a", "b", "c"], Dict(:internal => ["c"]))
    chn1 = Chains(rand(10, 2, 2), ["d", "e"], Dict(:internal => ["d"]))

    # incorrect ranges
    @test_throws ArgumentError hcat(chn, Chains(rand(10, 2, 2); start = 2))
    @test_throws ArgumentError hcat(chn, Chains(rand(10, 2, 2); thin = 2))

    # non-unique names
    @test_throws ArgumentError hcat(chn, chn)

    # incorrect number of chains
    @test_throws ArgumentError hcat(chn, Chains(rand(10, 2, 3)))

    # concatenate two chains
    chn2 = hcat(chn, chn1)
    @test chn2.value.data == hcat(chn.value.data, chn1.value.data)
    @test size(chn2) == (10, 5, 2)
    @test names(chn2) == vcat(names(chn), names(chn1))
    @test range(chn2) == 1:10
    @test chn2.name_map == (parameters = [:a, :b, :e], internal = [:c, :d])
    
    chn2a = cat(chn, chn1; dims = Val(2))
    @test chn2a.value == chn2.value
    @test chn2a.name_map == chn2.name_map
    @test chn2a.info == chn2.info

    chn2b = cat(chn, chn1; dims = 2)
    @test chn2b.value == chn2.value
    @test chn2b.name_map == chn2.name_map
    @test chn2b.info == chn2.info
end

@testset "chainscat" begin
    chn = Chains(rand(10, 3, 3), ["a", "b", "c"], Dict(:internal => ["c"]))
    chn1 = Chains(rand(10, 3, 2), ["a", "b", "c"], Dict(:internal => ["a"]))

    # incorrect ranges
    @test_throws ArgumentError chainscat(chn, Chains(rand(10, 3, 2); start = 2))
    @test_throws ArgumentError chainscat(chn, Chains(rand(10, 3, 2); thin = 2))

    # incorrect names
    @test_throws ArgumentError chainscat(chn, Chains(rand(10, 3, 2), ["d", "e", "f"]))

    # concatenate two chains
    chn2 = chainscat(chn, chn1)
    @test chn2.value.data == cat(chn.value.data, chn1.value.data; dims = 3)
    @test size(chn2) == (10, 3, 5)
    @test names(chn2) == names(chn)
    @test range(chn2) == 1:10
    # just keep the name map of the first argument
    @test chn2.name_map == (parameters = [:a, :b], internal = [:c])
    
    chn2a = cat(chn, chn1; dims = Val(3))
    @test chn2a.value == chn2.value
    @test chn2a.name_map == chn2.name_map
    @test chn2a.info == chn2.info

    chn2b = cat(chn, chn1; dims = 3)
    @test chn2b.value == chn2.value
    @test chn2b.name_map == chn2.name_map
    @test chn2b.info == chn2.info
end