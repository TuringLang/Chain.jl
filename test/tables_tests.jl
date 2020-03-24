using Tables
using TableTraits
using IteratorInterfaceExtensions
using DataFrames

@testset "Tables interface tests" begin

    @testset "Chains" begin
        val = rand(1000, 8, 4)
        colnames = ["a", "b", "c", "d", "e", "f", "g", "h"]
        internal_colnames = ["c", "d", "e", "f", "g", "h"]
        chn = Chains(val, colnames, Dict(:internals => internal_colnames))

        @testset "Tables interface" begin
            @test Tables.istable(typeof(chn))
            @test Tables.columnaccess(typeof(chn))
            @test Tables.columns(chn) === chn
            @test Tables.columnnames(chn) ==
                  (:Iteration, :Chain, :a, :b, :c, :d, :e, :f, :g, :h)
            @test Tables.getcolumn(chn, :Iteration) == [1:1000; 1:1000; 1:1000; 1:1000]
            @test Tables.getcolumn(chn, :Chain) ==
                  [fill(1, 1000); fill(2, 1000); fill(3, 1000); fill(4, 1000)]
            @test Tables.getcolumn(chn, :a) == [
                vec(chn[:, :a, 1].value)
                vec(chn[:, :a, 2].value)
                vec(chn[:, :a, 3].value)
                vec(chn[:, :a, 4].value)
            ]
            @test_throws Exception Tables.getcolumn(chn, :j)
            @test Tables.getcolumn(chn, 1) == Tables.getcolumn(chn, :Iteration)
            @test Tables.getcolumn(chn, 2) == Tables.getcolumn(chn, :Chain)
            @test Tables.getcolumn(chn, 3) == Tables.getcolumn(chn, :a)
            @test_throws Exception Tables.getcolumn(chn, :i)
            @test_throws Exception Tables.getcolumn(chn, 11)
            @test Tables.rowaccess(typeof(chn))
            @test Tables.rows(chn) === chn
            @test length(Tables.rowtable(chn)) == 4000
            nt = Tables.rowtable(chn)[1]
            @test nt ==
                  (; (k => Tables.getcolumn(chn, k)[1] for k in Tables.columnnames(chn))...)
            @test nt == collect(Iterators.take(Tables.namedtupleiterator(chn), 1))[1]
            nt = Tables.rowtable(chn)[2]
            @test nt ==
                  (; (k => Tables.getcolumn(chn, k)[2] for k in Tables.columnnames(chn))...)
            @test nt == collect(Iterators.take(Tables.namedtupleiterator(chn), 2))[2]
            @test Tables.schema(chn) isa Tables.Schema
            @test Tables.schema(chn).names ===
                  (:Iteration, :Chain, :a, :b, :c, :d, :e, :f, :g, :h)
            @test Tables.schema(chn).types === (
                Int64,
                Int64,
                Float64,
                Float64,
                Float64,
                Float64,
                Float64,
                Float64,
                Float64,
                Float64,
            )
            @test Tables.matrix(chn[:, :, 1])[:, 3:end] ≈ chn[:, :, 1].value
            @test Tables.matrix(chn[:, :, 2])[:, 3:end] ≈ chn[:, :, 2].value
        end

        @testset "TableTraits interface" begin
            @test IteratorInterfaceExtensions.isiterable(chn)
            @test TableTraits.isiterabletable(chn)
            nt = collect(Iterators.take(IteratorInterfaceExtensions.getiterator(chn), 1))[1]
            @test nt ==
                  (; (k => Tables.getcolumn(chn, k)[1] for k in Tables.columnnames(chn))...)
            nt = collect(Iterators.take(IteratorInterfaceExtensions.getiterator(chn), 2))[2]
            @test nt ==
                  (; (k => Tables.getcolumn(chn, k)[2] for k in Tables.columnnames(chn))...)
        end

        @testset "DataFrames.DataFrame constructor" begin
            @inferred DataFrame(chn)
            @test DataFrame(chn) isa DataFrame
            df = DataFrame(chn)
            @test Tables.columntable(df) == Tables.columntable(chn)
        end
    end

    @testset "ChainDataFrames" begin
        val = rand(1000, 8, 4)
        colnames = ["a", "b", "c", "d", "e", "f", "g", "h"]
        internal_colnames = ["c", "d", "e", "f", "g", "h"]
        chn = Chains(val, colnames, Dict(:internals => internal_colnames))
        cdf = describe(chn)[1]

        @testset "Tables interface" begin
            @test Tables.istable(typeof(cdf))
            @test Tables.columnaccess(typeof(cdf))
            @test Tables.columns(cdf) === cdf
            @test Tables.columnnames(cdf) == keys(cdf.nt)
            @testset for k in keys(cdf.nt)
                @test Tables.getcolumn(cdf, k) == getproperty(cdf.nt, k)
            end
            @test Tables.getcolumn(cdf, 1) == Tables.getcolumn(cdf, keys(cdf.nt)[1])
            @test Tables.getcolumn(cdf, 2) == Tables.getcolumn(cdf, keys(cdf.nt)[2])
            @test_throws Exception Tables.getcolumn(cdf, :blah)
            @test_throws Exception Tables.getcolumn(cdf, length(keys(cdf.nt)) + 1)
            @test Tables.rowaccess(typeof(cdf))
            @test Tables.rows(cdf) === cdf
            @test length(Tables.rowtable(cdf)) == length(values(cdf.nt)[1])
            @test Tables.columntable(cdf) == cdf.nt
            nt = Tables.rowtable(cdf)[1]
            @test nt == (; (k => getproperty(cdf.nt, k)[1] for k in keys(cdf.nt))...)
            @test nt == collect(Iterators.take(Tables.namedtupleiterator(cdf), 1))[1]
            nt = Tables.rowtable(cdf)[2]
            @test nt == (; (k => getproperty(cdf.nt, k)[2] for k in keys(cdf.nt))...)
            @test nt == collect(Iterators.take(Tables.namedtupleiterator(cdf), 2))[2]
            @test Tables.schema(cdf) isa Tables.Schema
            @test Tables.schema(cdf).names === keys(cdf.nt)
            @test Tables.schema(cdf).types === eltype.(values(cdf.nt))
        end

        @testset "TableTraits interface" begin
            @test IteratorInterfaceExtensions.isiterable(cdf)
            @test TableTraits.isiterabletable(cdf)
            nt = collect(Iterators.take(IteratorInterfaceExtensions.getiterator(cdf), 1))[1]
            @test nt == (; (k => getproperty(cdf.nt, k)[1] for k in keys(cdf.nt))...)
            nt = collect(Iterators.take(IteratorInterfaceExtensions.getiterator(cdf), 2))[2]
            @test nt == (; (k => getproperty(cdf.nt, k)[2] for k in keys(cdf.nt))...)
        end

        @testset "DataFrames.DataFrame constructor" begin
            @inferred DataFrame(cdf)
            df = DataFrame(cdf)
            @test df isa DataFrame
            @test Tables.columntable(df) == cdf.nt
        end
    end
end
