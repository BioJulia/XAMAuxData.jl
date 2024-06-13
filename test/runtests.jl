module XAMAuxDataTests

using XAMAuxData: SAM, BAM, AuxTag, Hex
using Test
using MemViews: MemView
using FormatSpecimens

@testset "SAM" begin

@testset "Reading" begin
    str = "AN:A:z\ta1:Z:abc def \tbc:H:a4e9\tkv:i:-25234\tzz:f:-14.466e-3\tAC:B:c3,2,-34,25,62,-123"
    aux = SAM.Auxiliary(str)
    d = Dict(aux)

    # The array
    @test d == Dict(
        AuxTag("AN") => 'z',
        AuxTag("a1") => "abc def ",
        AuxTag("bc") => UInt8[0xa4, 0xe9],
        AuxTag("zz") => -14.466f-3,
        AuxTag("AC") => Int8[3, 2, -34, 25, 62, -123],
        AuxTag("kv") => -25234,
    )

    for (uint_tag, et) in [('C', UInt8), ('S', UInt16), ('I', UInt32)]
        s = "ab:B:$(uint_tag)47,1,252,44,11"
        (k, v) = only(SAM.Auxiliary(s))
        @test k === AuxTag("ab")
        @test v == [47, 1, 252, 44, 11]
        @test eltype(v) == et
    end

    for (int_tag, et) in [('c', Int8), ('s', Int16), ('i', Int32)]
        s = "ab:B:$(int_tag)-13,12,-128,127,15"
        (k, v) = only(SAM.Auxiliary(s))
        @test k === AuxTag("ab")
        @test v == [-13, 12, -128, 127, 15]
        @test eltype(v) == et
    end
end

@testset "Mutating" begin
    aux = SAM.Auxiliary(UInt8[], 1)
    d = Dict{AuxTag, Any}()
    for (k, v) in ["L1" => 15, "SO" => Int8[-5, 15, 18], "Xa" => "some string!", "f1" => Float32(-191.0001)]
        d[k] = v
        aux[k] = v
    end
    @test Dict(aux) == d

    # Mutating existing dicts
    for (k, v) in ["kk" => Float32(-14.24), "SO" => "abc ", "L1" => -11223344, "Xa" => UInt32[3, typemax(UInt32), 0]]
        aux[k] = v
        d[k] = v
    end
    @test Dict(aux) == d

    # Hex
    aux["kM"] = Hex([0x7a, 0xf1, 0x38])
    @test aux["kM"] == [0x7a, 0xf1, 0x38]
end

@testset "Writing" begin
    s = "AN:A:z\ta1:Z:abc def \tbc:H:a4e9\tkv:i:-25234\tzz:f:-14.466e-3\tAC:B:c3,2,-34,25,62,-123"
    cu = codeunits(s)
    mv = MemView(s)
    @test mv == cu
    mem = Memory{UInt8}(undef, length(cu))
    copyto!(mem, 1, mv, 1, length(mv))
    @test mem == cu == mv
end

@testset "Existing good files" begin
    all_is_valid = true
    for file_metadata in list_valid_specimens("SAM")
        path = joinpath(path_of_format("SAM"), file_metadata["filename"])
        lines = collect(eachline(path))
        filter!(!startswith('@'), lines)
        for line in lines
            fields = split(line, '\t')
            auxs = join(fields[12:end], '\t')
            aux = SAM.Auxiliary(auxs)
            all_is_valid &= isvalid(aux)
            if !all_is_valid
                error("Invalid: " * line)
            end
            for (k, v) in aux
                nothing
            end
        end
    end
    @test all_is_valid
end

end

@testset "BAM" begin

@testset "Reading" begin
    str = "ANAza1Zabc def \0bcHa4e9\0kvin\x9d\xff\xffzzf\xcd\x02m\xbcACBc\6\0\0\0\x03\x02\xde\x19>\x85"
    aux = BAM.Auxiliary(str)
    d = Dict(aux)

    # The array
    @test d == Dict(
        AuxTag("AN") => 'z',
        AuxTag("a1") => "abc def ",
        AuxTag("bc") => UInt8[0xa4, 0xe9],
        AuxTag("zz") => -14.466f-3,
        AuxTag("AC") => Int8[3, 2, -34, 25, 62, -123],
        AuxTag("kv") => -25234,
    )

    for (uint_tag, et) in [('C', UInt8), ('S', UInt16), ('I', UInt32)]
        buf = IOBuffer()
        write(buf, "abB", UInt8(uint_tag), Int32(5))
        arr = et.([47, 1, 252, 44, 11])
        write(buf, arr)
        s = String(take!(buf))
        (k, v) = only(BAM.Auxiliary(s))
        @test k === AuxTag("ab")
        @test v == arr
        @test eltype(v) == et
    end

    for (int_tag, et) in [('c', Int8), ('s', Int16), ('i', Int32)]
        buf = IOBuffer()
        write(buf, "abB", UInt8(int_tag), Int32(5))
        arr = et.([-13, 12, -128, 127, 15])
        write(buf, arr)
        s = String(take!(buf))
        (k, v) = only(BAM.Auxiliary(s))
        @test k === AuxTag("ab")
        @test v == arr
        @test eltype(v) == et
    end
end

@testset "Mutating" begin
    aux = BAM.Auxiliary(UInt8[], 1)
    d = Dict{AuxTag, Any}()
    for (k, v) in ["L1" => 15, "SO" => Int8[-5, 15, 18], "Xa" => "some string!", "f1" => Float32(-191.0001)]
        d[k] = v
        aux[k] = v
    end
    @test Dict(aux) == d

    # Mutating existing dicts
    for (k, v) in ["kk" => Float32(-14.24), "SO" => "abc ", "L1" => -11223344, "Xa" => UInt32[3, typemax(UInt32), 0]]
        aux[k] = v
        d[k] = v
    end
    @test Dict(aux) == d

    # Hex
    aux["kM"] = Hex([0x7a, 0xf1, 0x38])
    @test aux["kM"] == [0x7a, 0xf1, 0x38]

    # Correct integer type
    for (k, v) in Any[
        "B4" => 0xa5,
        "B6" => Int8(-9),
        "kN" => 0x41f8,
        "nP" => Int16(-5124),
        "Ww" => 0x01020304,
        "Jj" => Int32(2053352342),
    ]
        aux[k] = v
        @test aux[k] === v
    end
end

end

end # module