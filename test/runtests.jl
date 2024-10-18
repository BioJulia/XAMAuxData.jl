module XAMAuxDataTests

using XAMAuxData: XAMAuxData, SAM, BAM, AuxTag, Hex, DelimitedIterator, Errors, try_auxtag
using Test
using MemoryViews: MemoryView
using FormatSpecimens
using StringViews: StringView

INT_TYPE_TO_CHAR = Dict(
    UInt8 => 'C',
    Int8 => 'c',
    UInt16 => 'S',
    Int16 => 's',
    UInt32 => 'I',
    Int32 => 'i'
)

@testset "DelimitedIterator" begin
    it = DelimitedIterator("abceacdga", UInt8('a'))
    @test collect(it) == [
        b"",
        b"bce",
        b"cdg",
        b"",
    ]
end

# This method relies on Base internals that are very unlikely to change,
# but we should test it nonetheless
@testset "Unsafe u8" begin
    @test all(XAMAuxData.unsafe_u8(c) == UInt8(c) for c in '\0':'\x7f')
end

@testset "AuxTag" begin
    @test AuxTag("ac") == AuxTag('a', 'c') == AuxTag(0x61, 0x63)
    for good in ["AB", "A1", "a9", "w5"]
        @test try_auxtag(good) isa AuxTag
    end
    for bad in ["", "A", "ABC", "A11", "5A", "AØ", "a!"]
        @test try_auxtag(bad) === nothing
    end

    @test sort!(AuxTag.(["XA", "X1", "ab", "AB"])) == AuxTag.(["AB", "X1", "XA", "ab"])
end

@testset "SAM" begin

@testset "Reading" begin
    str = "AN:A:z\ta1:Z:abc def \tbc:H:a4e9\tkv:i:-25234\tzz:f:-14.466e-3\tAC:B:c,3,2,-34,25,62,-123"
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

    @test get(aux, "kv", nothing) == -25234
    @test get(aux, "aN", 0x99) === 0x99 

    for (uint_tag, et) in [('C', UInt8), ('S', UInt16), ('I', UInt32)]
        s = "ab:B:$(uint_tag),47,1,252,44,11"
        (k, v) = only(SAM.Auxiliary(s))
        @test k === AuxTag("ab")
        @test v == [47, 1, 252, 44, 11]
        @test eltype(v) == et
    end

    for (int_tag, et) in [('c', Int8), ('s', Int16), ('i', Int32)]
        s = "ab:B:$(int_tag),-13,12,-128,127,15"
        (k, v) = only(SAM.Auxiliary(s))
        @test k === AuxTag("ab")
        @test v == [-13, 12, -128, 127, 15]
        @test eltype(v) == et
    end

    aux = SAM.Auxiliary("ka:i:-13\thc:w:1234")
    @test aux["hc"] == Errors.InvalidTypeTag
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
    s = "AN:A:z\ta1:Z:abc def \tbc:H:a4e9\tkv:i:-25234\tzz:f:-14.466e-3\tAC:B:c,3,2,-34,25,62,-123"
    cu = codeunits(s)
    mv = MemoryView(s)
    @test mv == cu
    mem = Memory{UInt8}(undef, length(cu))
    copyto!(mem, 1, mv, 1, length(mv))
    @test mem == cu == mv

    # Test encoding of various different types

    @testset "Integers" begin
        aux = SAM.Auxiliary(rand(UInt8, 19), 20)
        for i in Any[
            UInt8(55),
            Int8(12),
            UInt16(22),
            Int16(-99),
            UInt32(3444251),
            Int32(43094234),
            UInt64(21345353),
            Int64(231435511),
            UInt128(34311),
            Int128(-391239931),
            BigInt(2033843841),
        ]
            empty!(aux)
            aux["AB"] = i
            @test String(MemoryView(aux)) == "AB:i:" * string(i)
            @test aux["AB"] === Int(i)
        end

        for bad_data in [
            "--1",
            "1e9",
            "0xaedc",
            "\tKL:A:p",
            "111122223333444455556666777788889999",
        ]
            aux = SAM.Auxiliary("AB:i:" * bad_data)
            @test aux["AB"] == Errors.InvalidInt
        end
    end

    @testset "Floats" begin
        aux = SAM.Auxiliary(rand(UInt8, 19), 20)
        aux["xa"] = "some string!  "
        for i in Any[
            Float16(19.9),
            Float32(-2034.22342),
            Float64(-5.353e-11),
            BigFloat("231.23234e22"),
            pi,
        ]
            aux["FL"] = i
            @test String(MemoryView(aux)) == "xa:Z:some string!  \tFL:f:" * string(Float32(i))
            @test aux["FL"] === Float32(i)
        end

        for bad_data in [
            "Ke300",
            "++++9.0",
            "1.3e2F0",
            "\tKL:A:p",
        ]
            aux = SAM.Auxiliary("AB:f:" * bad_data)
            @test aux["AB"] == Errors.InvalidFloat
        end
    end

    @testset "Strings" begin
        aux = SAM.Auxiliary(rand(UInt8, 4), 5)
        for (n,i) in enumerate(Any[
            "some content",
            view("another content", 2:10),
            Test.GenericString("content"),
            StringView(collect(codeunits("lkwjdlkd"))),
        ])
            aux["S" * string(n)] = i
            @test String(MemoryView(aux)) == "S" * string(n) * ":Z:" * String(i)
            @test aux["S" * string(n)] == i
            delete!(aux, "S" * string(n))
        end

        aux = SAM.Auxiliary(rand(UInt8, 4), 5)
        @test_throws Exception aux["KL"] = "a\tb"

        for bad_data in [
            "abcde\ra",
            "a\0\0a",
            "æøå",
        ]
            aux = SAM.Auxiliary("AB:Z:" * bad_data)
            @test aux["AB"] == Errors.InvalidString
        end
    end

    @testset "Chars" begin
        aux = SAM.Auxiliary(rand(UInt8, 3), 4)
        for c in Any[
            'w',
            '1',
        ]
            aux["k1"] = c
            @test String(MemoryView(aux)) == "k1:A:" * Char(c)
            @test aux["k1"] === Char(c)
            empty!(aux)
        end

        for bad_data in [
            "Æ",
            "\f",
            "\x7f",
            "ab",
            "\tKV:Z:A",
            "ab\tKV:Z:A",
            "",
        ]
            aux = SAM.Auxiliary("AB:A:" * bad_data)
            @test aux["AB"] == Errors.InvalidChar
        end
    end

    @testset "Roundtrip" begin
        d1 = Dict{AuxTag, Any}("AZ" => 15, "pU" => "hello, world!", "xx" => Float32(1.1322))
        aux = SAM.Auxiliary(UInt8[], 1)
        merge!(aux, d1)
        str = String(MemoryView(aux))
        d2 = Dict(SAM.Auxiliary(str))
        @test d1 == d2
    end
end

@testset "Zero-length strings and hex" begin
    # Reading zero-length strings and hex
    aux = SAM.Auxiliary(collect(codeunits("ABCDEFGHAB:Z:\tKL:i:-5\tAZ:H:\tKM:f:-19.1")), 9)
    str = aux["AB"]
    @test str isa StringView
    @test isempty(str)
    hex = aux["AZ"]
    @test hex isa Memory{UInt8}
    @test isempty(hex)

    # Writing: Can't write zero-length hex, but can write zero-length string
    delete!(aux, "AZ")
    delete!(aux, "AB")
    aux["AB"] = ""
    @test_throws Exception aux["AZ"] = Hex(UInt8[])
    str = aux["AB"]
    @test str isa StringView
    @test isempty(str)
    @test_throws Exception aux["AZ"]
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

@testset "Deletion" begin
    s = "PG:Z:slfda\tkj:i:-2443\tHJ:A:p"
    aux = SAM.Auxiliary(collect(codeunits(s)), 1)
    delete!(aux, "kj")
    @test String(MemoryView(aux)) == "PG:Z:slfda\tHJ:A:p"
    delete!(aux, "HJ")
    @test String(MemoryView(aux)) == "PG:Z:slfda"
    delete!(aux, "PG")
    @test isempty(MemoryView(aux))
end

end # SAM

@testset "BAM" begin

@testset "Reading" begin
    str = "ANAza1Zabc def \0bcHa4e9\0kvin\x9d\xff\xffzzf\xcd\x02m\xbcACBc\6\0\0\0\x03\x02\xde\x19>\x85"
    aux = BAM.Auxiliary(str)
    d = Dict(aux)

    @test get(aux, "AN", nothing) === 'z'
    @test get(aux, "zz", 5) === -14.466f-3
    @test get(aux, "AB", 0x98) === 0x98
    @test get(aux, "kN", missing) === missing

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

@testset "Writing" begin
    @testset "Integers" begin
        aux = BAM.Auxiliary(rand(UInt8, 10), 11)
        for n in Any[
            UInt8(62),
            Int8(-112),
            UInt16(500),
            Int16(-31001),
            UInt32(1_000_000_102),
            Int32(500203030)
        ]
            aux["Xx"] = n
            @test aux["Xx"] === n
            buf = IOBuffer()
            write(buf, "Xx", INT_TYPE_TO_CHAR[typeof(n)], n)
            @test String(MemoryView(aux)) == String(take!(buf))
        end
        empty!(aux)
        aux["KX"] = BigInt(340932)
        @test aux["KX"] === Int32(340932)
    end

    @testset "Floats and reals" begin
        aux = BAM.Auxiliary(UInt8[], 1)
        for i in [
            Float16(-2843.23),
            Float32(323.235),
            Float64(NaN),
            Float64(-Inf),
            Float64(123.5),
            pi,
            MathConstants.e
        ]
            x = Float32(i)
            aux["KA"] = i
            @test aux["KA"] === x
            buf = IOBuffer()
            write(buf, "KAf", x)
            @test String(take!(buf)) == String(MemoryView(aux))
        end
    end

    @testset "Chars" begin
        aux = BAM.Auxiliary(rand(UInt8, 3), 4)
        for c in Any[
            'w',
            '1',
        ]
            aux["k1"] = c
            @test aux["k1"] === Char(c)
            @test String(MemoryView(aux)) == "k1A" * Char(c)
            empty!(aux)
        end

        for bad_data in [
            "Æ",
            "\f",
            "\x7f",
        ]
            aux = BAM.Auxiliary("ABA" * bad_data)
            @test aux["AB"] == Errors.InvalidChar
        end
    end

    @testset "Strings" begin
        aux = BAM.Auxiliary(rand(UInt8, 4), 5)
        for (n,i) in enumerate(Any[
            "some content",
            view("another content", 2:10),
            Test.GenericString("content"),
            StringView(collect(codeunits("lkwjdlkd"))),
            ""
        ])
            aux["S" * string(n)] = i
            @test String(MemoryView(aux)) == "S" * string(n) * "Z" * String(i) * '\0'
            @test aux["S" * string(n)] == i
            delete!(aux, "S" * string(n))
        end

        @test BAM.Auxiliary("ABZab\x7fa\0")["AB"] == Errors.InvalidString
        @test BAM.Auxiliary("ABZab\r\0")["AB"] == Errors.InvalidString
        @test BAM.Auxiliary("ABZ\na\0")["AB"] == Errors.InvalidString
    end

    @testset "Hex" begin
        aux = BAM.Auxiliary(UInt8[], 1)
        aux["LA"] = Hex([0x01, 0x02, 0x03, 0xaf])
        @test aux["LA"] == [0x01, 0x02, 0x03, 0xaf]
        @test String(MemoryView(aux)) == "LAH010203AF\0"

        aux = SAM.Auxiliary(b"JK:H:ae1bM9")
        @test aux["JK"] == Errors.InvalidHex
    end

    @testset "Round trip" begin
        d1 = Dict{AuxTag, Any}(
            "k1" => UInt8(5),
            "S9" => Int16(-199),
            "fa" => Float32(919.3221e3),
            "ch" => 'W',
            "st" => "some dD@!~~]]  ",
            "hx" => Int16[1, 3, -24221],
            "lM" => UInt32[2314122, 2423, 92933, 0],
        )
        aux = BAM.Auxiliary(rand(UInt8, 19), 20)
        merge!(aux, d1)
        d2 = Dict(d1)
        @test d1 == d2
    end

    @testset "Deletion" begin
        s = "pLdsryh352PGZslfda\0kji\4\0\0\0HJAp"
        aux = BAM.Auxiliary(collect(codeunits(s)), 11)
        @test String(MemoryView(aux)) == "PGZslfda\0kji\4\0\0\0HJAp"
        delete!(aux, "kj")
        @test String(MemoryView(aux)) == "PGZslfda\0HJAp"
        delete!(aux, "HJ")
        @test String(MemoryView(aux)) == "PGZslfda\0"
        delete!(aux, "PG")
        @test isempty(MemoryView(aux))
    end
end

end # BAM

@testset "Arrays" begin
    for T in [SAM.Auxiliary, BAM.Auxiliary]
        aux = T(UInt8[], 1)
        for v in Any[
            [0x01, 0x03, 0x05],
            Int8[4,1,-2],
            UInt16[5, 2, 9, 12, 45221],
            Int16[-2000, 4020, 2],
            UInt32[4, 134, 9],
            Int32[-324932, 2313, 2]
        ]
            aux["A1"] = view(v, 1:lastindex(v))
            @test aux["A1"] == v
            if T == SAM.Auxiliary
                @test String(MemoryView(aux)) == "A1:B:" * INT_TYPE_TO_CHAR[eltype(v)] * ',' * join(map(string, v), ',')
                @test aux["A1"] isa Memory{eltype(v)}
            elseif T == BAM.Auxiliary
                @test aux["A1"] isa AbstractVector{}
                buf = IOBuffer()
                write(buf, "A1B", INT_TYPE_TO_CHAR[eltype(v)], UInt32(length(v)))
                write(buf, v)
                @test String(MemoryView(aux)) == String(take!(buf))
            else
                error()
            end
            empty!(aux)
        end

        # Unsupported types
        aux = T(rand(UInt8, 10), 11)
        v = BigInt[-234432212, 2441, 2423133]
        aux["P9"] = v
        if T == SAM.Auxiliary
            @test String(MemoryView(aux)) == "P9:B:i,-234432212,2441,2423133"
            @test aux["P9"] isa Memory{Int32}
        elseif T == BAM.Auxiliary
            buf = IOBuffer()
            write(buf, "P9Bi", Int32(length(v)), Int32.(v))
            @test String(MemoryView(aux)) == String(take!(buf))
        else
            error()
        end
        empty!(aux)
        mem = Memory{Irrational}(undef, 2)
        copyto!(mem, Irrational[pi, MathConstants.e])
        aux["Kk"] = view(mem, 1:lastindex(mem))
        @test aux["Kk"] isa AbstractVector{Float32}
        @test aux["Kk"] == Float32.(mem)
    end

    @testset "Arrays of non-basic types" begin
        for T in [SAM.Auxiliary, BAM.Auxiliary]
            aux = T(rand(UInt8, 10), 11)
            inp = BigInt[20, -400, 102]
            aux["AB"] = inp
            outp = aux["AB"]
            @test outp isa (T == SAM.Auxiliary ? Memory{Int32} : AbstractVector{Int32})
            @test outp == inp

            empty!(aux)
            inp = [pi, MathConstants.e]
            aux["KV"] = inp
            outp = aux["KV"]
            @test outp isa (T == SAM.Auxiliary ? Memory{Float32} : AbstractVector{Float32})
            @test outp ≈ inp
        end
    end

    @testset "Bad arrays" begin
        for bad_eltype in [
            "W,1,2,3",
            "F,1.0,2.0",
            "b,0,1"
        ]
            aux = SAM.Auxiliary("AB:B:" * bad_eltype)
            @test aux["AB"] == Errors.InvalidArrayEltype
        end

        for bad_array in [
            "C,1,-2,3",
            "c,1,2,3,-128,128",
            "s,39000,1",
            "C1,2,3",
            "I,1,2,3,",
            "I,1,-2,3",
            "i,",
        ]
            aux = SAM.Auxiliary("AB:B:" * bad_array)
            @test aux["AB"] == Errors.InvalidArray
        end

        for T in [SAM.Auxiliary, BAM.Auxiliary]
            aux = T(UInt8[], 1)
            @test_throws Exception aux["AB"] = Vector{Union{}}(undef, 3)
        end
    end
end

@testset "Various operations" begin
    for T in [SAM.Auxiliary, BAM.Auxiliary]
        aux = T(UInt8[], 1)
        @test Dict(empty(T)) == Dict(aux)
        @test length(aux) == 0
        @test isempty(aux)
        aux["k1"] = 1234
        @test length(aux) == 1
        @test !isempty(aux)
        aux["v9"] = 'w'
        @test length(aux) == 2
        aux["AA"] = Int16[-500, 205, 1]
        @test length(aux) == 3

        for k in ["k1", "v9", "AA"]
            @test haskey(aux, k)
        end
        @test collect(keys(aux)) == AuxTag.(["k1", "v9", "AA"])
        for k in ["k2", "v7", "BB", "c1"]
            @test !haskey(aux, k)
        end

        cp = copy(aux)
        d = Dict(aux)
        @test d == Dict(cp)
        len = 3
        for k in ["k1", "v9", "AA"]
            delete!(aux, k)
            len -= 1
            @test length(aux) == len
        end
        @test Dict(cp) == d
    end
end

end # module
