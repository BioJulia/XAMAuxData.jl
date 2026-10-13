module XAMAuxDataTests

using XAMAuxData: XAMAuxData, SAM, BAM, AuxTag, Hex, Errors, try_auxtag, try_get, is_well_formed
using Test
using MemoryViews: MemoryView
using FormatSpecimens
using StringViews: StringView
using Aqua

INT_TYPE_TO_CHAR = Dict(
    UInt8 => 'C',
    Int8 => 'c',
    UInt16 => 'S',
    Int16 => 's',
    UInt32 => 'I',
    Int32 => 'i'
)

# This method relies on Base internals that are very unlikely to change,
# but we should test it nonetheless
@testset "Unsafe u8" begin
    @test all(XAMAuxData.unsafe_u8(c) == UInt8(c) for c in '\0':'\x7f')
end

@testset "AuxTag" begin
    @test AuxTag("ac") == AuxTag('a', 'c') == AuxTag(0x61, 0x63)
    for good in ["AB", "A1", "a9", "w5"]
        @test try_auxtag(good) isa AuxTag
        @test try_auxtag(Test.GenericString(good)) isa AuxTag
    end
    for bad in ["", "A", "ABC", "A11", "5A", "AØ", "a!"]
        @test try_auxtag(bad) === nothing
        @test try_auxtag(Test.GenericString(bad)) === nothing
    end

    @test sort!(AuxTag.(["XA", "X1", "ab", "AB"])) == AuxTag.(["AB", "X1", "XA", "ab"])

    @testset "AbstractString interface" begin
        t = AuxTag("Ab")
        @test AuxTag(t) === t
        @test String(t) == "Ab"
        @test collect(t) == ['A', 'b']
        @test length(t) == ncodeunits(t) == 2
        @test codeunits(t) == codeunits("Ab")
        @test isvalid(t)
        @test [isvalid(t, i) for i in 0:3] == [false, true, true, false]
        @test_throws BoundsError codeunit(t, 3)
        @test t * "c" == "Abc"
        @test sprint(print, t) == "Ab"
        @test sprint(show, t) == "AuxTag(\"Ab\")"
        @test convert(AuxTag, StringView(b"Ab")) === t
    end

    @testset "isequal and hash" begin
        t = AuxTag("Ab")
        for s in ["Ab", SubString("xAbx", 2, 3), StringView(b"Ab"), Test.GenericString("Ab"), AuxTag('A', 'b')]
            @test t == s
            @test s == t
            @test isequal(t, s)
            @test isequal(s, t)
            @test hash(t) == hash(s)
        end
        for s in ["AB", "ab", "Abc", "", AuxTag("AB")]
            @test !isequal(t, s)
            @test !isequal(s, t)
        end
        @test Dict(t => 1)["Ab"] == 1
        @test Dict("Ab" => 1)[t] == 1
    end
end

@testset "SAM" begin

    @testset "Reading" begin
        str = "AN:A:z\ta1:Z:abc def \tbc:H:A4E9\tkv:i:-25234\tzz:f:-14.466e-3\tAC:B:c,3,2,-34,25,62,-123"
        aux = SAM.Auxiliary(str)
        d = Dict(aux)

        # The array
        @test d == Dict(
            AuxTag("AN") => 'z',
            AuxTag("a1") => "abc def ",
            AuxTag("bc") => Hex(UInt8[0xa4, 0xe9]),
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

    @testset "Display" begin
        aux = SAM.Auxiliary("AB:i:242\tkV:Z:abcdef gh \tJJ:A:!\tzZ:f:-1.2")
        buf = IOBuffer()
        show(buf, MIME"text/plain"(), aux)
        str = String(take!(buf))
        @test str == """
            4-element XAMAuxData.SAM.Auxiliary{MemoryViews.ImmutableMemoryView{UInt8}}:
              "AB" => 242
              "kV" => "abcdef gh "
              "JJ" => '!'
              "zZ" => -1.2f0"""
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
        @test aux["kM"] == Hex([0x7a, 0xf1, 0x38])
    end

    @testset "Writing" begin
        s = "AN:A:z\ta1:Z:abc def \tbc:H:A4E9\tkv:i:-25234\tzz:f:-14.466e-3\tAC:B:c,3,2,-34,25,62,-123"
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
                    typemin(Int32),
                    typemax(UInt32),
                ]
                empty!(aux)
                aux["AB"] = i
                @test String(MemoryView(aux)) == "AB:i:" * string(i)
                @test aux["AB"] === Int64(i)
            end

            # Only values in [-2^31, 2^32) can be written
            for i in Any[Int64(typemin(Int32)) - 1, Int64(typemax(UInt32)) + 1, typemax(UInt64)]
                @test_throws XAMAuxData.AuxException aux["AB"] = i
            end

            # Only bit integers can be written
            for i in Any[true, BigInt(5)]
                @test_throws ArgumentError aux["AB"] = i
            end

            for bad_data in [
                    "--1",
                    "1e9",
                    "0xaedc",
                    "\tKL:A:p",
                    "111122223333444455556666777788889999",
                    "",
                    "-",
                    "+",
                    "+-1",
                    " 1",
                    "1 ",
                    "- 1",
                    "1 2",
                    "\u00a01",
                    "1_000",
                    "0b1",
                    "1.0",
                    "9223372036854775808",
                    "-9223372036854775809",
                    "18446744073709551616",
                ]
                aux = SAM.Auxiliary("AB:i:" * bad_data)
                @test aux["AB"] == Errors.InvalidInt
                @test !isvalid(aux)
            end

            for (good_data, x) in [
                    "0" => 0,
                    "-0" => 0,
                    "+0" => 0,
                    "+5" => 5,
                    "007" => 7,
                    "-123" => -123,
                    "9223372036854775807" => typemax(Int64),
                    "-9223372036854775808" => typemin(Int64),
                ]
                aux = SAM.Auxiliary("AB:i:" * good_data)
                @test aux["AB"] === x
                @test isvalid(aux)
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
                    " 1.0",
                    "1.0 ",
                    "- 1.0",
                    "nan",
                    "NaN",
                    "inf",
                    "-Inf",
                    "Infinity",
                    "0x1p3",
                    "1.",
                    ".",
                    "",
                    "-",
                    "1e",
                    "1e+",
                    "e5",
                    "1.0e5.0",
                    "1e300",
                ]
                aux = SAM.Auxiliary("AB:f:" * bad_data)
                @test aux["AB"] == Errors.InvalidFloat
                @test !isvalid(aux)
            end

            for (good_data, x) in [
                    "1" => 1.0f0,
                    "-1.5" => -1.5f0,
                    "+.5" => 0.5f0,
                    "12.25e2" => 1225.0f0,
                    "1E-2" => 0.01f0,
                    "-3e+1" => -30.0f0,
                ]
                aux = SAM.Auxiliary("AB:f:" * good_data)
                @test aux["AB"] === x
                @test isvalid(aux)
            end

            for x in [NaN, -Inf, Inf32, NaN16, 1.0e300, big"1e300"]
                aux = SAM.Auxiliary(UInt8[], 1)
                @test_throws XAMAuxData.AuxException aux["AB"] = x
                @test isempty(aux)
            end
        end

        @testset "Strings" begin
            aux = SAM.Auxiliary(rand(UInt8, 4), 5)
            for (n, i) in enumerate(
                    Any[
                        "some content",
                        view("another content", 2:10),
                        Test.GenericString("content"),
                        StringView(collect(codeunits("lkwjdlkd"))),
                    ]
                )
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
        @test hex isa Hex{Memory{UInt8}}
        @test isempty(hex)

        # Writing zero-length strings and hex
        delete!(aux, "AZ")
        delete!(aux, "AB")
        aux["AB"] = ""
        aux["AZ"] = Hex(UInt8[])
        str = aux["AB"]
        @test str isa StringView
        @test isempty(str)
        @test aux["AZ"] == Hex(UInt8[])
        @test String(MemoryView(aux)) == "KL:i:-5\tKM:f:-19.1\tAB:Z:\tAZ:H:"

        # Same in BAM
        aux = BAM.Auxiliary(UInt8[], 1)
        aux["AZ"] = Hex(UInt8[])
        @test aux["AZ"] == Hex(UInt8[])
        @test String(MemoryView(aux)) == "AZH\0"
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

    @testset "isvalid" begin
        s = "ab:i:4\tka:A:m\ttp:f:0.2123\tmn:H:13AF\tAC:B:c,3,2,-34,25,62,-123"
        aux = SAM.Auxiliary(collect(codeunits(s)), 1)
        @test isvalid(aux)
    end

end # SAM

@testset "BAM" begin

    @testset "Reading" begin
        str = "ANAza1Zabc def \0bcHA4E9\0kvin\x9d\xff\xffzzf\xcd\x02m\xbcACBc\6\0\0\0\x03\x02\xde\x19>\x85"
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
            AuxTag("bc") => Hex(UInt8[0xa4, 0xe9]),
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

    @testset "Display" begin
        aux = BAM.Auxiliary("ABS\0\xf2kVZabcdef gh \0JJA!zZf\x9a\x99\x99\xbf")
        buf = IOBuffer()
        show(buf, MIME"text/plain"(), aux)
        str = String(take!(buf))
        @test str == """
            4-element XAMAuxData.BAM.Auxiliary{MemoryViews.ImmutableMemoryView{UInt8}}:
              "AB" => 0xf200
              "kV" => "abcdef gh "
              "JJ" => '!'
              "zZ" => -1.2f0"""
    end

    @testset "Mutating" begin
        aux = BAM.Auxiliary(UInt8[], 1)
        d = Dict{AuxTag, Any}()
        for (k, v) in ["L1" => Int32(15), "SO" => Int8[-5, 15, 18], "Xa" => "some string!", "f1" => Float32(-191.0001)]
            d[k] = v
            aux[k] = v
        end
        @test Dict(aux) == d

        # Mutating existing dicts
        for (k, v) in ["kk" => Float32(-14.24), "SO" => "abc ", "L1" => Int32(-11223344), "Xa" => UInt32[3, typemax(UInt32), 0]]
            aux[k] = v
            d[k] = v
        end
        @test Dict(aux) == d

        # Hex
        aux["kM"] = Hex([0x7a, 0xf1, 0x38])
        @test aux["kM"] == Hex([0x7a, 0xf1, 0x38])

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
                    Int32(500203030),
                ]
                aux["Xx"] = n
                @test aux["Xx"] === n
                buf = IOBuffer()
                write(buf, "Xx", INT_TYPE_TO_CHAR[typeof(n)], n)
                @test String(MemoryView(aux)) == String(take!(buf))
            end
            empty!(aux)

            # Only integer types natively supported by BAM can be written
            for n in Any[Int64(5), UInt64(5), Int128(5), BigInt(5), true]
                @test_throws ArgumentError aux["KX"] = n
            end
            @test isempty(aux)
        end

        @testset "Floats and reals" begin
            aux = BAM.Auxiliary(UInt8[], 1)
            for i in [
                    Float16(-2843.23),
                    Float32(323.235),
                    Float64(123.5),
                    pi,
                    MathConstants.e,
                ]
                x = Float32(i)
                aux["KA"] = i
                @test aux["KA"] === x
                buf = IOBuffer()
                write(buf, "KAf", x)
                @test String(take!(buf)) == String(MemoryView(aux))
            end

            for x in [NaN, -Inf, Inf32, NaN16, 1.0e300]
                aux = BAM.Auxiliary(UInt8[], 1)
                @test_throws XAMAuxData.AuxException aux["AB"] = x
                @test isempty(aux)
            end

            for x in [NaN32, Inf32, -Inf32]
                buf = IOBuffer()
                write(buf, "KAf", x)
                aux = BAM.Auxiliary(take!(buf))
                @test aux["KA"] == Errors.InvalidFloat
                @test !isvalid(aux)
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
            for (n, i) in enumerate(
                    Any[
                        "some content",
                        view("another content", 2:10),
                        Test.GenericString("content"),
                        StringView(collect(codeunits("lkwjdlkd"))),
                        "",
                    ]
                )
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
            @test aux["LA"] == Hex([0x01, 0x02, 0x03, 0xaf])
            @test String(MemoryView(aux)) == "LAH010203AF\0"

            aux = SAM.Auxiliary(b"JK:H:ae1bM9")
            @test aux["JK"] == Errors.InvalidHex

            # Lowercase hex is not permitted by the SAM specs
            for bad in ["ab", "0a", "Ab", "0G", "G0", " 0", "0 ", "ABC", "/0", ":0", "@0"]
                aux = SAM.Auxiliary("JK:H:" * bad)
                @test aux["JK"] == Errors.InvalidHex
                @test !isvalid(aux)
                aux = BAM.Auxiliary("JKH" * bad * "\0")
                @test aux["JK"] == Errors.InvalidHex
                @test !isvalid(aux)
            end

            aux = SAM.Auxiliary("JK:H:0123456789ABCDEF")
            @test aux["JK"] == Hex([0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef])
            @test isvalid(aux)

            # Test all bytes, and a length where the decoding loop SIMDs.
            # Tabs are skipped, since they separate SAM fields.
            for byte in 0x00:0xff
                byte == UInt8('\t') && continue
                expected = byte in codeunits("0123456789ABCDEF") ? Hex(fill(UInt8(parse(Int, Char(byte); base = 16) * 17), 32)) : Errors.InvalidHex
                aux = SAM.Auxiliary([codeunits("JK:H:"); fill(byte, 64)])
                @test aux["JK"] == expected
                @test isvalid(aux) == !(expected isa Errors.Error)
            end
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
            d2 = Dict(aux)
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

    @testset "isvalid" begin
        s = "abAcbbHAFDE01\0kaZabcdefgh\0kkC8mnBC\3\0\0\0abc"
        aux = BAM.Auxiliary(collect(codeunits(s)), 1)
        @test isvalid(aux)
    end

end # BAM

@testset "Arrays" begin
    for T in [SAM.Auxiliary, BAM.Auxiliary]
        aux = T(UInt8[], 1)
        for v in Any[
                [0x01, 0x03, 0x05],
                Int8[4, 1, -2],
                UInt16[5, 2, 9, 12, 45221],
                Int16[-2000, 4020, 2],
                UInt32[4, 134, 9],
                Int32[-324932, 2313, 2],
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

        # Unsupported integer element types
        aux = T(rand(UInt8, 10), 11)
        for v in Any[BigInt[-234432212, 2441], Int64[1, 2], UInt64[1], Bool[true], Union{Int8, Int16}[Int8(1)]]
            @test_throws ArgumentError aux["P9"] = v
        end
        @test isempty(aux)

        # Non-Float32 reals are converted to Float32
        mem = Memory{Irrational}(undef, 2)
        copyto!(mem, Irrational[pi, MathConstants.e])
        aux["Kk"] = view(mem, 1:lastindex(mem))
        @test aux["Kk"] isa AbstractVector{Float32}
        @test aux["Kk"] == Float32.(mem)
    end

    @testset "Arrays of non-basic types" begin
        for T in [SAM.Auxiliary, BAM.Auxiliary]
            aux = T(rand(UInt8, 10), 11)
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
                "b,0,1",
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
                "C,0x10",
                "f,1.0,nan",
                "f,inf",
                "f, 1.0",
                "f,0x1p3",
                "f,1e300",
                "C, 1",
                "C,1 ",
                "i,1,,2",
                "c,128",
                "c,-129",
                "I,-1",
                "i,0x10",
            ]
            aux = SAM.Auxiliary("AB:B:" * bad_array)
            @test aux["AB"] == Errors.InvalidArray
            @test !isvalid(aux)
        end

        for (good_array, expected) in Any[
                "C,+5,-0,255" => UInt8[5, 0, 255],
                "c,-128,127" => Int8[-128, 127],
                "I,4294967295" => UInt32[typemax(UInt32)],
                "i,-2147483648" => Int32[typemin(Int32)],
                "f,-1.5,+.5,1E2" => Float32[-1.5, 0.5, 100],
            ]
            aux = SAM.Auxiliary("AB:B:" * good_array)
            @test aux["AB"] == expected
            @test eltype(aux["AB"]) == eltype(expected)
            @test isvalid(aux)
        end

        for x in [NaN32, Inf32]
            buf = IOBuffer()
            write(buf, "ABBf", UInt32(2), 1.0f0, x)
            aux = BAM.Auxiliary(take!(buf))
            @test aux["AB"] == Errors.InvalidArray
            @test !isvalid(aux)
        end

        for T in [SAM.Auxiliary, BAM.Auxiliary]
            for v in [[1.0, NaN], Float32[Inf32], [1.0e300], Real[1, Inf]]
                aux = T(UInt8[], 1)
                @test_throws XAMAuxData.AuxException aux["AB"] = v
                @test isempty(aux)
            end
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
        aux["k1"] = Int32(1234)
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

@testset "try_get" begin
    # SAM
    aux = SAM.Auxiliary("KM:i:252\tAB:Z:αβγδ\tAK::C\tXY:i:1")
    @test try_get(aux, "KM") == 252
    @test try_get(aux, AuxTag("KM")) == 252
    @test aux[AuxTag("KM")] == 252
    @test try_get(aux, "AB") == Errors.InvalidString
    @test try_get(aux, "XY") == Errors.NoColons
    @test try_get(aux, "QQ") == Errors.NoColons
    @test try_get(SAM.Auxiliary("KM:i:252"), "QQ") === nothing
    @test try_get(SAM.Auxiliary(""), "QQ") === nothing
    @test try_get(SAM.Auxiliary("\t"), "QQ") == Errors.TooShortMemory
    @test try_get(SAM.Auxiliary("11:i:1"), "QQ") == Errors.InvalidAuxTag

    # BAM
    aux = BAM.Auxiliary("KMC\x05ABAWXYW123")
    @test try_get(aux, "KM") === 0x05
    @test try_get(aux, "AB") === 'W'
    @test try_get(aux, "XY") == Errors.InvalidTypeTag
    @test try_get(aux, "QQ") == Errors.InvalidTypeTag
    @test try_get(BAM.Auxiliary("KMC\x05"), "QQ") === nothing
    @test try_get(BAM.Auxiliary("ABA\n"), "AB") == Errors.InvalidChar
    @test try_get(BAM.Auxiliary("ABZabc"), "AB") == Errors.NoNullByte
    @test try_get(BAM.Auxiliary("ABBW1234"), "AB") == Errors.InvalidArrayEltype
    @test try_get(BAM.Auxiliary("AB"), "AB") == Errors.TooShortMemory

    # Bad keys still throw
    @test_throws Exception try_get(aux, "1A")
end

@testset "BAM setindex! with start offset" begin
    prefix = collect(b"prefix bytes")
    for (old, new) in Any[
            ("a long string" => "short"),
            ("short" => "a much longer string"),
            (Int32(5) => 0x01),
            (0x01 => Int32(-5)),
            (Int16[1, 2, 3] => UInt8[1]),
        ]
        v = copy(prefix)
        aux = BAM.Auxiliary(v, length(v) + 1)
        aux["AB"] = 'x'
        aux["CD"] = old
        aux["EF"] = Int16(-3)
        aux["CD"] = new
        @test v[1:length(prefix)] == prefix
        @test isvalid(aux)
        @test Dict(aux) == Dict(AuxTag("AB") => 'x', AuxTag("CD") => new, AuxTag("EF") => Int16(-3))
    end
end

@testset "Iterating past malformed data terminates" begin
    # Previously, SAM iteration looped forever after a malformed field
    for str in ["AB:i:1\t", "AB:i:1\tbad\tCD:i:2", "\t", "bad"]
        aux = SAM.Auxiliary(str)
        @test collect(XAMAuxData.iter_encodings(aux))[end] isa XAMAuxData.Error
        @test_throws XAMAuxData.AuxException haskey(aux, "QQ")
    end
end

@testset "length and haskey throw on malformed data" begin
    for aux in [SAM.Auxiliary("AB:i:1\tbad\tEF:i:2"), BAM.Auxiliary("ABAxCDq\x01EFAy")]
        @test !isempty(aux)
        @test_throws XAMAuxData.AuxException length(aux)
        # Keys before the malformed field are found, like with `get`
        @test haskey(aux, "AB")
        @test_throws XAMAuxData.AuxException haskey(aux, "EF")
        @test_throws XAMAuxData.AuxException haskey(aux, "QQ")
    end

    # `isempty` only checks for zero bytes, so it never throws
    aux = SAM.Auxiliary("\t")
    @test !isempty(aux)
    @test_throws XAMAuxData.AuxException length(aux)
end

@testset "BAM iteration with invalid values yields pairs" begin
    aux = BAM.Auxiliary("ABA\nCDC\x05")
    @test collect(aux) == [AuxTag("AB") => Errors.InvalidChar, AuxTag("CD") => 0x05]
    @test Dict(aux) == Dict(AuxTag("AB") => Errors.InvalidChar, AuxTag("CD") => 0x05)
    @test sprint(show, MIME"text/plain"(), aux) isa String
end

@testset "SAM setindex! keeps old value on invalid input" begin
    aux = SAM.Auxiliary(collect(codeunits("AB:i:1\tCD:Z:keep")), 1)
    @test_throws XAMAuxData.AuxException aux["CD"] = "bad\tstring"
    @test_throws XAMAuxData.AuxException aux["AB"] = typemax(Int64)
    @test Dict(aux) == Dict(AuxTag("AB") => 1, AuxTag("CD") => "keep")
end

@testset "setindex_nonexisting! has methods for SAM and BAM" begin
    for T in (SAM.Auxiliary, BAM.Auxiliary)
        aux = T(UInt8[], 1)
        @test XAMAuxData.setindex_nonexisting!(aux, 'k', "BA") === aux
        @test aux["BA"] == 'k'
    end
end

@testset "Mutating after vector is shrunk below start errors" begin
    for T in (SAM.Auxiliary, BAM.Auxiliary)
        v = collect(b"prefix bytes")
        aux = T(v, length(v) + 1)
        aux["AB"] = 'x'
        resize!(v, 3)
        @test isempty(aux)
        @test_throws ErrorException empty!(aux)
        @test_throws ErrorException delete!(aux, "AB")
        @test_throws ErrorException aux["AB"] = 'y'
        @test_throws ErrorException XAMAuxData.setindex_nonexisting!(aux, 'y', "AB")
        @test_throws ErrorException copy!(aux, (T === SAM.Auxiliary ? BAM : SAM).Auxiliary(UInt8[], 1))
        @test v == b"pre"

        # Shrinking to exactly `start - 1` elements is fine
        v = collect(b"prefix")
        aux = T(v, length(v) + 1)
        aux["AB"] = 'x'
        resize!(v, 6)
        aux["CD"] = 'y'
        @test collect(aux) == [AuxTag("CD") => 'y']
    end
end

@testset "Writing values aliasing the auxiliary" begin
    # Whether aliasing values are corrupted when mutating depends on how
    # Base moves data in the vector, so we test both that aliasing is detected,
    # and some cases which are known to have corrupted data before.
    @testset "Detection" begin
        v = collect(codeunits("prefix"))
        aux = BAM.Auxiliary(v, length(v) + 1)
        aux["AB"] = "abc"
        aux["CD"] = Int16[1, -2, 3]
        aux["EF"] = UInt8[4, 5, 6]
        values = Any[
            aux["AB"],
            SubString(aux["AB"], 2, 3),
            aux["CD"],
            view(aux["CD"], 2:3),
            aux["EF"],
            StringView(MemoryView(v)[1:6]),
            v,
        ]
        for x in values
            y = XAMAuxData.unalias_value(aux.x, x)
            @test y !== x
            @test y == x
        end
        h = Hex(aux["EF"])
        y = XAMAuxData.unalias_value(aux.x, h)
        @test parent(y) !== parent(h)
        @test y isa Hex && y == h

        aux = SAM.Auxiliary(collect(codeunits("AB:Z:abc")), 1)
        x = aux["AB"]
        @test XAMAuxData.unalias_value(aux.x, x) !== x

        # Values that do not alias are not copied
        for x in Any["abc", Int16[1, 2], UInt8[], 'a', 0x01, Hex([0x01])]
            @test XAMAuxData.unalias_value(aux.x, x) === x
        end
    end

    @testset "Values are not corrupted" begin
        prefix = join(["P$i:i:$i" for i in 1:9], '\t')
        # SAM: deleting AB moves CD, which the value of CD is a view into
        aux = SAM.Auxiliary(collect(codeunits(prefix * "\tAB:Z:x\tCD:Z:hello world")), 1)
        aux["AB"] = aux["CD"]
        @test aux["AB"] == aux["CD"] == "hello world"

        # BAM: Make sure growing does not reallocate, such that data is moved in-place
        aux = BAM.Auxiliary(sizehint!(UInt8[], 1000), 1)
        for i in 1:9
            aux["P$i"] = Int32(i)
        end
        aux["AB"] = "x"
        aux["CD"] = "hello world"
        aux["EF"] = Int16[1, -2, 3]
        aux["GH"] = UInt8[4, 5, 6, 7]
        aux["AB"] = aux["CD"] # grow
        @test aux["AB"] == aux["CD"] == "hello world"
        aux["AB"] = aux["EF"] # shrink
        @test aux["AB"] == aux["EF"] == [1, -2, 3]
        aux["CD"] = aux["GH"] # shrink
        @test aux["CD"] == aux["GH"] == [4, 5, 6, 7]
        aux["P1"] = aux["CD"] # grow
        @test aux["P1"] == [4, 5, 6, 7]

        # Values may also alias the data before the start of the Auxiliary,
        # which may be moved when the auxiliary is mutated
        for T in (SAM.Auxiliary, BAM.Auxiliary)
            v = collect(codeunits("a fairly long read name"))
            n = length(v)
            aux = T(v, n + 1)
            aux["AB"] = Int32(1)
            for i in Int32(1):Int32(20)
                aux["X" * string(i % 10)] = i
            end
            aux["AB"] = StringView(MemoryView(v)[1:n])
            @test aux["AB"] == "a fairly long read name"
            @test String(v[1:n]) == "a fairly long read name"
        end
    end
end

@testset "Duplicate keys are malformed" begin
    letters = [UInt8(c) for c in ['A':'Z'; 'a':'z']]
    all_tags = [AuxTag(a, b) for a in letters for b in [letters; [UInt8(c) for c in '0':'9']]]
    function make(T, tags)
        aux = T(UInt8[], 1)
        for tag in tags
            XAMAuxData.setindex_nonexisting!(aux, 'x', tag)
        end
        return aux
    end
    for T in (SAM.Auxiliary, BAM.Auxiliary)
        @test is_well_formed(make(T, all_tags))
        @test isvalid(make(T, all_tags))
        # Duplicates both among the first keys, and later
        for n in [2, 5, 16, 17, 40, length(all_tags)], (i, j) in [(1, 2), (1, n), (n - 1, n)]
            tags = all_tags[1:n]
            tags[j] = tags[i]
            aux = make(T, tags)
            @test !is_well_formed(aux)
            @test !isvalid(aux)
            @test startswith(sprint(show, MIME"text/plain"(), aux), "Malformed")
        end
    end
    @test !is_well_formed(SAM.Auxiliary("AB:i:1\tCD:Z:abc\tAB:i:1"))
    @test !is_well_formed(BAM.Auxiliary("ABC\x01CDC\x02ABC\x01"))

    # Opting out of checking for duplicates
    for aux in [SAM.Auxiliary("AB:i:1\tCD:Z:abc\tAB:i:1"), BAM.Auxiliary("ABC\x01CDC\x02ABC\x01")]
        @test is_well_formed(aux; check_duplicates = false)
        @test isvalid(aux; check_duplicates = false)
        @test !is_well_formed(aux; check_duplicates = true)
        @test !isvalid(aux; check_duplicates = true)
    end
    # Other malformations and invalid values are still detected
    for aux in [SAM.Auxiliary("AB:i:1\tAB:i:1\tX"), BAM.Auxiliary("ABC\x01ABC\x01XXq")]
        @test !is_well_formed(aux; check_duplicates = false)
        @test !isvalid(aux; check_duplicates = false)
    end
    for aux in [SAM.Auxiliary("AB:i:1\tAB:i:x"), BAM.Auxiliary("ABC\x01ABA\n")]
        @test is_well_formed(aux; check_duplicates = false)
        @test !isvalid(aux; check_duplicates = false)
    end

    # Checking for duplicates does not allocate.
    # We only optimize codegen for Julia 1.13 and later.
    if VERSION >= v"1.13"
        for T in (SAM.Auxiliary, BAM.Auxiliary)
            aux = make(T, all_tags[1:20])
            is_well_formed(aux)
            isvalid(aux)
            for check_duplicates in (true, false)
                @test @allocations(is_well_formed(aux; check_duplicates)) == 0
                @test @allocations(isvalid(aux; check_duplicates)) == 0
            end
        end
    end
end

@testset "Exception types" begin
    @test_throws(
        "Not enough bytes remaining in Aux data to encode a full field",
        first(values(SAM.Auxiliary("\t"))),
    )

    @test_throws(
        "Invalid AuxTag. Tags must conform to r\"^[A-Za-z][A-Za-z0-9]\$\".",
        first(values(SAM.Auxiliary("11:i:1"))),
    )

    @test_throws(
        "Invalid SAM tag header. Expected <AuxTag>:<type tag>:, but found no colons.",
        first(values(SAM.Auxiliary("ABi111"))),
    )

    @test_throws(
        "BAM string or Hex type not terminated by null byte",
        BAM.Auxiliary("ABZabc")["AB"],
    )

    @test_throws(
        "BAM string or Hex type not terminated by null byte",
        BAM.Auxiliary("ABH1234")["AB"],
    )

    @test_throws(
        "Unknown type tag in aux value.",
        BAM.Auxiliary("ABW123")["AB"],
    )

    # Setting a value into a malformed aux throws the same exception in SAM and BAM
    for (T, s) in [(SAM.Auxiliary, "AB:i:1\t"), (BAM.Auxiliary, "ABC\x01C")]
        aux = T(collect(codeunits(s)), 1)
        @test_throws XAMAuxData.AuxException(Errors.TooShortMemory) aux["CD"] = 'a'
    end

    @test_throws(
        "Invalid array element type tag. Valid values are CcSsIif.",
        BAM.Auxiliary("ABBW1234")["AB"],
    )

    @test_throws(
        "Auxiliary Char (type 'A') must be in '!':'~'.",
        SAM.Auxiliary(UInt8[], 1)["AB"] = '\n',
    )

    @test_throws(
        "SAM integer cannot be parsed as base-ten Int64, or is not in [-2^31, 2^32) when writing or converting to BAM.",
        SAM.Auxiliary(UInt8[], 1)["AB"] = typemax(Int64)
    )

    @test_throws(
        "Float is NaN, infinite or not representable as a Float32, or is not a decimal number in SAM.",
        SAM.Auxiliary(UInt8[], 1)["AB"] = 22.0e304,
    )

    @test_throws(
        "Auxiliary String (type 'Z') can only contain bytes in re\"[ !-~]\".",
        SAM.Auxiliary(UInt8[], 1)["AB"] = "Rødgrød med fløde"
    )
end

@testset "Conversion from SAM to BAM" begin
    kvs = Any[
        AuxTag("VN") => 'i',
        AuxTag("ab") => [0x02, 0x7a, 0xf1],
        AuxTag("F1") => 1.5f0,
        AuxTag("Ni") => Int32(55),
    ]
    for (T1, T2) in [
            (BAM.Auxiliary, SAM.Auxiliary),
            (SAM.Auxiliary, BAM.Auxiliary),
        ]
        dst = T1(UInt8[], 1)
        src = T2(UInt8[], 1)
        dst["KA"] = Int32(55)
        for (k, v) in kvs
            src[k] = v
        end
        copy!(dst, src)

        @test dst == Dict(kvs)
    end

    # SAM integers are converted to the BAM type Int32 or UInt32
    dst = BAM.Auxiliary(UInt8[], 1)
    copy!(dst, SAM.Auxiliary("AB:i:-2147483648\tCD:i:2147483648\tEF:i:4294967295"))
    @test dst["AB"] === typemin(Int32)
    @test dst["CD"] === UInt32(2147483648)
    @test dst["EF"] === typemax(UInt32)
    for bad in ["AB:i:-2147483649", "AB:i:4294967296"]
        @test_throws XAMAuxData.AuxException copy!(BAM.Auxiliary(UInt8[], 1), SAM.Auxiliary(bad))
    end

    # BAM UInt32 values can be converted to SAM
    src = BAM.Auxiliary(UInt8[], 1)
    src["AB"] = typemax(UInt32)
    dst = copy!(SAM.Auxiliary(UInt8[], 1), src)
    @test dst["AB"] === Int64(typemax(UInt32))

    # Existing fields are replaced, and data before the start is kept
    for (T1, T2) in [(SAM.Auxiliary, BAM.Auxiliary), (BAM.Auxiliary, SAM.Auxiliary)]
        src = T2(UInt8[], 1)
        src["AB"] = Int32(1)
        src["CD"] = "abc"
        for old_fields in [[], ["XX" => 'x'], ["XX" => 'x', "AB" => 'y']]
            v = collect(codeunits("prefix"))
            dst = T1(v, 7)
            for (k, x) in old_fields
                dst[k] = x
            end
            copy!(dst, src)
            @test String(v[1:6]) == "prefix"
            @test dst == Dict(AuxTag("AB") => 1, AuxTag("CD") => "abc")
            @test isvalid(dst)
            copy!(dst, T2(UInt8[], 1))
            @test isempty(dst)
            @test v == codeunits("prefix")
        end
    end

    # If copying fails, the destination is unchanged
    bad_sources = [
        (SAM.Auxiliary("AB:i:1\tCD:i:x"), Errors.InvalidInt), # invalid value
        (SAM.Auxiliary("AB:i:1\tCD:Z:abc\tbad"), Errors.TooShortMemory), # malformed
        (SAM.Auxiliary("AB:i:1\tCD:i:4294967296"), Errors.InvalidInt), # out of range
        (BAM.Auxiliary("ABC\x01CDA\n"), Errors.InvalidChar), # invalid value
        (BAM.Auxiliary("ABC\x01CDq\x01"), Errors.InvalidTypeTag), # malformed
    ]
    for (src, err) in bad_sources
        T = src isa SAM.Auxiliary ? BAM.Auxiliary : SAM.Auxiliary
        v = collect(codeunits("prefix"))
        dst = T(v, 7)
        dst["XX"] = 'x'
        dst["YY"] = Int32(5)
        old = copy(v)
        @test try
            copy!(dst, src)
            false
        catch e
            e isa XAMAuxData.AuxException && e.error == err
        end
        @test v == old
    end
end

Aqua.test_all(XAMAuxData; ambiguities = false)

@testset "H values are loaded as Hex" begin
    @test Hex([0x01, 0x02]) == Hex(Memory{UInt8}([0x01, 0x02]))
    @test isequal(Hex([0x01, 0x02]), Hex(Memory{UInt8}([0x01, 0x02])))
    @test hash(Hex([0x01, 0x02])) == hash(Hex(Memory{UInt8}([0x01, 0x02])))
    @test Hex([0x01, 0x02]) != Hex([0x01, 0x03])

    # Hex is an AbstractVector with the same content as its parent
    v = [0x01, 0x02, 0x03]
    h = Hex(v)
    @test parent(h) === v
    @test h == v
    @test isequal(h, v)
    @test hash(h) == hash(v)
    @test length(h) == 3
    @test h[2] == 0x02
    @test_throws BoundsError h[4]
    @test collect(h) == v
    @test isempty(Hex(UInt8[]))

    for (T, s) in [(SAM.Auxiliary, "AB:H:01AF\tCD:H:"), (BAM.Auxiliary, "ABH01AF\0CDH\0")]
        v = collect(codeunits(s))
        aux = T(v, 1)
        @test aux["AB"] isa Hex{Memory{UInt8}}
        @test aux["AB"] == Hex([0x01, 0xaf])

        # Writing a loaded value back preserves the H type, also when empty
        aux["AB"] = aux["AB"]
        aux["CD"] = aux["CD"]
        @test v == codeunits(T == SAM.Auxiliary ? "AB:H:01AF\tCD:H:" : "ABH01AF\0CDH\0")

        # Converting between SAM and BAM also preserves the H type
        U = T == SAM.Auxiliary ? BAM.Auxiliary : SAM.Auxiliary
        dst = copy!(U(UInt8[], 1), aux)
        @test dst["AB"] == Hex([0x01, 0xaf])
        @test dst["CD"] isa Hex && isempty(dst["CD"])
    end
end

@testset "Encoded iterator eltype" begin
    for aux in [SAM.Auxiliary("AB:i:1\tCD:Z:abc\tEF"), BAM.Auxiliary("ABC\1CDZabc\0EF")]
        it = XAMAuxData.iter_encodings(aux)
        T = eltype(it)
        @test T == Union{Errors.Error, Tuple{AuxTag, UInt8, UnitRange{Int}}}
        elements = collect(it)
        @test length(elements) == 3
        @test all(x -> x isa T, elements)
        @test last(elements) isa Errors.Error
    end
end

@testset "Auxiliary type bounds" begin
    for A in [SAM.Auxiliary, BAM.Auxiliary]
        @test A{Vector{UInt8}} isa DataType
        @test_throws TypeError A{Vector{Int}}
        @test_throws TypeError A{String}
    end
end

if VERSION >= v"1.14.0-DEV.1695"
    @test isempty(Test.detect_closure_boxes(XAMAuxData))
end

end # module
