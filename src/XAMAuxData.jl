module XAMAuxData

using MemViews: MemView, MutableMemView, ImmutableMemView
using StringViews: StringView

struct Unsafe end
const unsafe = Unsafe()

export Hex, AuxTag, SAM, BAM

#=
TODO:
Review use of errors. I.e. they must be used consistently.
Maybe InvalidAuxTag is one error and InvalidAuxHeader is another. These are the only ones that
can be emitted from iterating encodings.

# TODO: What about GFA? Maybe make a new module which uses
functions from SAM, and maybe even forwards some methods.
=#

# These are the numerical types supported by the BAM format.
const AUX_NUMBER_TYPES = Union{Int8, UInt8, Int16, UInt16, Int32, UInt32, Float32}

# TODO: Replace with trait
const ByteString = Union{String, SubString{String}, StringView{ImmutableMemView{UInt8}}}

include("auxtag.jl")

"""
    Hex(v::AbstractVector{UInt8})

Wrapper type around a byte vector. When this type is assigned to an Auxiliary,
it is encoded as a hex value (type tag `H`) as opposed to a raw byte array (
type tag `B`).
"""
struct Hex{V <: AbstractVector{UInt8}}
    x::V

    function Hex(v::AbstractVector{UInt8})
        isempty(v) && error("Hex values cannot be empty according to SAM specs")
        new{typeof(v)}(v)
    end
end

# Must encode to uppercase A-F
hexencode_nibble(u::UInt8)::UInt8 = u < 0x0a ? UInt8('0') + u : UInt8('A') - 0x0a + u

function hexencode!(mem::MutableMemView, hex::Hex)
    @inbounds for (byte_no, byte) in enumerate(hex.x)
        mem[2 * byte_no - 1] = hexencode_nibble(byte >> 4)
        mem[2 * byte_no] = hexencode_nibble(byte & 0x0f)
    end
    mem
end

# These are the type tags use to determine the value of the serialized data
# the input type must already be known to be appropriate for the given format,
# i.e. UInt8's are serialized as `i` in SAM, because they should be converted
# to Int32 first.
get_type_tag(::Type{UInt8}) = UInt8('C')
get_type_tag(::Type{Int8}) = UInt8('c')
get_type_tag(::Type{UInt16}) = UInt8('S')
get_type_tag(::Type{Int16}) = UInt8('s')
get_type_tag(::Type{UInt32}) = UInt8('I')
get_type_tag(::Type{Int32}) = UInt8('i')
get_type_tag(::Type{Float32}) = UInt8('f')
get_type_tag(::Type{Char}) = UInt8('A')
get_type_tag(::Type{<:AbstractString}) = UInt8('Z')
get_type_tag(::Type{<:AbstractVector}) = UInt8('B')
get_type_tag(::Type{<:Hex}) = UInt8('H')

function as_aux_type(T::Type{<:AUX_NUMBER_TYPES})
    T != Union{} ? T : error("Cannot convert Union{} to XAM-compatible type")
end

as_aux_type(::Type{<:Real}) = Float32
as_aux_type(::Type{<:Integer}) = Int32

function as_aux_value end
as_sam_aux_value(x) = as_aux_value(x)
as_bam_aux_value(x) = as_aux_value(x)

as_aux_value(x::Real) = Float32(x)::Float32

function as_aux_value(x::AbstractChar)::Union{Char, Error}
    c = Char(x)::Char
    isascii(c) ? c : error("AUX chars must be in '!':'~'")
end

as_aux_value(x::Hex) = x
as_aux_value(s::AbstractString)::ByteString = as_aux_value(String(s))

function as_aux_value(
    s::Union{String, SubString{String}, StringView{ImmutableMemView{UInt8}}}
)::Union{ByteString, Error}
    mem = ImmutableMemView(codeunits(s))
    if is_printable(mem)
        s
    else
        error("AUX string can only contain chars in [ !-~]")
    end
end

# Returns an AbstractVector{<:AUX_NUMBER_TYPES}
function as_aux_value(v::AbstractVector{<:Real})
    eltype(v) != Union{} && eltype(v) <: AUX_NUMBER_TYPES && return v
    E = as_aux_type(eltype(v))
    mem = Memory{E}(undef, length(v))
    for i in eachindex(v, mem)
        mem[i] = E(v[i])::E
    end
    mem
end

function is_printable(v::AbstractVector{UInt8})
    res = true
    for b in v
        res &= (b == UInt8(' ')) | is_printable_char(b)
    end
    res
end

abstract type AbstractAuxiliary{T} <: AbstractDict{AuxTag, Any} end

function (T::Type{<:AbstractAuxiliary})(itr)
    y = empty(T)
    for (k, v) in itr
        tag = convert(AuxTag, k)
        setindex!(y, v, tag)
    end
    y
end

function Base.isvalid(aux::AbstractAuxiliary)
    all(i -> !isa(i, Error), iter_encodings(aux))
end

function Base.copy(aux::AbstractAuxiliary)
    x = aux.x
    v = if x isa Vector{UInt8}
        x[aux.start:end]
    else
        copy(MemView(aux))
    end
    typeof(aux)(v, 1)
end

Base.empty(T::Type{<:AbstractAuxiliary{V}}) where V = T(empty(V), 1)

function Base.length(aux::AbstractAuxiliary)::Int
    n = 0
    for val in iter_encodings(aux)
        val isa Error && return n
        n += 1
    end
    n
end

function Base.haskey(aux::AbstractAuxiliary, k)::Bool
    key = convert(AuxTag, k)
    any(iter_encodings(aux)) do i
        !(i isa Error) && first(i) == key
    end
end

function Base.keys(aux::AbstractAuxiliary)
    Iterators.map(iter_encodings(aux)) do val
        if val isa Error
            error("Bad AuxTag")
        else
            first(val)
        end
    end
end

function setindex_nonexisting! end

const ELTYPE_DICT = Dict(
    UInt8('C') => UInt8,
    UInt8('c') => Int8,
    UInt8('S') => UInt16,
    UInt8('s') => Int16,
    UInt8('I') => UInt32,
    UInt8('i') => Int32,
    UInt8('f') => Float32,
)

is_printable_char(x::UInt8) = in(x, UInt8('!'):UInt8('~'))

module Errors

@enum Error::Int32 begin
    # These four can be emitted from the EncodedIterator
    TooShortMemory # too short memory for any data, or for data type specifically for BAM
    InvalidAuxTag # Tags must conform to r\"[A-Za-z][A-Za-z0-9]
    NoColons # not AB:C:x, SAM only
    NoNullByte # BAM only

    # These are emitted from EncodedIterator in BAM's case,
    # but when loading the aux value in SAM's case.
    InvalidTypeTag
    InvalidArrayEltype

    # Only emitted when loading the value
    InvalidChar # not printable
    InvalidInt # SAM only
    InvalidFloat # SAM only
    InvalidString
    InvalidArray
    InvalidHex
end

end # module errors

using .Errors: Errors, Error

abstract type AbstractEncodedIterator end
function iter_encodings end

Base.IteratorSize(::Type{<:AbstractEncodedIterator}) = Base.SizeUnknown()
Base.eltype(::Type{AbstractEncodedIterator}) = Union{Error, Tuple{AuxTag, UInt8, UnitRange{Int}}}

function load_hex(mem::ImmutableMemView)::Union{Memory{UInt8}, Error}
    len = length(mem)
    # Note: According to specs, Hex can't be empty, but we load it anyway
    # because we should be generous in what we accept
    isodd(len) && return Errors.InvalidHex
    reslen = div(len % UInt, 2) % Int
    hex = Memory{UInt8}(undef, reslen)
    @inbounds for i in 1:reslen
        a = pack_hex(mem[2i - 1])
        b = pack_hex(mem[2i])
        (a | b) == 0xff && return Errors.InvalidHex
        hex[i] = (a << 4) | b
    end
    hex
end

function pack_hex(a::UInt8)::UInt8 # 0xff for bad hex
    if a ∈ 0x30:0x39
        a - 0x30
    elseif a ∈ UInt8('A'):UInt8('F')
        a - UInt8('A') + UInt8(10)
    elseif a ∈ UInt8('a'):UInt8('f')
        a - UInt8('a') + UInt8(10)
    else
        0xff
    end
end

# For this constructor, we can rely on the keys being unique, so we can
# use the more efficient setindex_nonexisting!
function (T::Type{<:AbstractAuxiliary})(d::AbstractDict)
    y = empty(T)
    for (k, v) in pairs(d)
        tag = convert(AuxTag, k)
        setindex_nonexisting!(y, v, tag)
    end
    y
end

function Base.show(io::IO, ::MIME"text/plain", x::AbstractAuxiliary)
    if isvalid(x)
        buf = IOBuffer()
        println(buf, length(x), "-element ", typeof(x), ':')
        content = IOContext(buf, :limit => true, :compact => true)
        for (n, (key, value)) in enumerate(x)
            if n > 20
                println(content, "⋮")
                break
            end
            # TODO: Content length limitation doesn't actually work
            println(content, "  \"", string(key), "\" => ", repr(value))
        end
        write(io, take!(buf)[1:end-1]) # remove trailing newline
    else
        print(io, "Invalid ", typeof(x))
    end
end

# Because checking the size is O(N).
Base.IteratorSize(::Type{<:AbstractAuxiliary}) = Base.SizeUnknown()
Base.eltype(::Type{<:AbstractAuxiliary}) = Pair{AuxTag, Any}

include("delimited.jl")
include("bam.jl")
include("sam.jl")

end # module XAMAuxData
