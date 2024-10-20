module XAMAuxData

using MemoryViews: MemoryViews, MemoryView, MutableMemoryView, ImmutableMemoryView
using StringViews: StringView

struct Unsafe end
const unsafe = Unsafe()

export Hex, AuxTag, SAM, BAM, is_well_formed
public Errors, Error

# These are the numerical types supported by the BAM format.
const AUX_NUMBER_TYPES = Union{Int8, UInt8, Int16, UInt16, Int32, UInt32, Float32}

include("auxtag.jl")

"""
    Hex(v::AbstractVector{UInt8})

Wrapper type around a byte vector. When this type is assigned to an Auxiliary,
it is encoded as a hex value (type tag `H`) as opposed to a raw byte array (
type tag `B`).

```jldoctest
julia> aux = SAM.Auxiliary(UInt8[], 1);

julia> aux["AB"] = Hex([0xae, 0xf8, 0x6c]);

julia> String(MemoryView(aux))
"AB:H:AEF86C"

julia> aux = BAM.Auxiliary(UInt8[], 1);

julia> aux["AB"] = Hex([0xae, 0xf8, 0x6c]);

julia> String(MemoryView(aux))
"ABHAEF86C\\0"
```
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

function hexencode!(mem::MutableMemoryView, hex::Hex)
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

as_aux_type(::Type{<:Real}) = Float32
as_aux_type(::Type{<:Integer}) = Int32

function as_aux_value end
as_sam_aux_value(x) = as_aux_value(x)
as_bam_aux_value(x) = as_aux_value(x)

as_aux_value(x::Real) = Float32(x)::Float32

function as_aux_value(x::AbstractChar)::Union{Char, Error}
    c = Char(x)::Char
    isascii(c) ? c : throw(AuxException(Errors.InvalidChar))
end

as_aux_value(x::Hex) = x

function as_aux_value(s::AbstractString)
    cu = codeunits(s)
    auxs = if MemoryViews.MemoryKind(typeof(cu)) isa MemoryViews.IsMemory
        s
    else
        String(s)
    end
    # Take view of codeunits, because StringViews' codeunits return
    # the underlying array, so this makes it work for StringViews,
    # without having to implement MemoryView(::StringView), which would
    # be piracy in this package.
    mem = ImmutableMemoryView(codeunits(auxs))
    if is_printable(mem)
        auxs
    else
        throw(AuxException(Errors.InvalidString))
    end
    auxs
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

function striptype end
function Base.copy(aux::AbstractAuxiliary)
    x = aux.x
    v = if x isa Vector{UInt8}
        x[aux.start:end]
    else
        copy(MemoryView(aux))
    end
    striptype(typeof(aux))(v, 1)
end

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
            throw(AuxException(val))
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

"""
    Errors

Module containing the error values of XAMAuxData.

See also: [`Error`](@ref)
"""
module Errors

# See the related showerror method for descriptions of these values
@enum Error::Int32 begin
    # These four can be emitted from the EncodedIterator
    TooShortMemory
    InvalidAuxTag
    NoColons # not AB:C:x, SAM only
    NoNullByte # BAM only

    # These are emitted from EncodedIterator in BAM's case,
    # but when loading the aux value in SAM's case.
    InvalidTypeTag
    InvalidArrayEltype

    # Only emitted when loading the value
    InvalidChar
    InvalidInt # SAM only
    InvalidFloat # SAM only
    InvalidString
    InvalidArray
    InvalidHex
end

end # module errors

using .Errors: Errors, Error

"""
    Error

Enum type representing errors returned when loading invalid auxiliary data
values.
The errors are _nonexhausitve_ - more might be added in a non-breaking release.

The following errors may contained in `AbstractAuxiliary` instead of the real
value:
* `InvalidTypeTag` (SAM only): An aux value with an unknown type
* `InvalidArrayEltype` (SAM only): A `B` value with an unknown element type
* `InvalidInt` (SAM only): An integer that can't be parsed to an `Int32`
* `InvalidFloat` (SAM only): A float that can't be parsed to a `Float32`
* `InvalidChar`: Loading a `Char` not in `'!':'~'`
* `InvalidString`: A string that contains a character not in `re"[ !-~]"`
* `InvalidHex`: a `H` value with an odd number of symbols, or symbol not in `re"[0-9A-F]"`
* `InvalidArray`: A malformed `B` value

See also: [`Errors`](@ref)
"""
Error

struct AuxException <: Exception
    error::Error
end

function Base.showerror(io::IO, error::AuxException)
    inner = error.error
    s = if inner == Errors.TooShortMemory
        "Not enough bytes remaining in Aux data to encode a full field"
    elseif inner == Errors.InvalidAuxTag
        "Invalid AuxTag. Tags must conform to r\"^[A-Za-z][A-Za-z0-9]\$\"."
    elseif inner == Errors.NoColons
        "Invalid SAM tag header. Expected <AuxTag>:<type tag>:, but found no colons."
    elseif inner == Errors.NoNullByte
        "BAM string or Hex type not terminated by null byte"
    elseif inner == Errors.InvalidTypeTag
        "Unknown type tag in aux value."
    elseif inner == Errors.InvalidArrayEltype
        "Invalid array element type tag. Valid values are CcSsIif."
    elseif inner == Errors.InvalidChar
        "Auxiliary Char (type 'A') must be in '!':'~'."
    elseif inner == Errors.InvalidInt
        "Data in SAM Auxiliary cannot be parsed as base-ten Int32."
    elseif inner == Errors.InvalidFloat
        "Data in SAM Auxiliary cannot be parsed as Float32."
    elseif inner == Errors.InvalidString
        "Auxiliary String (type 'Z') can only contain bytes in re\"[ !-~]\"."
    elseif inner == Errors.InvalidArray
        "Invalid array value. Not enough bytes available, or could not parse data to array element type."
    elseif inner == Errors.InvalidHex
        "Auxiliary Hex (type 'H') contains data not parsable as hexidecimal bytes."
    else
        @assert false # unreachable
    end
    print(io, s)
end

abstract type AbstractEncodedIterator end
function iter_encodings end

Base.IteratorSize(::Type{<:AbstractEncodedIterator}) = Base.SizeUnknown()
Base.eltype(::Type{AbstractEncodedIterator}) = Union{Error, Tuple{AuxTag, UInt8, UnitRange{Int}}}

function load_hex(mem::ImmutableMemoryView)::Union{Memory{UInt8}, Error}
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

"""
    is_well_formed(aux::AbstractAuxiliary) -> Bool

Check if the auxiliary is well formed.
`AbstractAuxiliary` distinguishes two ways their data can be incorrect:
If the data is corrupted such that it is not possible to identify the key, or
start/end of the encoded data of the value for a key/value pair, we say the
auxiliary is malformed.
Accessing values, or iterating such an auxiliary may or may not throw an exception.
This notion of malformedness is what this function checks for.

Alternatively, if the keys and encoded values _can_ be identified, but the encoded
values are corrupt and the value cannot be loaded, we say the auxiliary is invalid.
Such auxiliaries can be iterated over, and the values can be loaded, although the
loaded value will be of type [`XAMAuxData.Error`](@ref).
This can be checked for with `isvalid`.
Valid auxiliaries are always well-formed.

# Examples
```jldoctest
julia> aux = SAM.Auxiliary("KM:i:252\\tAK::C"); # bad key

julia> is_well_formed(aux)
false

julia> aux["KM"] # loading a value may error
252

julia> aux["AK"] # loading a value may error
ERROR: Invalid SAM tag header. Expected <AuxTag>:<type tag>:, but found no colons.
[...]

julia> aux = SAM.Auxiliary("AB:Z:αβγδ\\tCD:f:-1.2"); # bad value encoding

julia> (is_well_formed(aux), isvalid(aux))
(true, false)

julia> aux["AB"] # note: numerical value of enum is not API
InvalidString::Error = 9
```
"""
function is_well_formed(it::AbstractAuxiliary)
    all(iter_encodings(it)) do i
        !isa(i, Error)
    end
end

"""
    isvalid(aux::AbstractAuxiliary) -> Bool

Check if the auxiliary is well formed, and that all the values can be loaded
without returning an [`XAMAuxData.Error`](@ref).

# Examples
```jldoctest
julia> aux = SAM.Auxiliary("AB:i:not an integer");

julia> is_well_formed(aux)
true

julia> isvalid(aux)
false
```

See also: [`is_well_formed`](@ref)
"""
Base.isvalid(x::AbstractAuxiliary) = throw(MethodError(isvalid, (x,)))

function validate_hex(mem::ImmutableMemoryView{UInt8})::Bool
    len = length(mem)
    iseven(len) || return false
    good = true
    for byte in mem
        good &= byte in UInt8('0'):UInt8('9') ||
        byte in UInt8('a'):UInt8('h') ||
        byte in UInt8('A'):UInt8('H')
    end
    good
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

function Base.show(io::IO, ::MIME"text/plain", x::AbstractAuxiliary)
    if is_well_formed(x)
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
        print(io, "Malformed ", typeof(x))
    end
end

# Because checking the size is O(N).
Base.IteratorSize(::Type{<:AbstractAuxiliary}) = Base.SizeUnknown()
Base.eltype(::Type{<:AbstractAuxiliary}) = Pair{AuxTag, Any}

include("delimited.jl")
include("bam.jl")
include("sam.jl")

end # module XAMAuxData
