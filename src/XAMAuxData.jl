module XAMAuxData

using MemoryViews: MemoryViews, MemoryView, MutableMemoryView, ImmutableMemoryView, RefVector
using StringViews: StringView

struct Unsafe end
const unsafe = Unsafe()

export Hex, AuxTag, try_auxtag, try_get, SAM, BAM, is_well_formed, setindex_nonexisting!
public Errors, Error

# These are the numerical types supported by the BAM format, and by SAM arrays.
const AUX_INTEGER_TYPES = Union{Int8, UInt8, Int16, UInt16, Int32, UInt32}
const AUX_NUMBER_TYPES = Union{AUX_INTEGER_TYPES, Float32}
const BIT_INTEGERS = Union{UInt8, Int8, UInt16, Int16, UInt32, Int32, UInt64, Int64, UInt128, Int128}

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
        return new{typeof(v)}(v)
    end
end

# Must encode to uppercase A-F
hexencode_nibble(u::UInt8)::UInt8 = u < 0x0a ? UInt8('0') + u : UInt8('A') - 0x0a + u

function hexencode!(mem::MutableMemoryView, hex::Hex)
    @inbounds for (byte_no, byte) in enumerate(hex.x)
        mem[2 * byte_no - 1] = hexencode_nibble(byte >> 4)
        mem[2 * byte_no] = hexencode_nibble(byte & 0x0f)
    end
    return mem
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
get_type_tag(::Type{Int64}) = UInt8('i') # SAM only
get_type_tag(::Type{Float32}) = UInt8('f')
get_type_tag(::Type{Char}) = UInt8('A')
get_type_tag(::Type{<:AbstractString}) = UInt8('Z')
get_type_tag(::Type{<:AbstractVector}) = UInt8('B')
get_type_tag(::Type{<:Hex}) = UInt8('H')

function as_aux_value end
as_sam_aux_value(x) = as_aux_value(x)
as_bam_aux_value(x) = as_aux_value(x)

# NaN and Inf are not permitted by the SAM specs, and floats too large to be
# represented as a Float32 would be converted to Inf.
function as_aux_value(x::Real)::Float32
    f = Float32(x)::Float32
    isfinite(f) || throw(AuxException(Errors.InvalidFloat))
    return f
end

function as_aux_value(x::AbstractChar)::Char
    c = Char(x)::Char
    return is_printable_char(c) ? c : throw(AuxException(Errors.InvalidChar))
end

as_aux_value(x::Hex) = x

as_aux_value(s::AbstractString) = as_aux_value(String(s))

function as_aux_value(
        s::Union{String, SubString{String}, StringView, SubString{<:StringView}},
    )
    mem = ImmutableMemoryView(s)
    is_printable(mem) || throw(AuxException(Errors.InvalidString))
    return s
end

# Returns an AbstractVector{<:AUX_NUMBER_TYPES}
function as_aux_value(v::AbstractVector{<:Real})
    E = eltype(v)
    if E === Float32
        all(isfinite, v) || throw(AuxException(Errors.InvalidFloat))
        return v
    end
    E <: AUX_NUMBER_TYPES && isconcretetype(E) && return v
    # We do not choose an integer type on behalf of the user, since the
    # element type is part of the array's encoding.
    if E <: Integer
        throw(
            ArgumentError(
                "Integer arrays must have element type Int8, UInt8, Int16, UInt16, Int32 or UInt32",
            ),
        )
    end
    mem = Memory{Float32}(undef, length(v))
    for i in eachindex(v, mem)
        mem[i] = as_aux_value(v[i]::Real)::Float32
    end
    return mem
end

# Values loaded from an Auxiliary may be views into its data. If such a value is
# written into a mutable Auxiliary, mutating its data before the value is written
# would corrupt the value. To prevent this, we copy values that may share memory
# with the Auxiliary. This is done by finding a memory view which contains all the
# data of the value, and checking whether the memory regions overlap.
# Values for which no such memory view can be found are assumed not to alias.

# Return a `MemoryView` containing all memory of `x`, or `nothing`
memory_of(x::Union{StringView, SubString{<:StringView}}) = ImmutableMemoryView(x)
memory_of(x::Union{Vector, Memory, MemoryView}) = ImmutableMemoryView(x)
memory_of(x::Hex) = memory_of(x.x)

# Note: The parent may contain more memory than `x` itself, so this may produce
# false positives, which only cause a needless copy.
function memory_of(x::AbstractArray)
    p = parent(x)
    return p === x ? nothing : memory_of(p)
end

memory_of(::Any) = nothing

function overlaps(a::MemoryView, b::MemoryView)::Bool
    (isempty(a) | isempty(b)) && return false
    return GC.@preserve a b begin
        pa = UInt(pointer(a))
        pb = UInt(pointer(b))
        (pa < pb + sizeof(b)) & (pb < pa + sizeof(a))
    end
end

overlaps(::MemoryView, ::Nothing) = false

copy_value(x::AbstractString) = String(x)
copy_value(x::AbstractVector) = collect(x)
copy_value(x::Hex) = Hex(collect(x.x))

# Return `v`, or a copy of `v` if it may share memory with `data`.
# We check against all of `data`, not just the part used by the Auxiliary,
# because mutating the Auxiliary may move data before its start.
function unalias_value(data::Vector{UInt8}, v)
    return overlaps(ImmutableMemoryView(data), memory_of(v)) ? copy_value(v) : v
end

function is_printable(v::AbstractVector{UInt8})
    res = true
    for b in v
        res &= (b == UInt8(' ')) | is_printable_char(b)
    end
    return res
end

abstract type AbstractAuxiliary{T} <: AbstractDict{AuxTag, Any} end

# Shrinking the vector of a mutable Auxiliary to below `start - 1` elements is
# not supported, but must cause an error, not out-of-bounds access, when mutating.
function check_not_shrunk(aux::AbstractAuxiliary{Vector{UInt8}})
    aux.start > length(aux.x) + 1 && throw_shrunk_error()
    return nothing
end

@noinline function throw_shrunk_error()
    error("Vector of mutable Auxiliary has been shrunk to below its start index")
end

function striptype end

function Base.copy(aux::AbstractAuxiliary)
    x = aux.x
    T = striptype(typeof(aux))
    return if x isa Vector{UInt8}
        T(x[aux.start:end], 1)
    else
        T(copy(MemoryView(aux)))
    end
end

function Base.length(aux::AbstractAuxiliary)::Int
    n = 0
    for val in iter_encodings(aux)
        val isa Error && throw(AuxException(val))
        n += 1
    end
    return n
end

# Like `get`, this throws if `aux` is malformed before `k` is found.
function Base.haskey(aux::AbstractAuxiliary, k)::Bool
    key = convert(AuxTag, k)
    for i in iter_encodings(aux)
        i isa Error && throw(AuxException(i))
        first(i) == key && return true
    end
    return false
end

function Base.keys(aux::AbstractAuxiliary)
    return Iterators.map(iter_encodings(aux)) do val
        if val isa Error
            throw(AuxException(val))
        else
            first(val)
        end
    end
end

"""
    setindex_nonexisting!(dst::MutableAuxiliary, val, key) -> dst

Same as `dst[key] = val`, but assumes that `key` is not present
in `dst`.

!!! warning
    If `key` is in `dst`, this will result in a malformed `dst` with a duplicate key,
    which may cause invalid behaviour.

# Examples
```jldoctest
julia> v = SAM.Auxiliary(UInt8[], 1);

julia> setindex_nonexisting!(v, 'k', "BA")
1-element XAMAuxData.SAM.Auxiliary{Vector{UInt8}}:
  "BA" => 'k'
```
"""
function setindex_nonexisting! end

is_printable_char(x::UInt8) = in(x, UInt8('!'):UInt8('~'))
is_printable_char(c::Char) = '!' ≤ c ≤ '~'

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
        InvalidFloat
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
The errors are _nonexhaustive_ - more might be added in a non-breaking release.

The following errors may contained in `AbstractAuxiliary` instead of the real
value:
* `InvalidTypeTag` (SAM only): An aux value with an unknown type
* `InvalidArrayEltype` (SAM only): A `B` value with an unknown element type
* `InvalidInt` (SAM only): An integer that is not a decimal number in the range of `Int64`
* `InvalidFloat`: A float that is NaN or infinite, or, in SAM, which is not a decimal
  number representable by a `Float32`
* `InvalidChar`: Loading a `Char` not in `'!':'~'`
* `InvalidString`: A string that contains a character not in `re"[ !-~]"`
* `InvalidHex`: a `H` value with an odd number of symbols, or symbol not in `re"[0-9A-F]"`
* `InvalidArray`: A malformed `B` value, or one containing an invalid element

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
        "SAM integer cannot be parsed as base-ten Int64, or is not in [-2^31, 2^32) when writing or converting to BAM."
    elseif inner == Errors.InvalidFloat
        "Float is NaN, infinite or not representable as a Float32, or is not a decimal number in SAM."
    elseif inner == Errors.InvalidString
        "Auxiliary String (type 'Z') can only contain bytes in re\"[ !-~]\"."
    elseif inner == Errors.InvalidArray
        "Invalid array value. Not enough bytes available, or could not parse data to array element type."
    elseif inner == Errors.InvalidHex
        "Auxiliary Hex (type 'H') contains data not parsable as hexadecimal bytes."
    else
        @assert false # unreachable
    end
    return print(io, s)
end

abstract type AbstractEncodedIterator end
function iter_encodings end

Base.IteratorSize(::Type{<:AbstractEncodedIterator}) = Base.SizeUnknown()
Base.eltype(::Type{<:AbstractEncodedIterator}) = Union{Error, Tuple{AuxTag, UInt8, UnitRange{Int}}}

function load_hex(mem::ImmutableMemoryView)::Union{Memory{UInt8}, Error}
    len = length(mem)
    # Note: According to specs, Hex can't be empty, but we load it anyway
    # because we should be generous in what we accept
    isodd(len) && return Errors.InvalidHex
    reslen = div(len % UInt, 2) % Int
    hex = Memory{UInt8}(undef, reslen)
    # Do not return early on bad bytes, such that this loop is branchless and SIMDs
    bad = false
    @inbounds for i in 1:reslen
        (a, bad_a) = unhex_nibble(mem[2i - 1])
        (b, bad_b) = unhex_nibble(mem[2i])
        bad |= bad_a | bad_b
        hex[i] = (a << 4) | b
    end
    return bad ? Errors.InvalidHex : hex
end

"""
    is_well_formed(aux::AbstractAuxiliary; check_duplicates::Bool = true) -> Bool

Check if the auxiliary is well formed.
`AbstractAuxiliary` distinguishes two ways their data can be incorrect:
If the data is corrupted such that no loaded value can be trusted without checking
the whole auxiliary, we say the auxiliary is malformed.
This is the case if it is not possible to identify the key, or start/end of the
encoded data of the value for a key/value pair, or if a key occurs more than once.
Accessing values, or iterating such an auxiliary _may_ throw an exception, but is not
guaranteed to. For example, loading a duplicated key may return any of its values.
This notion of malformedness is what this function checks for.

Checking for duplicate keys has a significant cost relative to the rest of this
function. If `check_duplicates` is `false`, duplicate keys are not checked for.
Only do this if you know `aux` has no duplicate keys, e.g. because its data was
produced by a trusted source.

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

julia> aux = SAM.Auxiliary("AB:i:1\\tAB:i:2"); # duplicate key

julia> is_well_formed(aux)
false

julia> is_well_formed(aux; check_duplicates = false)
true
```
"""
function is_well_formed(aux::AbstractAuxiliary; check_duplicates::Bool = true)
    # Note: This is faster than calling `all_encodings` with `check_duplicates = false`
    check_duplicates || return all(i -> !isa(i, Error), iter_encodings(aux))
    return all_encodings((_, _) -> true, aux, true)
end

# Check that `aux` is well-formed, i.e. that iterating its encodings yields no
# `Error`, and, if `check_duplicates`, that it contains no duplicate keys,
# and that `f(it, encoding)` returns `true` for all encodings, where `it` is the
# `iter_encodings` iterator.
# Duplicate keys are detected using a bitset with one bit per possible key.
# The bitset is small enough to be stack allocated.
function all_encodings(f, aux::AbstractAuxiliary, check_duplicates::Bool)::Bool
    it = iter_encodings(aux)
    if check_duplicates
        bitset = RefVector{UInt64}(undef, N_TAG_BITSET_WORDS)
        fill!(bitset, zero(UInt64))
        for encoding in it
            encoding isa Error && return false
            add_to_tag_bitset!(bitset, first(encoding)) || return false
            f(it, encoding) || return false
        end
    else
        for encoding in it
            encoding isa Error && return false
            f(it, encoding) || return false
        end
    end
    return true
end

# Map each valid tag to a distinct integer in 0:(59 * 75 - 1).
# The first byte is a letter, and `& 0x3f` maps letters to distinct values in 1:58.
# The second byte is in '0':'z', which subtracting '0' maps to 0:74.
tag_bitset_index(tag::AuxTag) = Int(tag.x[1] & 0x3f) * 75 + Int(tag.x[2] - UInt8('0'))
const N_TAG_BITSET_WORDS = cld(59 * 75, 64)

# Add `tag` to the bitset, returning `false` if it was already present.
# The bitset must not escape, since that prevents stack allocation. Therefore,
# this is inlined, and uses pointers, since a bounds check (which may be forced
# by `--check-bounds=yes`) would make the bitset escape in its error path.
# This is safe, since all valid tags map to an index in the bitset.
@inline function add_to_tag_bitset!(bitset::RefVector{UInt64}, tag::AuxTag)::Bool
    index = tag_bitset_index(tag)
    word_index = (index >>> 6) + 1
    mask = UInt64(1) << (index & 63)
    return GC.@preserve bitset begin
        ptr = pointer(bitset)
        word = unsafe_load(ptr, word_index)
        if iszero(word & mask)
            unsafe_store!(ptr, word | mask, word_index)
            true
        else
            false
        end
    end
end

"""
    try_get(
        aux::AbstractAuxiliary,
        key
    )::Union{Nothing, XAMAuxData.Error, Any}

Convert `key` to [`AuxTag`](@ref), then get the value associated with
`key` in `aux`, without throwing on corrupt data.
Return `nothing` if `key` is not present, and an [`XAMAuxData.Error`](@ref) if
the value is invalid, or if `aux` is malformed before `key` is found.

Unlike checking [`is_well_formed`](@ref) before indexing, this only scans the
auxiliary data until `key` is found.
Hence, if `key` is found, it is not checked whether data after `key` is malformed,
for example, whether `key` occurs again.

# Examples
```jldoctest
julia> aux = SAM.Auxiliary("KM:i:252\tAB:Z:αβγδ\tAK::C");

julia> try_get(aux, "KM")
252

julia> try_get(aux, "AB") # invalid value
InvalidString::Error = 9

julia> try_get(aux, "XX") # malformed before key is found
NoColons::Error = 2

julia> try_get(SAM.Auxiliary("KM:i:252"), "XX") === nothing # key not present
true
```

See also: [`is_well_formed`](@ref)
"""
function try_get end

"""
    isvalid(aux::AbstractAuxiliary; check_duplicates::Bool = true) -> Bool

Check if the auxiliary is well formed, and that all the values can be loaded
without returning an [`XAMAuxData.Error`](@ref).
If `check_duplicates` is `false`, duplicate keys are not checked for.
See [`is_well_formed`](@ref) for details.

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
Base.isvalid(x::AbstractAuxiliary; check_duplicates::Bool = true) = throw(MethodError(isvalid, (x,)))

function validate_hex(mem::ImmutableMemoryView{UInt8})::Bool
    iseven(length(mem)) || return false
    bad = false
    for byte in mem
        bad |= last(unhex_nibble(byte))
    end
    return !bad
end

# Decode a byte in `[0-9A-F]` to its value in 0:15, and return whether the byte
# was bad. Lowercase is not permitted by the SAM specs.
# This is branchless, such that loops using it can SIMD.
@inline function unhex_nibble(u::UInt8)::Tuple{UInt8, Bool}
    u -= UInt8('0')
    # After subtracting '0', 'A':'F' becomes 0x11:0x16
    letter = u > 0x09
    bad = (u > 0x16) | (letter & (u < 0x11))
    return (u - 0x07 * letter, bad)
end

function Base.show(io::IO, ::MIME"text/plain", x::AbstractAuxiliary)
    return if is_well_formed(x)
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
        write(io, take!(buf)[1:(end - 1)]) # remove trailing newline
    else
        print(io, "Malformed ", typeof(x))
    end
end

# Because checking the size is O(N).
Base.IteratorSize(::Type{<:AbstractAuxiliary}) = Base.SizeUnknown()
Base.eltype(::Type{<:AbstractAuxiliary}) = Pair{AuxTag, Any}

include("bam.jl")
include("sam.jl")

# SAM integers are loaded as `Int64`, which cannot be written to BAM
sam_to_bam_value(x) = x

function sam_to_bam_value(x::Int64)::Union{Int32, UInt32}
    typemin(Int32) ≤ x ≤ typemax(Int32) && return x % Int32
    0 ≤ x ≤ typemax(UInt32) && return x % UInt32
    throw(AuxException(Errors.InvalidInt))
end

function Base.copy!(dst::BAM.MutableAuxiliary, src::SAM.Auxiliary)
    return copy_converted!(sam_to_bam_value, dst, src)
end

function Base.copy!(dst::SAM.MutableAuxiliary, src::BAM.Auxiliary)
    return copy_converted!(identity, dst, src)
end

# Replace the content of `dst` with the content of `src`, converting values with `f`.
# To leave `dst` unchanged if an exception is thrown, we append the new fields
# after the existing ones, and only remove the existing ones when all new fields
# have been written.
function copy_converted!(f, dst::AbstractAuxiliary{Vector{UInt8}}, src::AbstractAuxiliary)
    check_not_shrunk(dst)
    data = dst.x
    oldlen = length(data)
    try
        for (auxtag, v) in src
            v isa Error && throw(AuxException(v))
            setindex_nonexisting!(dst, f(v), auxtag)
        end
    catch
        resize!(data, oldlen)
        rethrow()
    end
    if oldlen >= dst.start
        # In SAM, a tab separates the existing and the new fields
        has_tab = dst isa SAM.MutableAuxiliary && length(data) > oldlen
        deleteat!(data, dst.start:(oldlen + has_tab))
    end
    return dst
end

end # module XAMAuxData
