module SAM

# Default Julia methods throw with bad keys, but not bad values (which are simply an error value)

import ..AuxTag, ..AbstractAuxiliary, ..is_printable_char, ..is_printable, ..Hex, ..setindex_nonexisting!
import ..get_type_tag, ..Error, ..Errors, ..load_hex, ..validate_hex
import ..try_auxtag, ..as_sam_aux_value, ..AUX_NUMBER_TYPES, ..hexencode!
import ..iter_encodings, ..AbstractEncodedIterator, ..AuxException, ..striptype
import ..is_well_formed, ..BIT_INTEGERS, ..load_encoded, ..unalias_value, ..all_encodings, ..check_not_shrunk

public Auxiliary, AuxTag, Hex, Errors, Error

using MemoryViews: MemoryViews, MemoryView, ImmutableMemoryView, DelimitedIterator, Immutable, split_each
using StringViews: StringView

struct EncodedIterator <: AbstractEncodedIterator
    x::DelimitedIterator{UInt8, Immutable, UInt8}
end

function Base.iterate(it::EncodedIterator, state::Int = 1)
    itval = iterate(it.x, state)
    isnothing(itval) && return nothing
    (mem, newstate) = itval
    parsed = parse_encoded_aux(mem)
    # Note: `it.x` yields a trailing empty element at state `length + 1`,
    # so we must return a state beyond that to stop iteration.
    parsed isa Error && return (parsed, length(it.x.v) + 2)
    (tag, eltype) = parsed
    span = (5 + state):(lastindex(mem) + state - 1)
    return ((tag, eltype, span), newstate)
end

function parse_encoded_aux(
        mem::ImmutableMemoryView{UInt8},
    )::Union{Tuple{AuxTag, UInt8}, Error}
    length(mem) < 5 && return Errors.TooShortMemory
    t1 = @inbounds mem[1]
    t2 = @inbounds mem[2]
    tag = try_auxtag(t1, t2)
    isnothing(tag) && return Errors.InvalidAuxTag
    (@inbounds mem[3] == mem[5] == UInt8(':')) || return Errors.NoColons
    eltype = @inbounds mem[4]
    return (tag, eltype)
end

"""
    SAM.Auxiliary{T <: AbstractVector{UInt8}} <: AbstractDict{AuxTag, Any}

Lazily loaded `AbstractDict` representing the auxiliary data fields of a SAM
record. Immutable aux's can be constructed with `Auxiliary(x)` for any `x`
with `MemoryView(x)` defined.
Mutable aux data is constructed with `Auxiliary(x::Vector{UInt8}, start::Int)`,
where `start` gives the first index of the used data in `x` - all data before
`start` will be ignored and never modified.

# Examples
```jldoctest
julia> immut = SAM.Auxiliary("KJ:i:-1\tAB:Z:abc");

julia> immut["KJ"]
-1

julia> haskey(immut, "AB")
true
```

# Extended help
Since fields of `Auxiliary` are lazily loaded, auxiliaries may contain invalid
data even after successful construction.
`Auxiliary` operates with two distinct notions of invalid data - malformedness
and invalidness. See [`is_well_formed`](@ref) for their definitions.

Values loaded from a mutable `Auxiliary` may be views into its data, and are
invalidated by mutating the `Auxiliary`.
See the section "Values may be views into the auxiliary" in the documentation.
"""
struct Auxiliary{T <: AbstractVector{UInt8}} <: AbstractAuxiliary{T}
    x::T
    start::Int

    function Auxiliary(v::Vector{UInt8}, i::Integer)
        i = Int(i)::Int
        ((i - 1) % UInt) > (length(v) % UInt) && error("Start index must be in 1:length(vector) + 1")
        return new{Vector{UInt8}}(v, Int(i)::Int)
    end

    function Auxiliary(x)
        mem = ImmutableMemoryView(x)
        eltype(mem) == UInt8 || error("Must construct Auxiliary from MemoryView{UInt8}")
        return new{ImmutableMemoryView{UInt8}}(mem, 1)
    end
end

striptype(::Type{<:Auxiliary}) = Auxiliary
Base.empty(::Type{Auxiliary}) = Auxiliary(UInt8[], 1)

const MutableAuxiliary = Auxiliary{Vector{UInt8}}

function iter_encodings(aux::Auxiliary)
    return EncodedIterator(split_each(ImmutableMemoryView(aux), UInt8('\t')))
end

MemoryViews.MemoryView(x::Auxiliary) = @inbounds MemoryView(x.x)[x.start:end]

function Base.empty!(x::MutableAuxiliary)
    check_not_shrunk(x)
    resize!(x.x, x.start - 1)
    return x
end

Base.isempty(x::Auxiliary) = x.start > length(x.x)

function Base.iterate(aux::Auxiliary, state::Int = 1)
    it = iter_encodings(aux)
    itval = iterate(it, state)
    itval === nothing && return nothing
    (val, new_state) = itval
    val isa Error && throw(AuxException(val))
    (key, typetag, span) = val
    value = load_auxvalue(typetag, @inbounds it.x.v[span])
    return (key => value, new_state)
end

function load_encoded(aux::Auxiliary, (_, typetag, span)::Tuple{AuxTag, UInt8, UnitRange{Int}})
    return load_auxvalue(typetag, @inbounds ImmutableMemoryView(aux)[span])
end

function Base.isvalid(aux::Auxiliary; check_duplicates::Bool = true)
    return all_encodings(aux, check_duplicates) do it, (_, eltype, span)
        validate_encoding(eltype, @inbounds it.x.v[span])
    end
end

# Call `f(T)`, where `T` is the array element type given by `eltype_tag`,
# or return `default` if the tag is invalid.
# Branch on the tag instead of looking up the element type in a dict,
# such that the element type is known statically, which is required
# for compiling with `juliac --trim`.
@inline function with_array_eltype(f, eltype_tag::UInt8, default)
    return if eltype_tag == UInt8('C')
        f(UInt8)
    elseif eltype_tag == UInt8('c')
        f(Int8)
    elseif eltype_tag == UInt8('S')
        f(UInt16)
    elseif eltype_tag == UInt8('s')
        f(Int16)
    elseif eltype_tag == UInt8('I')
        f(UInt32)
    elseif eltype_tag == UInt8('i')
        f(Int32)
    elseif eltype_tag == UInt8('f')
        f(Float32)
    else
        default
    end
end

# Julia's integer parser accepts more than the SAM specs, e.g. surrounding whitespace,
# and is slower. Hence, we parse integers matching the regex `[-+]?[0-9]+` by hand.
function parse_number(
        ::Type{T}, mem::ImmutableMemoryView{UInt8}
    )::Union{T, Nothing} where {T <: Union{Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64}}
    len = length(mem)
    iszero(len) && return nothing
    b = @inbounds mem[1]
    neg = b == UInt8('-')
    i = 1 + (neg | (b == UInt8('+')))
    i > len && return nothing
    # Accumulate the absolute value, which fits in an UInt64 for all `T`
    n = UInt64(0)
    overflow = false
    @inbounds while i <= len
        d = mem[i] - UInt8('0')
        d > 0x09 && return nothing
        (n, o1) = Base.mul_with_overflow(n, UInt64(10))
        (n, o2) = Base.add_with_overflow(n, d % UInt64)
        overflow |= o1 | o2
        i += 1
    end
    overflow && return nothing
    return if neg
        # The absolute value of typemin(T) is typemax(T) + 1
        n > (T <: Signed ? (typemax(T) % UInt64) + 1 : UInt64(0)) && return nothing
        (UInt64(0) - n) % T
    else
        n > typemax(T) % UInt64 && return nothing
        n % T
    end
end

# Julia's float parser accepts more than the SAM specs, e.g. surrounding whitespace,
# hex floats, NaN and Inf. Hence, check the input is permitted by the SAM specs' regex
# before parsing.
function parse_number(::Type{Float32}, mem::ImmutableMemoryView{UInt8})::Union{Float32, Nothing}
    is_sam_float(mem) || return nothing
    # This may still return nothing, if the value is not representable as a Float32
    return tryparse(Float32, StringView(mem))
end

# Check if `mem` matches the regex `[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?`
function is_sam_float(mem::ImmutableMemoryView{UInt8})::Bool
    i = skip_sign(mem, 1)
    (i, n_int) = skip_digits(mem, i)
    if i <= length(mem) && @inbounds(mem[i]) == UInt8('.')
        (i, n_frac) = skip_digits(mem, i + 1)
        iszero(n_frac) && return false
    else
        iszero(n_int) && return false
    end
    i > length(mem) && return true
    @inbounds(mem[i]) in (UInt8('e'), UInt8('E')) || return false
    (i, n_exp) = skip_digits(mem, skip_sign(mem, i + 1))
    return !iszero(n_exp) && i > length(mem)
end

function skip_sign(mem::ImmutableMemoryView{UInt8}, i::Int)::Int
    i > length(mem) && return i
    b = @inbounds mem[i]
    return i + ((b == UInt8('-')) | (b == UInt8('+')))
end

# Return the index after the digits starting at `i`, and the number of digits
function skip_digits(mem::ImmutableMemoryView{UInt8}, i::Int)::Tuple{Int, Int}
    start = i
    while i <= length(mem) && @inbounds(mem[i]) in UInt8('0'):UInt8('9')
        i += 1
    end
    return (i, i - start)
end

function load_array(mem::ImmutableMemoryView{UInt8})::Union{Memory, Error}
    isempty(mem) && return Errors.InvalidArrayEltype
    eltype_tag = @inbounds mem[1]
    data = @inbounds mem[2:end]
    return with_array_eltype(T -> load_array(T, data), eltype_tag, Errors.InvalidArrayEltype)
end

function load_array(::Type{T}, mem::ImmutableMemoryView{UInt8})::Union{Memory, Error} where {T}
    isempty(mem) && return Memory{T}()
    length(mem) == 1 && return Errors.InvalidArray
    len = count(==(UInt8(',')), mem)
    @inbounds(mem[1]) == UInt8(',') || return Errors.InvalidArray
    res = Memory{T}(undef, len)
    n = 0
    for ele_mem in split_each(@inbounds(mem[2:end]), UInt8(','))
        n += 1
        val = parse_number(T, ele_mem)
        val === nothing && return Errors.InvalidArray
        @inbounds res[n] = val

    end
    return res
end

function validate_array(mem::ImmutableMemoryView{UInt8})::Bool
    isempty(mem) && return false
    eltype_tag = @inbounds mem[1]
    data = @inbounds mem[2:end]
    return with_array_eltype(T -> validate_array(T, data), eltype_tag, false)
end

function validate_array(::Type{T}, mem::ImmutableMemoryView{UInt8})::Bool where {T}
    isempty(mem) && return true
    length(mem) == 1 && return false
    @inbounds(mem[1]) == UInt8(',') || return false
    return all(split_each(@inbounds(mem[2:end]), UInt8(','))) do bytes
        !isnothing(parse_number(T, bytes))
    end
end

function validate_encoding(type_tag::UInt8, mem::ImmutableMemoryView{UInt8})::Bool
    return if type_tag == UInt8('A')
        length(mem) == 1 || return false
        b = @inbounds mem[1]
        is_printable_char(b)
    elseif type_tag == UInt8('i')
        parse_number(Int64, mem) !== nothing
    elseif type_tag == UInt8('f')
        parse_number(Float32, mem) !== nothing
    elseif type_tag == UInt8('Z')
        is_printable(mem)
    elseif type_tag == UInt8('H')
        validate_hex(mem)
    elseif type_tag == UInt8('B')
        validate_array(mem)
    else
        false
    end
end

function load_auxvalue(type_tag::UInt8, mem::ImmutableMemoryView{UInt8})
    return if type_tag == UInt8('A')
        length(mem) == 1 || return Errors.InvalidChar
        b = @inbounds mem[1]
        is_printable_char(b) ? Char(b) : Errors.InvalidChar
    elseif type_tag == UInt8('i')
        # Note: We only allow writing Int32 to align with SAM specs,
        # but we allow reading Int64, in order to be liberal with what
        # we parse
        n = parse_number(Int64, mem)
        n === nothing ? Errors.InvalidInt : n
    elseif type_tag == UInt8('f')
        n = parse_number(Float32, mem)
        n === nothing ? Errors.InvalidFloat : n
    elseif type_tag == UInt8('Z')
        is_printable(mem) ? StringView(mem) : Errors.InvalidString
    elseif type_tag == UInt8('H')
        load_hex(mem)
    elseif type_tag == UInt8('B')
        load_array(mem)
    else
        Errors.InvalidTypeTag
    end
end

function as_sam_aux_value(::Integer)
    throw(
        ArgumentError(
            "SAM integers must be of type Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64, Int128 or UInt128",
        ),
    )
end

# The SAM specs do not limit the size of integers, but recommends the range
# [-2^31, 2^32), since only this range can be represented in BAM.
function as_sam_aux_value(x::BIT_INTEGERS)::Int64
    return if (x < typemin(Int32)) | (x > typemax(UInt32))
        throw(AuxException(Errors.InvalidInt))
    else
        x % Int64
    end
end

function setindex_nonexisting!(aux::MutableAuxiliary, val, key)
    check_not_shrunk(aux)
    key = convert(AuxTag, key)::AuxTag
    return append_converted!(aux, unalias_value(aux.x, as_sam_aux_value(val)), key)
end

# `v` must already have been converted with `as_sam_aux_value`
function append_converted!(aux::MutableAuxiliary, v, key::AuxTag)
    type_tag = get_type_tag(typeof(v))
    value_bytes = as_serialized_bytes(v)
    write_tab = !isempty(aux)
    oldlen = length(aux.x)
    resize!(aux.x, oldlen + 5 + write_tab + length(value_bytes))
    data = aux.x
    @inbounds begin
        if write_tab
            data[oldlen + 1] = UInt8('\t')
            oldlen += 1
        end
        data[oldlen + 1] = key.x[1]
        data[oldlen + 2] = key.x[2]
        data[oldlen + 3] = UInt8(':')
        data[oldlen + 4] = type_tag
        data[oldlen + 5] = UInt8(':')
    end
    unsafe_copyto!(MemoryView(data)[(oldlen + 6):end], MemoryView(value_bytes)[1:length(value_bytes)])
    return aux
end

function as_serialized_bytes(val::Union{Char, Real, AbstractString})
    buf = IOBuffer()
    print(buf, val)
    return take!(buf)
end

function as_serialized_bytes(hex::Hex)
    m = Memory{UInt8}(undef, 2 * length(hex.x))
    hexencode!(MemoryView(m), hex)
    return m
end

function as_serialized_bytes(v::AbstractVector)
    buf = IOBuffer()
    write(buf, get_type_tag(eltype(v)))
    for i in v
        print(buf, ',')
        print(buf, i)
    end
    return take!(buf)
end

function Base.delete!(aux::MutableAuxiliary, k)
    check_not_shrunk(aux)
    key = convert(AuxTag, k)
    for v in iter_encodings(aux)
        v isa Error && throw(AuxException(v))
        (tag, _, span) = v
        if tag == key
            # Delete also an adjacent tab if it's not the only element
            # And the leading 5 bytes e.g. "AB:i:"
            offset = aux.start - 1
            start_delete = first(span) - 5 + offset
            stop_delete = last(span) + offset
            to_delete = if first(span) == 6
                if last(span) + offset == lastindex(aux.x)
                    # If sole element, delete no tabs (because there are none)
                    start_delete:stop_delete
                else
                    # Otherwise, if first element, delete trailing tab
                    start_delete:(stop_delete + 1)
                end
            else
                # Otherwise, delete leading tab
                (start_delete - 1):stop_delete
            end
            deleteat!(aux.x, to_delete)
            break
        end
    end
    return aux
end

function Base.setindex!(aux::MutableAuxiliary, val, k)
    check_not_shrunk(aux)
    key = convert(AuxTag, k)
    # Convert before deleting, such that the existing value is kept if `val` is invalid,
    # and unalias, such that deleting does not corrupt `val`.
    v = unalias_value(aux.x, as_sam_aux_value(val))
    delete!(aux, key)
    return append_converted!(aux, v, key)
end

end # module
