"""
    SAM

Submodule of `XAMAuxData` with functionality for reading an writing SAM
auxiliary data.
See the public and exported names in this module for its API.
"""
module SAM

# Default Julia methods throw with bad keys, but not bad values (which are simply an error value)

import ..AuxTag, ..AbstractAuxiliary, ..ELTYPE_DICT, ..is_printable_char, ..is_printable, ..Hex, ..setindex_nonexisting!
import ..DelimitedIterator, ..get_type_tag, ..Error, ..Errors, ..load_hex, ..validate_hex
import ..try_auxtag, ..Unsafe, ..as_sam_aux_value, ..AUX_NUMBER_TYPES, ..hexencode!
import ..iter_encodings, ..AbstractEncodedIterator, ..AuxException, ..striptype
import ..is_well_formed

public Auxiliary, AuxTag, Hex, Errors, Error

using MemoryViews: MemoryView, ImmutableMemoryView
using StringViews: StringView

struct EncodedIterator <: AbstractEncodedIterator
    x::DelimitedIterator{UInt8}
end

function Base.iterate(it::EncodedIterator, state::Int=1)
    itval = iterate(it.x, state)
    isnothing(itval) && return nothing
    (mem, newstate) = itval
    parsed = parse_encoded_aux(mem)
    parsed isa Error && return (parsed, length(it.x.v) + 1)
    (tag, eltype) = parsed
    span = 5 + state:lastindex(mem) + state - 1
    ((tag, eltype, span), newstate)
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
    (tag, eltype)
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
"""
struct Auxiliary{T <: AbstractVector{UInt8}} <: AbstractAuxiliary{T}
    x::T
    start::Int

    function Auxiliary(v::Vector{UInt8}, i::Integer)
        i = Int(i)::Int
        ((i - 1) % UInt) > (length(v) % UInt) + 1 && error("Start index must be in 1:length(vector) + 1")
        new{Vector{UInt8}}(v, Int(i)::Int)
    end

    function Auxiliary(x)
        mem = ImmutableMemoryView(x)
        eltype(mem) == UInt8 || error("Must construct Auxiliary from MemoryView{UInt8}")
        new{ImmutableMemoryView{UInt8}}(mem, 1)
    end
end

striptype(::Type{<:Auxiliary}) = Auxiliary
Base.empty(::Type{Auxiliary}) = Auxiliary(UInt8[], 1)

const MutableAuxiliary = Auxiliary{Vector{UInt8}}

function iter_encodings(aux::Auxiliary)
    EncodedIterator(DelimitedIterator(ImmutableMemoryView(aux), UInt8('\t')))
end

MemoryView(x::Auxiliary) = @inbounds MemoryView(x.x)[x.start:end]

function Base.empty!(x::MutableAuxiliary)
    resize!(x.x, x.start - 1)
    x
end

Base.isempty(x::Auxiliary) = x.start > length(x.x)

function Base.iterate(aux::Auxiliary, state::Int=1)
    it = iter_encodings(aux)
    itval = iterate(it, state)
    itval === nothing && return nothing
    (val, new_state) = itval
    val isa Error && throw(AuxException(val))
    (key, typetag, span) = val
    value = load_auxvalue(typetag, @inbounds it.x.v[span])
    (key => value, new_state)
end

function Base.get(aux::Auxiliary, k, default)
    key = AuxTag(k)
    it = iter_encodings(aux)
    for i in it
        i isa Error && throw(AuxException(i))
        (auxtag, typetag, span) = i
        if auxtag == key
            return load_auxvalue(typetag, @inbounds it.x.v[span])
        end
    end
    default
end

function Base.isvalid(aux::Auxiliary)
    it = iter_encodings(aux)
    for i in it
        i isa Error && return false
        (_, eltype, span) = i
        mem = it.x.v[span]
        validate_encoding(eltype, mem) || return false
    end
    true
end

function load_array(mem::ImmutableMemoryView{UInt8})::Union{Memory, Error}
    isempty(mem) && return Errors.InvalidArrayEltype
    eltype_tag = @inbounds mem[1]
    eltype = get(ELTYPE_DICT, eltype_tag, nothing)
    isnothing(eltype) && return Errors.InvalidArrayEltype
    load_array(eltype, @inbounds mem[2:end])
end

function load_array(T::Type, mem::ImmutableMemoryView{UInt8})::Union{Memory, Error}
    isempty(mem) && return Memory{T}()
    length(mem) == 1 && return Errors.InvalidArray
    len = count(==(UInt8(',')), mem)
    @inbounds(mem[1]) == UInt8(',') || return Errors.InvalidArray
    res = Memory{T}(undef, len)
    n = 0
    for ele_mem in DelimitedIterator(@inbounds(mem[2:end]), UInt8(','))
        n += 1
        val = tryparse(T, StringView(ele_mem))
        val === nothing && return Errors.InvalidArray
        @inbounds res[n] = val

    end
    res
end

function validate_array(mem::ImmutableMemoryView{UInt8})::Bool
    isempty(mem) && return false
    eltype_tag = @inbounds mem[1]
    eltype = get(ELTYPE_DICT, eltype_tag, nothing)
    isnothing(eltype) && return false
    validate_array(eltype, @inbounds mem[2:end])
end

function validate_array(T::Type, mem::ImmutableMemoryView{UInt8})::Bool
    isempty(mem) && return true
    length(mem) == 1 && return false
    @inbounds(mem[1]) == UInt8(',') || return false
    all(DelimitedIterator(@inbounds(mem[2:end]), UInt8(','))) do bytes
        !isnothing(tryparse(T, bytes))
    end
end

function validate_encoding(type_tag::UInt8, mem::ImmutableMemoryView{UInt8})::Bool
    if type_tag == UInt8('A')
        length(mem) == 1 || return false
        b = @inbounds mem[1]
        is_printable_char(b)
    elseif type_tag == UInt8('i')
        tryparse(Int, StringView(mem); base=10) !== nothing
    elseif type_tag == UInt8('f')
        tryparse(Float32, StringView(mem)) !== nothing
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
        isempty(mem) && return Errors.InvalidInt
        n = tryparse(Int, StringView(mem); base=10)
        n === nothing ? Errors.InvalidInt : n
    elseif type_tag == UInt8('f')
        n = tryparse(Float32, StringView(mem))
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

as_sam_aux_value(x::Integer) = Int32(x)::Int32

function setindex_nonexisting!(aux::MutableAuxiliary, val, key::AuxTag)
    v = as_sam_aux_value(val)
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
    unsafe_copyto!(MemoryView(data)[oldlen+6:end], MemoryView(value_bytes)[1:length(value_bytes)])
    aux
end

function as_serialized_bytes(val::Union{Char, Real, AbstractString})
    buf = IOBuffer()
    print(buf, val)
    take!(buf)
end

function as_serialized_bytes(hex::Hex)
    m = Memory{UInt8}(undef, 2 * length(hex.x))
    hexencode!(MemoryView(m), hex)
    m
end 

function as_serialized_bytes(v::AbstractVector)
    buf = IOBuffer()
    write(buf, get_type_tag(eltype(v)))
    for i in v
        print(buf, ',')
        print(buf, i)
    end
    take!(buf)
end

function Base.delete!(aux::MutableAuxiliary, k)
    key = convert(AuxTag, k)
    for v in iter_encodings(aux)
        v isa Error && throw(AuxException(v))
        (tag, _, span) = v
        if tag == key
            # Delete also an adjecent tab if it's not the only element
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
                    start_delete:stop_delete + 1
                end
            else
                # Otherwise, delete leading tab
                start_delete-1:stop_delete
            end
            deleteat!(aux.x, to_delete)
            break
        end
    end
    aux
end

function Base.setindex!(aux::MutableAuxiliary, val, k)
    key = convert(AuxTag, k)
    delete!(aux, key)
    setindex_nonexisting!(aux, val, key)
end

end # module
