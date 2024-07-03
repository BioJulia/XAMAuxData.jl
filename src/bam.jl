module BAM

import ..AuxTag, ..AbstractAuxiliary, ..Hex
import ..AUX_NUMBER_TYPES, ..try_auxtag, ..Error, ..Errors
import ..is_printable, ..ELTYPE_DICT, ..load_hex, ..iter_encodings, ..AbstractEncodedIterator
import ..is_printable_char, ..as_bam_aux_value, ..get_type_tag, ..hexencode!, ..AuxException
import ..striptype

public Auxiliary, AuxTag, Error, Errors

using MemoryViews: ImmutableMemoryView, MutableMemoryView, MemoryView
using StringViews: StringView

"""
    BAM.Auxiliary{T <: AbstractVector{UInt8}} <: AbstractDict{AuxTag, Any}

Lazily loaded `AbstractDict` representing the auxiliary data fields of a BAM
record. Immutable aux's can be constructed with `Auxiliary(x)` for any `x`
with `MemoryView(x)` defined.
Mutable aux data is constructed with `Auxiliary(x::Vector{UInt8}, start::Int)`,
where `start` gives the first index of the used data in `x` - all data before
`start` will be ignored and never modified.

# Examples
```jldoctest
julia> immut = BAM.Auxiliary("KJS\7\\0ABZabc\\0");

julia> immut["KJ"]
0x0007

julia> haskey(immut, "AB")
true
```

# Extended help
Since fields of `Auxiliary` are lazily loaded, it may contain invalid data.
This package distinguishes two ways of being invalid: If the data does not conform
to tab-separated fields of `AB:X:Z`, iterating the `Auxiliary` or accessing fields
will error, and the whole `Auxiliary` is considered invalid.
Use isvalid to check if this is the case.
Alternatively, if the _value_ `Z` is invalid, the `Auxiliary` can be iterated and
indexed, but that value will be returned as a value of [`Error`](@ref).

See also: [`Error`](@ref)

# Examples
```jldoctest
julia> invalid = BAM.Auxiliary("KVXA");

julia> isvalid(invalid)
false

julia> invalid["KV"]
ERROR: Unknown type tag in aux value.
[...]

julia> valid_badvalue = BAM.Auxiliary("KVA\f");

julia> isvalid(valid_badvalue)
true

julia> valid_badvalue["KV"] isa Error
true
```
"""
struct Auxiliary{T} <: AbstractAuxiliary{T}
    x::T
    start::Int

    function Auxiliary(v::Vector{UInt8}, i::Integer)
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
    value = load_auxvalue(typetag, @inbounds it.mem[span])
    value isa Error && return (value, new_state)
    (key => value, new_state)
end

struct EncodedIterator <: AbstractEncodedIterator
    mem::ImmutableMemoryView{UInt8}
end

iter_encodings(aux::Auxiliary) = EncodedIterator(ImmutableMemoryView(aux))

function Base.iterate(it::EncodedIterator, state::Int=1)
    mem = it.mem
    state > length(mem) && return nothing
    state > length(mem) - 3 && return (Errors.TooShortMemory, length(mem) + 1)
    t1 = @inbounds mem[state]
    t2 = @inbounds mem[state + 1]
    tag = @something try_auxtag(t1, t2) return (Errors.InvalidAuxTag, length(mem) + 1)
    type_tag = @inbounds mem[state + 2]
    start = state + 3
    data_length = length(mem) - start + 1
    # One byte values
    stop = if type_tag in (UInt8('C'), UInt8('c'), UInt8('A'))
        start
    # Two byte values
    elseif type_tag in (UInt8('S'), UInt8('s'))
        data_length < 2 && return (Errors.TooShortMemory, length(mem)+1)
        start + 1
    # Four byte values
    elseif type_tag in (UInt8('I'), UInt8('i'), UInt8('f'))
        data_length < 4 && return (Errors.TooShortMemory, length(mem)+1)
        start + 3
    # Null-terminated values
    elseif type_tag in (UInt8('Z'), UInt8('H'))
        zeropos = findnext(iszero, mem, start)
        isnothing(zeropos) && return (Errors.NoNullByte, length(mem) + 1)
        zeropos
    # Arrays
    elseif type_tag == UInt8('B')
        # Minimum data length for empty array:
        # Array element type byte plus 4 for array length
        data_length < 5 && return (Errors.TooShortMemory, length(mem) + 1)
        eltype_tag = @inbounds mem[start]
        # Don't use dispatch here for efficiency to avoid type instability
        eltype_size = if eltype_tag in (UInt8('C'), UInt8('c'))
            1
        elseif eltype_tag in (UInt8('S'), UInt8('s'))
            2
        elseif eltype_tag in (UInt8('I'), UInt8('i'), UInt8('f'))
            4
        else
            return (Errors.InvalidArrayEltype, length(mem) + 1)
        end
        # Note: All BAM integers are little endian so we can do this
        n_elements = @inbounds begin
            mem[start + 1] % UInt32 |
            (mem[start + 2] % UInt32) << 8 |
            (mem[start + 3] % UInt32) << 16 |
            (mem[start + 4] % UInt32) << 24
        end
        len = n_elements * eltype_size
        data_length < len + 5 && return (Errors.TooShortMemory, length(mem) + 1)
        start + 4 + len
    else
        return (Errors.InvalidTypeTag, length(mem) + 1)
    end
    ((tag, type_tag, start:stop), stop + 1)
end

function load_array(mem::ImmutableMemoryView{UInt8})
    # This might not be possible to hit in practise.
    length(mem) < 5 && return Errors.InvalidArray
    @inbounds begin
        # The correctness of this byte has already been validated in the EncodedIterator
        eltype = ELTYPE_DICT[mem[1]]
            n_elements = mem[2] % UInt32 |
            (mem[3] % UInt32) << 8 |
            (mem[4] % UInt32) << 16 |
            (mem[5] % UInt32) << 24
    end
    load_array(eltype, n_elements, @inbounds mem[6:end])
end

function load_array(T::Type, n_elements::UInt32, mem::ImmutableMemoryView{UInt8})
    res = reinterpret(T, mem)
    # Should not be possible, since the number of elements is used to determine
    # the memory size
    length(res) == n_elements || return Errors.InvalidArray
    res
end

function load_auxvalue(type_tag::UInt8, mem::ImmutableMemoryView{UInt8})
    if type_tag == UInt8('C')
        @inbounds mem[1]
    elseif type_tag == UInt8('c')
        @inbounds mem[1] % Int8
    elseif type_tag == UInt8('A')
        c = @inbounds mem[1]
        is_printable_char(c) || return Errors.InvalidChar
        Char(c)
    else
        GC.@preserve mem begin
            ptr = pointer(mem)
            if type_tag == UInt8('s')
                ltoh(unsafe_load(Ptr{Int16}(ptr)))
            elseif type_tag == UInt8('S')
                ltoh(unsafe_load(Ptr{UInt16}(ptr)))
            elseif type_tag == UInt8('i')
                ltoh(unsafe_load(Ptr{Int32}(ptr)))
            elseif type_tag == UInt8('I')
                ltoh(unsafe_load(Ptr{UInt32}(ptr)))
            elseif type_tag == UInt8('f')
                ltoh(unsafe_load(Ptr{Float32}(ptr)))
            elseif type_tag == UInt8('Z')
                # Compensate for null terminator byte
                sm = @inbounds mem[1:end-1]
                is_printable(sm) ? StringView(sm) : Errors.InvalidString
            elseif type_tag == UInt8('H')
                mem = @inbounds mem[1:end-1]
                load_hex(mem)
            elseif type_tag == UInt8('B')
                load_array(mem)
            else
                # should be unreachable, has been validated in EncodedIterator
                Errors.InvalidTypeTag
            end
        end
    end
end

function Base.delete!(aux::MutableAuxiliary, k)
    key = convert(AuxTag, k)
    for v in iter_encodings(aux)
        v isa Error && throw(AuxException(v))
        (tag, _, span) = v
        if tag == key
            offset = aux.start - 1
            deleteat!(aux.x, first(span)-3+offset:last(span)+offset)
            break
        end
    end
    aux
end

bytes_needed(x::Union{Int8, UInt8, Char}) = 1
bytes_needed(x::Union{Int16, UInt16}) = 2
bytes_needed(x::Union{Int32, UInt32, Float32}) = 4

bytes_needed(x::AbstractString) = length(MemoryView(codeunits(x))) + 1 # null byte
bytes_needed(x::Hex) = 2 * length(x.x) + 1 # null byte

function bytes_needed(x::AbstractVector{<:AUX_NUMBER_TYPES})
    base = 1 + 4
    isempty(x) ? base : base + bytes_needed(@inbounds x[1]) * length(x) 
end

as_bam_aux_value(x::AUX_NUMBER_TYPES) = x
as_bam_aux_value(x::Signed) = Int32(x)
as_bam_aux_value(x::Unsigned) = UInt32(x)

function Base.setindex!(aux::MutableAuxiliary, val, k)
    key = convert(AuxTag, k)
    for v in iter_encodings(aux)
        v isa Error && error("Cannot set value into invalid Auxiliary")
        (tag, type_tag, span) = v
        if tag == key
            bam_val = as_bam_aux_value(val)
            n_bytes_needed = bytes_needed(bam_val)
            if length(span) > n_bytes_needed
                leftshift = length(span) - n_bytes_needed
                deleteat!(aux.x, first(span)-1:first(span)-2+leftshift)
            # Shift right
            elseif length(span) < n_bytes_needed
                rightshift = n_bytes_needed - length(span)
                oldsize = length(aux.x)
                resize!(aux.x, oldsize + rightshift)
                mem = MemoryView(aux)
                @inbounds for i in length(mem):-1:last(span)+rightshift+1
                    mem[i] = mem[i - rightshift]
                end
            end
            mem = @inbounds MemoryView(aux)[first(span)-1:first(span) + n_bytes_needed - 1]
            write_auxvalue_typetag!(mem, bam_val)
            return aux
        end
    end
    setindex_nonexisting!(aux, val, key)
    aux
end

function setindex_nonexisting!(aux::MutableAuxiliary, val, k)
    key = convert(AuxTag, k)
    bam_val = as_bam_aux_value(val)
    n_bytes_needed = bytes_needed(bam_val)
    data = aux.x
    old_len = length(data)
    resize!(data, old_len + 3 + n_bytes_needed)
    @inbounds data[old_len + 1] = key.x[1]
    @inbounds data[old_len + 2] = key.x[2]
    mem = @inbounds MemoryView(aux)[end - n_bytes_needed: end]
    write_auxvalue_typetag!(mem, bam_val)
    aux
end

function write_auxvalue_typetag!(mem::MutableMemoryView{UInt8}, bam_val)
    type_tag = get_type_tag(typeof(bam_val))
    @inbounds mem[1] = type_tag
    if type_tag in (UInt8('C'), UInt8('c'))
        @inbounds mem[2] = reinterpret(UInt8, bam_val)
    elseif type_tag == UInt8('A')
        @inbounds mem[2] = (reinterpret(UInt32, bam_val) >> 24) % UInt8
    else
        GC.@preserve mem begin
            ptr = pointer(mem) + 1
            if type_tag in (UInt8('S'), UInt8('s'))
                unsafe_store!(Ptr{UInt16}(ptr), htol(reinterpret(UInt16, bam_val)))
            elseif type_tag in (UInt8('i'), UInt8('I'), UInt8('f'))
                unsafe_store!(Ptr{UInt32}(ptr), htol(reinterpret(UInt32, bam_val)))
            elseif type_tag == UInt8('Z')
                unsafe_copyto!(mem[2:end], ImmutableMemoryView(codeunits(bam_val)))
                @inbounds mem[end] = 0x00
            elseif type_tag == UInt8('H')
                m = @inbounds mem[2:end-1]
                hexencode!(m, bam_val)
                @inbounds mem[end] = 0x00
            elseif type_tag == UInt8('B')
                eltype_tag = get_type_tag(eltype(bam_val))
                @inbounds mem[2] = eltype_tag
                unsafe_store!(Ptr{UInt32}(ptr) + 1, htol(length(bam_val) % UInt32))
                ptr = Ptr{eltype(bam_val)}(ptr + 5)
                for num in bam_val
                    unsafe_store!(ptr, htol(num))
                    ptr += sizeof(num)
                end
            else
                error("Unreachable")
            end
        end
    end
end

function Base.get(aux::Auxiliary, k, default)
    key = AuxTag(k)
    it = iter_encodings(aux)
    for i in it
        i isa Error && throw(AuxException(i))
        (auxtag, typetag, span) = i
        if auxtag == key
            return load_auxvalue(typetag, @inbounds it.mem[span])
        end
    end
    default
end

end # module
