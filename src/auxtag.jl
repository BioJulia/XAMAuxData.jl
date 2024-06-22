"""
    AuxTag(::Union{String, SubString{String}})
    AuxTag(::Char, ::Char)
    AuxTag(::UInt8, ::UInt8)

The type of a key of `Auxiliary`. This type is lightweight, and validates
that the key is of appropriate format.
`Auxiliary` can also be modified using strings as keys, which will be converted
to `AuxTag`, but this may cause a needless allocation of the string.

# Examples
```jldoctest
julia> AuxTag("AC")
AuxTag("AC")

julia> AuxTag('k', 'w')
AuxTag("kw")

julia> AuxTag("1C") # invalid tag
ERROR: Invalid AuxTag. Tags must conform to r"^[A-Za-z][A-Za-z0-9]\$".
[...]
```
"""
struct AuxTag
    x::Tuple{UInt8, UInt8}

    function AuxTag(::Unsafe, x::UInt8, y::UInt8)
        new((x, y))
    end
end

# Makes testing easier
Base.isless(x::AuxTag, y::AuxTag) = isless(x.x, y.x)
Base.isequal(x::AuxTag, y::Union{String, SubString{String}}) = x == try_auxtag(y)

AuxTag(x) = @something try_auxtag(x) throw(AuxException(Errors.InvalidAuxTag))
AuxTag(a, b) = @something try_auxtag(a, b) throw(AuxException(Errors.InvalidAuxTag))

# Per BAM specs, tags must conform to #[A-Za-z][A-Za-z0-9]
function is_valid_auxtag(x::UInt8, y::UInt8)
    digit = UInt8('0'):UInt8('9')
    upper = UInt8('A'):UInt8('Z')
    lower = UInt8('a'):UInt8('z')
    (in(x, upper) | in(x, lower)) & (in(y, digit) | in(y, upper) | in(y, lower))
end

function try_auxtag(a::UInt8, b::UInt8)::Union{Nothing, AuxTag}
    is_valid_auxtag(a, b) ? AuxTag(unsafe, a, b) : nothing
end

unsafe_u8(c::Char) = (reinterpret(UInt32, c) >> 24) % UInt8

# Convert a char to the corresponding UInt8 if ASCII, else to 0xff
# This is more efficient than the Base implementation.
function to_uint8(x::Char)::UInt8
    u = reinterpret(UInt32, x)
    iszero(u & 0x80ffffff) ? unsafe_u8(x) : 0xff
end

function try_auxtag(x::Char, y::Char)
    # If the chars are not bytes, to_uint8 will result in 0xff,
    # making the constructor error. 
    try_auxtag(to_uint8(x), to_uint8(y))
end

function try_auxtag(x::Union{String, SubString{String}})
    ncodeunits(x) == 2 || return nothing
    cu = codeunits(x)
    cu1 = @inbounds cu[1]
    cu2 = @inbounds cu[2]
    try_auxtag(cu1, cu2)
end

function try_auxtag(x::AbstractString)
    (x, s) = @something iterate(x) return nothing
    (y, s) = @something iterate(x, s) return nothing
    isnothing(iterate(x, s)) || return nothing
    try_auxtag(Char(x)::Char, Char(y)::Char)
end

# We don't want convert to allocate (since it's called implicitly), so
# we only accept string types that can efficiently be converted to tags
function Base.convert(::Type{AuxTag}, x::Union{String, SubString{String}})
    AuxTag(x)
end

Base.print(io::IO, x::AuxTag) = write(io, x.x[1], x.x[2])
Base.show(io::IO, x::AuxTag) = write(io, "AuxTag(\"", x.x[1], x.x[2], "\")")
