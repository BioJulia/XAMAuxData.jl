```@meta
CurrentModule = XAMAuxData
DocTestSetup = quote
    using XAMAuxData: BAM, SAM, AuxTag, Hex, Errors, Error, is_well_formed
    using MemoryViews: MemoryView
end
```

# XAMAuxData.jl
This package contains functionality to parse the auxiliary (optional) fields in the SAM and BAM format.
The formats PAF and GFA (and possibly others) share this same mini-format of SAM aux fields.
This package is intended to be used by other packages, such as XAM.jl and other parsing packages.

## Core concepts
XAMAuxData has two submodules - `SAM` and `BAM`. `SAM` is used for the text-based auxiliary format in SAM,
PAF and GFA files. `BAM` is used for the binary encoded aux format in BAM files.
Most examples in this documentation will use the SAM format, since it's more human readable.
Any differences to the BAM format will be explicitly mentioned.

!!! note
    Annoyingly, the specification of GFA auxiliary fields differ slightly from that of SAM
    auxiliary fields.
    Currently, this package implements only the SAM specification, and as such does not
    fully support GFA files. In the future, a dedicated GFA module may be introduced.

The single auxiliary field `AN:i:1234` is encoded as the key-value pair `AuxTag("AN") => 1234`.
A collection of aux fields are represented by a `SAM.Auxiliary` (or `BAM.Auxiliary`), which are subtypes of `AbstractDict{AuxTag, Any}`.

The package may be used like this:
```jldoctest
# Import the module you want to use
julia> using XAMAuxData: SAM

julia> data = "AN:A:z\ta1:Z:abc def \tkv:i:-25234\tzz:f:-14.466e-3";

julia> aux = SAM.Auxiliary(data)
4-element XAMAuxData.SAM.Auxiliary{MemoryViews.ImmutableMemoryView{UInt8}}:
  "AN" => 'z'
  "a1" => "abc def "
  "kv" => -25234
  "zz" => -0.014466f0
```

#### The `AuxTag` object
The keys of an `Auxiliary` are instances of `AuxTag`.
Why not simply have them be two-byte strings?
`AuxTags` are compact bitstypes and therefore not heap-allocated for extra performance.
Also, the construction of an `AugTax` validates that it conforms to the regex `[A-Za-z][A-Za-z0-9]` per the SAM specs.

Strings can be converted to `AuxTag` for convenience, as in below:

```jldoctest
julia> push!(AuxTag[], "AB") # implicit conversion
1-element Vector{AuxTag}:
 AuxTag("AB")
```

Attempting to construct an invalid `AuxTag` will error:
```jldoctest
julia> AuxTag("11")
ERROR: Invalid AuxTag. Tags must conform to r"^[A-Za-z][A-Za-z0-9]$".
[...]
```

## Constructing `Auxiliary` objects
`SAM.Auxiliary` and `BAM.Auxiliary` are constructed the same two ways.

Immutable auxiliaries are constructed from any bytes-like object which has a `MemoryView` method.
This may be a `String`, `SubString{String}`, `Memory{UInt8}` etc.
Auxiliary objects are constructed directly from these:

```jldoctest
# Make an IMMUTABLE Auxiliary
julia> aux = SAM.Auxiliary("AB:i:12\tKN:A:z")
2-element XAMAuxData.SAM.Auxiliary{MemoryViews.ImmutableMemoryView{UInt8}}:
  "AB" => 12
  "KN" => 'z'
```

In order to mutate auxiliary objects, their data need to be able to be resized, and therefore must be backed by a `Vector{UInt8}`.
Mutable auxiliaries can be constructed from a `Vector{UInt8}` and a starting index. Any data before the starting index is never touched by the `Auxiliary` object, even if e.g. the object is emptied.
This starting index enables `Auxiliary` objects to share the same `Vector` that is used to store other data of PAF, SAM or GFA records.

In the example below, the first 22 bytes of the vector (the `some not-aux data here` part)
corresponds to the data before the aux data, and hence the first index of the aux data in the vector is 23.

```jldoctest
julia> data = collect(codeunits("some not-aux data hereAB:i:12\tKN:A:z"));

julia> aux = SAM.Auxiliary(data, 23) # Make a MUTABLE Auxiliary
2-element XAMAuxData.SAM.Auxiliary{Vector{UInt8}}:
  "AB" => 12
  "KN" => 'z'
```

No matter whether constructed from a memory view or from a `Vector`, there cannot be any unused bytes at or after the starting index in an `Auxiliary`.
Any trailing bytes will be considered part of the auxiliary data, and may possibly be considered invalid:
```jldoctest
julia> bad_aux = SAM.Auxiliary("AB:A:p\t\t"); # trailing tabs

julia> isvalid(bad_aux)
false
```

## Manipulating `Auxiliary` objects
`Auxiliary`s can be read and written like a normal `AbstractDict{AuxTag, Any}`:

```jldoctest
julia> aux = SAM.Auxiliary(UInt8[], 1); # empty Auxiliary

julia> # Note: The strings are implicitly `convert`ed to AuxTag
       aux["AX"] = 'y'; aux["AX"]
'y': ASCII/Unicode U+0079 (category Ll: Letter, lowercase)

julia> aux["cm"] = 12.1; aux["cm"]
12.1f0

julia> aux["G1"] = [-1.24, 33.1]; aux["G1"]
2-element Memory{Float32}:
 -1.24
 33.1

julia> aux["AX"] = [0x01, 0x02]; aux["AX"] # overwrite AX key
2-element Memory{UInt8}:
 0x01
 0x02
```

Like `Dict`, the order of key/value pairs in auxiliaries is arbitrary and cannot be relied on.

## More details about element type
The SAM and BAM formats support only support a limited set of types, and the supported types differ between SAM and BAM.
When writing a value to an `Auxiliary` object, this package will attempt to convert your value to one of the supported types.
Hence, the value written to an `Auxiliary` may not be the same value when being read back.

### Reading: SAM element types
| SAM type | Julia type read  |
| -------- | ---------------- |
| `i`      | `Int`            |
| `A`      | `Char`           |
| `f`      | `Float32`        |
| `Z`      | `StringView`*    |
| `H`      | `Memory{UInt8}`  |
| `B:C`    | `Memory{UInt8}`  |
| `B:c`    | `Memory{Int8}`   |
| `B:S`    | `Memory{UInt16}` |
| `B:s`    | `Memory{Int8}`   |
| `B:I`    | `Memory{UInt32}` |
| `B:i`    | `Memory{Int32}`  |
| `B:f`    | `Memory{Float32}`|

* This is a parametric type. The precise concrete type depends on the type the view looks into,
  and is an implementation detail.

### Writing: SAM element types
| Input type                       | SAM type stored |
| -------------------------------- | --------------- |
| `Integer`*                        | `i`             |
| `AbstractChar`✝                   | `A`            |
| `AbstractFloat`‡                  | `f`             |
| `AbstractString`§                 | `Z`             |
| `AbstractVector{UInt8}`           | `B:C`           |
| `AbstractVector{Int8}`            | `B:c`           |
| `AbstractVector{UInt16}`          | `B:S`           |
| `AbstractVector{Int16}`           | `B:s`           |
| `AbstractVector{<:Signed}`¶       | `B:i`           |
| `AbstractVector{<:Unsigned}`¶     | `B:I`           |
| `AbstractVector{<:AbstractFloat}`‡| `B:f`           |
| `Hex`                             | `H`             |

- `*` Only values in `typemin(Int32):typemax(Int32)` are allowed.
  This is because the SAM specs recommend to limit yourselves to this range of values.
  However, when reading `i` fields, an `Int` is returned, such that values
  outside the recommended range can still be read on 64-bit systems.
- ✝ Permitted `Char` values are only those in `'!':'~'`.
- ‡ Only values representable by a `Float32` are allowed.
- § Only characters in `'!':'~'` and spaces (`' '`) are permitted in strings
- ¶ These are stored as `Int32` and `UInt32` for `Signed` and `Unsigned`, respectively.

### BAM element types
#### BAM integers
The main difference between SAM and BAM types is that the latter format permits different types of integers.
Hence, except the types mentioned below, all the SAM types in the table above are also supported in BAM,
with the same Julia <-> BAM type correspondance.
Further, reading a value of `i` will return an `Int32` instead of an `Int`.

| Input type | BAM type | Julia type read|
| -----------|----------|--------------- |
| `UInt8`    | `C`      | `UInt8`        |
| `Int8`     | `c`      | `Int8`         |
| `UInt16`   | `S`      | `UInt16`       |
| `Int16`    | `s`      | `Int16`        |
| `Unsigned` | `I`      | `UInt32`       |
| `Signed`   | `i`      | `Int32`        |

#### BAM Arrays
Whereas reading a SAM array (type `B`) will result in a `Memory`, BAM arrays only promise
that they return an `AbstractVector` of the correct element type:

```jldoctest
julia> aux = BAM.Auxiliary("ABBs\2\0\0\0\1\2\3\4");

julia> aux["AB"] isa AbstractVector{Int16}
true
```

### The `Hex` type
The SAM/BAM type `H` signifies a hex-encoded byte array.
Byte arrays are stored as the `B:C` type when written, since this type is more legible (in SAM),
and more space efficient in BAM.
Therefore, the `B:C` type is generally recommended over the `H` type when writing byte arrays.
However, if you want to write an `AbstractVector{UInt8}` value explicitly as an `H` type, wrap it in the `Hex` type:

```jldoctest
aux = SAM.Auxiliary(UInt8[], 1)
aux["AB"] = UInt8[0x01, 0x02]

using MemoryViews
# Print the memory content of the aux.
# The array was written as a value of the type B:c
println(String(MemoryView(aux)))

# Wrap input type in the Hex type
aux["AB"] = Hex(UInt8[0x01, 0x02])

# It is now written as a H instead
println(String(MemoryView(aux)))

# output
AB:B:C,1,2
AB:H:0102
```

## Writing `Auxiliary`s to files
Calling `MemoryView` on an `Auxiliary` will return a view of the underlying data.
This data is guaranteed to be valid SAM/BAM auxiliary data:

```jldoctest
(field1, field2) = ("PG:A:n", "ca:i:-241")
aux = SAM.Auxiliary(field1 * '\t' * field2)

# Get a view of the data underlying `aux`.
# This is guaranteed to be valid SAM data (and likewise for BAM)
using MemoryViews
mem = MemoryView(aux)

# We make no guarantees about which order the two fields are,
# but we DO guarantee the memory is a valid SAM aux data
# with these two fields
println(
  in(
    String(mem),
    (field1 * '\t' * field2, field2 * '\t' * field1)
  )
)

# output
true

```

## Invalid data in `Auxiliaries`
The elements of an `Auxiliary` are lazily loaded, and in the interest of speed,
there is no mandatory validation of the data done when constructing an `Auxiliary`.
Hence, they may contain invalid data.
This package distinguishes two different kinds of bad data:

1. If the data is malformed in such a way that it's not possible to identify
   the keys of the auxiliary, or the data segment of the corresponding values,
   we say the auxiliary is malformed.
   Loading keys or values from malformed auxiliaries _may_ throw an exception.
   The function [`is_well_formed`](@ref) can be used to check for malformed auxiliaries.

2. If the keys and the data segments corresponding to the values _can_ be identified,
   but the data itself is corrupt such that the values cannot be loaded, we instead
   say the auxiliary is invalid.
   Loading an invalid value return an object of type [`Error`](@ref).
   The validity of an auxiliary can be checked with `isvalid`.
   All valid records are also well-formed.

For example:
```jldoctest
julia> # This is completely mangled

julia> aux = SAM.Auxiliary("erwlifju093");

julia> (is_well_formed(aux), isvalid(aux))
(false, false)

julia> # Keys and values can be identified, but data can't be loaded as an integer

julia> aux = SAM.Auxiliary("AB:i:dslkjas");

julia> (is_well_formed(aux), isvalid(aux))
(true, false)

julia> only(values(aux)) isa Error
true
```
