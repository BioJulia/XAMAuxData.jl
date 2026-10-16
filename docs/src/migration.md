# Migrating from 0.1 to 0.2
```@meta
CurrentModule = XAMAuxData
```

This page lists the breaking changes in version 0.2, and how to update code written for 0.1.
For a full list of changes, including new features, see the changelog.

## Integers are no longer converted when written
Previously, integers were converted to a type the format could store. Now, only
integer types that can be stored without conversion are accepted, and other types
throw an `ArgumentError`:
* BAM accepts `Int8`, `UInt8`, `Int16`, `UInt16`, `Int32` and `UInt32`.
  In particular, `Int` (`Int64`) is no longer accepted.
* SAM accepts Base bit integers, such as `Int`, in the range `[-2^31, 2^32)`.
  `Bool` and `BigInt` are no longer accepted.
* Integer arrays must have one of the six element types accepted by BAM, in both formats.

To migrate, convert the value to the desired type before writing, e.g.
`aux["AB"] = Int32(x)` instead of `aux["AB"] = x`, or `aux["AB"] = Int32.(v)` for an array.

## `H` values are loaded as `Hex`
`H` values are now loaded as `Hex{Memory{UInt8}}` instead of `Memory{UInt8}`.
`Hex` is now an `AbstractVector{UInt8}`, so most code using the loaded value keeps working,
but generic operations on `Hex` may be slow.
Use `parent(aux["AB"])` to get the `Memory{UInt8}` returned in 0.1.

## `length` and `haskey` throw on malformed data
`length` and `haskey` now throw an [`AuxException`](@ref) when the auxiliary is malformed,
like `get` does. Previously, `length` returned the number of fields before the malformed
part, and `haskey` returned `false` if the key was not found before the malformed part.

To handle malformed data, check [`is_well_formed`](@ref) first.

## BAM iteration yields pairs for invalid values
Iterating a `BAM.Auxiliary` now yields `key => Error` for invalid values, like
`SAM.Auxiliary` does. Previously, it yielded a bare `Error`.
Code checking for errors during iteration should check the value of the pair, e.g.
`for (key, value) in aux; value isa Error && ...; end`.

## Stricter validation of values
The following values are now invalid, and are loaded as an [`Error`](@ref):
* SAM integers, floats and array elements not matching the regexes given in the SAM specs,
  e.g. with surrounding whitespace or in hexadecimal notation.
* SAM floats not representable as a finite `Float32`.
* NaN and infinite floats, in both SAM and BAM, including in arrays.
  These can also no longer be written.
* `H` values with lowercase hex digits, in both SAM and BAM.
