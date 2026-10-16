# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.2.0
### Added
- `try_get(aux, key)`, which returns the value of `key`, `nothing` if the key is absent,
  or an `Error` instead of throwing on malformed data.
- `is_well_formed` and `isvalid` take a `check_duplicates` keyword (default `true`).
  Auxiliaries with duplicate keys are now considered malformed.
- `get(f, aux, key)`, `pop!(aux, key)` and `pop!(aux, key, default)` for `Auxiliary`.
- `AuxException` is now public and documented.
- `setindex_nonexisting!` now has a method for `BAM.Auxiliary`. Previously, BAM
  defined its own, distinct function of the same name.
- XAMAuxData can now be compiled with `juliac --trim`.

### Changed
- **Breaking:** Integers are no longer converted to another integer type when written:
  - BAM only accepts `Int8`, `UInt8`, `Int16`, `UInt16`, `Int32` and `UInt32`.
    Other integers, such as `Int64` and `Bool`, throw an `ArgumentError` instead of being
    converted to `Int32`, `UInt32` or `Float32`.
  - SAM accepts Base bit integers in the range `[-2^31, 2^32)`, the range BAM can represent.
    Previously, only the `Int32` range was accepted. `Bool`, `BigInt` and other integer
    types throw an `ArgumentError`.
  - Integer arrays must have one of the six element types above, in both formats.
    Previously, other element types were converted to `Int32` or `UInt32`.
- **Breaking:** SAM values are parsed strictly according to the SAM specs:
  - Integers must match `[-+]?[0-9]+`, so e.g. surrounding whitespace and `0x10` are invalid.
  - Floats must match `[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?`, and be representable as
    a finite `Float32`.
  - Hex (`H`) values must be uppercase.
- **Breaking:** NaN and infinite floats are invalid in both SAM and BAM, including in arrays,
  and cannot be written. `InvalidFloat` can now also occur in BAM.
- **Breaking:** `H` values are loaded as `Hex{Memory{UInt8}}` instead of `Memory{UInt8}`,
  such that writing a loaded value back, or copying between SAM and BAM, preserves the type.
- **Breaking:** `Hex` is now an `AbstractVector{UInt8}`. Use `parent(hex)` to get the
  wrapped vector. Empty `Hex` values can now be written.
- **Breaking:** `length` and `haskey` throw an `AuxException` on malformed data,
  like `get` does, instead of silently returning a result.
- **Breaking:** BAM iteration yields `key => Error` for invalid values, like SAM,
  instead of a bare `Error`.
- **Breaking:** SAM `i` values are loaded as `Int64` on all platforms, instead of `Int`.
- `copy!` from SAM to BAM stores integers as `Int32` if possible, else `UInt32`,
  instead of throwing an `InexactError` for values above `typemax(Int32)`.
- `copy!` between SAM and BAM leaves the destination unchanged if copying fails, and throws
  an `AuxException` on invalid values in the source, instead of a `MethodError`.
- `setindex!` on a malformed `BAM.Auxiliary` throws an `AuxException` instead of an
  `ErrorException`.
- Mutating an `Auxiliary` whose underlying vector has been shrunk below its start index
  now throws an error instead of accessing memory out of bounds.
- Faster parsing of SAM integers and hex values.

### Fixed
- `haskey` looped forever on some malformed SAM data.
- BAM `setindex!` corrupted data when shrinking a value in an `Auxiliary` with a start
  index other than 1.
- SAM `setindex!` deleted the existing field when the new value was rejected.
- Writing a value loaded from the same `Auxiliary`, as in `aux["AB"] = aux["CD"]`,
  could silently corrupt the written value.
- BAM iteration yielding bare `Error`s broke `collect`, `Dict`, `show` and `copy!` on
  invalid auxiliaries.

[Unreleased]: https://github.com/BioJulia/XAMAuxData.jl/compare/v0.1.7...HEAD
