# AGENTS.md

This file provides guidance to LLM agents when working with code in this repository.

## Project Overview

XAMAuxData.jl is a Julia package for parsing and manipulating auxiliary (optional) data fields in SAM/BAM/PAF/GFA bioinformatics file formats. It is part of the BioJulia ecosystem and intended for use by other packages.

Requires Julia >= 1.11.

## Commands

**Run tests:**
```bash
JULIA_TEST_FAILFAST=true julia --project=. --startup=no -e 'using Pkg; Pkg.test()'
```

**Format code (uses Runic v1, not JuliaFormatter):**
```bash
julia -e 'using Runic; Runic.main(["--inplace", "src/", "test/"])'
```
Formatting is enforced in CI via `fredrikekre/runic-action@v1`.

**Build docs:**
```bash
julia --project=docs --startup=no docs/make.jl
```

## Architecture

The package provides two parallel submodules with nearly identical APIs:

- **`SAM`** (`src/sam.jl`) — text-based tab-delimited format (`TAG:TYPE:VALUE`)
- **`BAM`** (`src/bam.jl`) — binary format with type-specific byte encodings

Both expose an `Auxiliary` type that implements `AbstractDict{AuxTag, Any}` with lazy parsing. The mutable variant (`MutableAuxiliary`) wraps a `Vector{UInt8}` and supports mutation; the immutable variant wraps any `AbstractVector{UInt8}`.

**Key types defined in the main module** (`src/XAMAuxData.jl`):
- `AbstractAuxiliary{T} <: AbstractDict{AuxTag, Any}` — shared abstract type
- `AuxTag` (`src/auxtag.jl`) — 2-byte immutable key, validated against `[A-Za-z][A-Za-z0-9]`
- `Hex` — wrapper for hex-encoded byte arrays (type tag `H`)
- `Error` (enum in `Errors` module) — represents data corruption in values

**Data corruption model:** The package distinguishes "malformed" (can't identify key/value boundaries; checked with `is_well_formed()`) from "invalid" (boundaries found but values corrupt; returns `Error` enum values instead of throwing; checked with `isvalid()`).

**Cross-format conversion:** `Base.copy!` at the bottom of `XAMAuxData.jl` converts between SAM and BAM auxiliary data.

## Code Conventions

- Format with Runic (see `.JuliaFormatter.toml` for settings)
- Uses `MemoryViews` for zero-copy operations and `StringViews` for string views
- `@inbounds` used in performance-critical inner loops
- Internal sentinel type `Unsafe`/`unsafe` for bypassing validation in constructors
- Errors returned as values (`Error` enum) rather than thrown, except for malformed data
