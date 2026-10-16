# AGENTS.md
XAMAuxData.jl is a Julia package for parsing and manipulating auxiliary (optional) data fields in SAM/BAM/PAF/GFA bioinformatics file formats. It is part of the BioJulia ecosystem and intended for use by other packages.

## Design criteria
The package distinguishes "malformed" (can't identify key/value boundaries; checked with `is_well_formed()`) from "invalid" (boundaries found but values corrupt; returns `Error` enum values instead of throwing; checked with `isvalid()`).

## Code Conventions
- When making edits, format with `runic -i .`. If not installed, alert the user that the call failed,
  and do not attempt to install Runic.
- Uses `MemoryViews` for zero-copy operations and `StringViews` for string views
- `@inbounds` used in performance-critical inner loops
- Errors returned as values (`Error` enum) rather than thrown, except for malformed data

## Other guidelines
Do not add yourself (the AI agent) as a co-author to git commits
