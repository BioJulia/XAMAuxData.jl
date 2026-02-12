```@meta
CurrentModule = XAMAuxData
DocTestSetup = quote
    using XAMAuxData: BAM, SAM, AuxTag, Hex, Errors, Error, is_well_formed, setindex_nonexisting!
end
```

# Reference
```@docs
AuxTag
try_auxtag
Hex
Error
Errors
SAM.Auxiliary
BAM.Auxiliary
is_well_formed
Base.isvalid(::XAMAuxData.AbstractAuxiliary)
setindex_nonexisting!
```
