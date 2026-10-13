```@meta
CurrentModule = XAMAuxData
DocTestSetup = quote
    using XAMAuxData: BAM, SAM, AuxTag, Hex, Errors, Error, is_well_formed, setindex_nonexisting!, try_get
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
try_get
Base.isvalid(::XAMAuxData.AbstractAuxiliary)
setindex_nonexisting!
```
