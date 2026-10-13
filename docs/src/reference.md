```@meta
CurrentModule = XAMAuxData
DocTestSetup = quote
    using XAMAuxData: BAM, SAM, AuxTag, Hex, Errors, Error, AuxException, is_well_formed, setindex_nonexisting!, try_get
end
```

# Reference
```@docs
AuxTag
try_auxtag
Hex
Error
Errors
AuxException
SAM.Auxiliary
BAM.Auxiliary
is_well_formed
try_get
Base.isvalid(::XAMAuxData.AbstractAuxiliary)
setindex_nonexisting!
```
