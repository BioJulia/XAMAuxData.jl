```@meta
CurrentModule = XAMAuxData
DocTestSetup = quote
    using XAMAuxData: BAM, SAM, AuxTag, Hex, Errors, Error, is_well_formed
end
```

# Reference
```@docs
AuxTag
Hex
Error
Errors
SAM.Auxiliary
BAM.Auxiliary
is_well_formed
Base.isvalid(::XAMAuxData.AbstractAuxiliary)
```
