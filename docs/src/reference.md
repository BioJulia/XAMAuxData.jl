```@meta
CurrentModule = XAMAuxData
DocTestSetup = quote
    using XAMAuxData: BAM, SAM, AuxTag, Hex, Errors, Error
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
Base.isvalid(::XAMAuxData.AbstractAuxiliary)
```