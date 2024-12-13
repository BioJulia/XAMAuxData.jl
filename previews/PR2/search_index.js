var documenterSearchIndex = {"docs":
[{"location":"reference/","page":"Reference","title":"Reference","text":"CurrentModule = XAMAuxData\nDocTestSetup = quote\n    using XAMAuxData: BAM, SAM, AuxTag, Hex, Errors, Error, is_well_formed\nend","category":"page"},{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"AuxTag\nHex\nError\nErrors\nSAM.Auxiliary\nBAM.Auxiliary\nis_well_formed\nBase.isvalid(::XAMAuxData.AbstractAuxiliary)","category":"page"},{"location":"reference/#XAMAuxData.AuxTag","page":"Reference","title":"XAMAuxData.AuxTag","text":"AuxTag(::Union{String, SubString{String}})\nAuxTag(::Char, ::Char)\nAuxTag(::UInt8, ::UInt8)\n\nThe type of a key of Auxiliary. This type is lightweight, and validates that the key is of appropriate format. Auxiliary can also be modified using strings as keys, which will be converted to AuxTag, but this may cause a needless allocation of the string.\n\nExamples\n\njulia> AuxTag(\"AC\")\nAuxTag(\"AC\")\n\njulia> AuxTag('k', 'w')\nAuxTag(\"kw\")\n\njulia> AuxTag(\"1C\") # invalid tag\nERROR: Invalid AuxTag. Tags must conform to r\"^[A-Za-z][A-Za-z0-9]$\".\n[...]\n\n\n\n\n\n","category":"type"},{"location":"reference/#XAMAuxData.Hex","page":"Reference","title":"XAMAuxData.Hex","text":"Hex(v::AbstractVector{UInt8})\n\nWrapper type around a byte vector. When this type is assigned to an Auxiliary, it is encoded as a hex value (type tag H) as opposed to a raw byte array ( type tag B).\n\njulia> aux = SAM.Auxiliary(UInt8[], 1);\n\njulia> aux[\"AB\"] = Hex([0xae, 0xf8, 0x6c]);\n\njulia> String(MemoryView(aux))\n\"AB:H:AEF86C\"\n\njulia> aux = BAM.Auxiliary(UInt8[], 1);\n\njulia> aux[\"AB\"] = Hex([0xae, 0xf8, 0x6c]);\n\njulia> String(MemoryView(aux))\n\"ABHAEF86C\\0\"\n\n\n\n\n\n","category":"type"},{"location":"reference/#XAMAuxData.Errors.Error","page":"Reference","title":"XAMAuxData.Errors.Error","text":"Error\n\nEnum type representing errors returned when loading invalid auxiliary data values. The errors are nonexhausitve - more might be added in a non-breaking release.\n\nThe following errors may contained in AbstractAuxiliary instead of the real value:\n\nInvalidTypeTag (SAM only): An aux value with an unknown type\nInvalidArrayEltype (SAM only): A B value with an unknown element type\nInvalidInt (SAM only): An integer that can't be parsed to an Int32\nInvalidFloat (SAM only): A float that can't be parsed to a Float32\nInvalidChar: Loading a Char not in '!':'~'\nInvalidString: A string that contains a character not in re\"[ !-~]\"\nInvalidHex: a H value with an odd number of symbols, or symbol not in re\"[0-9A-F]\"\nInvalidArray: A malformed B value\n\nSee also: Errors\n\n\n\n\n\n","category":"type"},{"location":"reference/#XAMAuxData.Errors","page":"Reference","title":"XAMAuxData.Errors","text":"Errors\n\nModule containing the error values of XAMAuxData.\n\nSee also: Error\n\n\n\n\n\n","category":"module"},{"location":"reference/#XAMAuxData.SAM.Auxiliary","page":"Reference","title":"XAMAuxData.SAM.Auxiliary","text":"SAM.Auxiliary{T <: AbstractVector{UInt8}} <: AbstractDict{AuxTag, Any}\n\nLazily loaded AbstractDict representing the auxiliary data fields of a SAM record. Immutable aux's can be constructed with Auxiliary(x) for any x with MemoryView(x) defined. Mutable aux data is constructed with Auxiliary(x::Vector{UInt8}, start::Int), where start gives the first index of the used data in x - all data before start will be ignored and never modified.\n\nExamples\n\njulia> immut = SAM.Auxiliary(\"KJ:i:-1\tAB:Z:abc\");\n\njulia> immut[\"KJ\"]\n-1\n\njulia> haskey(immut, \"AB\")\ntrue\n\nExtended help\n\nSince fields of Auxiliary are lazily loaded, auxiliaries may contain invalid data even after successful construction. Auxiliary operates with two distinct notions of invalid data - malformedness and invalidness. See is_well_formed for their definitions.\n\n\n\n\n\n","category":"type"},{"location":"reference/#XAMAuxData.BAM.Auxiliary","page":"Reference","title":"XAMAuxData.BAM.Auxiliary","text":"BAM.Auxiliary{T <: AbstractVector{UInt8}} <: AbstractDict{AuxTag, Any}\n\nLazily loaded AbstractDict representing the auxiliary data fields of a BAM record. Immutable aux's can be constructed with Auxiliary(x) for any x with MemoryView(x) defined. Mutable aux data is constructed with Auxiliary(x::Vector{UInt8}, start::Int), where start gives the first index of the used data in x - all data before start will be ignored and never modified.\n\nExamples\n\njulia> immut = BAM.Auxiliary(\"KJS\u0007\\0ABZabc\\0\");\n\njulia> immut[\"KJ\"]\n0x0007\n\njulia> haskey(immut, \"AB\")\ntrue\n\nExtended help\n\nExtended help\n\nSince fields of Auxiliary are lazily loaded, auxiliaries may contain invalid data even after successful construction. Auxiliary operates with two distinct notions of invalid data - malformedness and invalidness. See is_well_formed for their definitions. ```\n\n\n\n\n\n","category":"type"},{"location":"reference/#XAMAuxData.is_well_formed","page":"Reference","title":"XAMAuxData.is_well_formed","text":"is_well_formed(aux::AbstractAuxiliary) -> Bool\n\nCheck if the auxiliary is well formed. AbstractAuxiliary distinguishes two ways their data can be incorrect: If the data is corrupted such that it is not possible to identify the key, or start/end of the encoded data of the value for a key/value pair, we say the auxiliary is malformed. Accessing values, or iterating such an auxiliary may or may not throw an exception. This notion of malformedness is what this function checks for.\n\nAlternatively, if the keys and encoded values can be identified, but the encoded values are corrupt and the value cannot be loaded, we say the auxiliary is invalid. Such auxiliaries can be iterated over, and the values can be loaded, although the loaded value will be of type XAMAuxData.Error. This can be checked for with isvalid. Valid auxiliaries are always well-formed.\n\nExamples\n\njulia> aux = SAM.Auxiliary(\"KM:i:252\\tAK::C\"); # bad key\n\njulia> is_well_formed(aux)\nfalse\n\njulia> aux[\"KM\"] # loading a value may error\n252\n\njulia> aux[\"AK\"] # loading a value may error\nERROR: Invalid SAM tag header. Expected <AuxTag>:<type tag>:, but found no colons.\n[...]\n\njulia> aux = SAM.Auxiliary(\"AB:Z:αβγδ\\tCD:f:-1.2\"); # bad value encoding\n\njulia> (is_well_formed(aux), isvalid(aux))\n(true, false)\n\njulia> aux[\"AB\"] # note: numerical value of enum is not API\nInvalidString::Error = 9\n\n\n\n\n\n","category":"function"},{"location":"reference/#Base.isvalid-Tuple{XAMAuxData.AbstractAuxiliary}","page":"Reference","title":"Base.isvalid","text":"isvalid(aux::AbstractAuxiliary) -> Bool\n\nCheck if the auxiliary is well formed, and that all the values can be loaded without returning an XAMAuxData.Error.\n\nExamples\n\njulia> aux = SAM.Auxiliary(\"AB:i:not an integer\");\n\njulia> is_well_formed(aux)\ntrue\n\njulia> isvalid(aux)\nfalse\n\nSee also: is_well_formed\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = XAMAuxData\nDocTestSetup = quote\n    using XAMAuxData: BAM, SAM, AuxTag, Hex, Errors, Error, is_well_formed\n    using MemoryViews: MemoryView\nend","category":"page"},{"location":"#XAMAuxData.jl","page":"Home","title":"XAMAuxData.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package contains functionality to parse the auxiliary (optional) fields in the SAM and BAM format. The formats PAF and GFA (and possibly others) share this same mini-format of SAM aux fields. This package is intended to be used by other packages, such as XAM.jl and other parsing packages.","category":"page"},{"location":"#Core-concepts","page":"Home","title":"Core concepts","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"XAMAuxData has two submodules - SAM and BAM. SAM is used for the text-based auxiliary format in SAM, PAF and GFA files. BAM is used for the binary encoded aux format in BAM files. Most examples in this documentation will use the SAM format, since it's more human readable. Any differences to the BAM format will be explicitly mentioned.","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Note\nAnnoyingly, the specification of GFA auxiliary fields differ slightly from that of SAM auxiliary fields. Currently, this package implements only the SAM specification, and as such does not fully support GFA files. In the future, a dedicated GFA module may be introduced.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The single auxiliary field AN:i:1234 is encoded as the key-value pair AuxTag(\"AN\") => 1234. A collection of aux fields are represented by a SAM.Auxiliary (or BAM.Auxiliary), which are subtypes of AbstractDict{AuxTag, Any}.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The package may be used like this:","category":"page"},{"location":"","page":"Home","title":"Home","text":"# Import the module you want to use\njulia> using XAMAuxData: SAM\n\njulia> data = \"AN:A:z\\ta1:Z:abc def \\tkv:i:-25234\\tzz:f:-14.466e-3\";\n\njulia> aux = SAM.Auxiliary(data)\n4-element XAMAuxData.SAM.Auxiliary{MemoryViews.ImmutableMemoryView{UInt8}}:\n  \"AN\" => 'z'\n  \"a1\" => \"abc def \"\n  \"kv\" => -25234\n  \"zz\" => -0.014466f0","category":"page"},{"location":"#The-AuxTag-object","page":"Home","title":"The AuxTag object","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The keys of an Auxiliary are instances of AuxTag. Why not simply have them be two-byte strings? AuxTags are compact bitstypes and therefore not heap-allocated for extra performance. Also, the construction of an AugTax validates that it conforms to the regex [A-Za-z][A-Za-z0-9] per the SAM specs.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Strings can be converted to AuxTag for convenience, as in below:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> push!(AuxTag[], \"AB\") # implicit conversion\n1-element Vector{AuxTag}:\n AuxTag(\"AB\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"Attempting to construct an invalid AuxTag will error:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> AuxTag(\"11\")\nERROR: Invalid AuxTag. Tags must conform to r\"^[A-Za-z][A-Za-z0-9]$\".\n[...]","category":"page"},{"location":"#Constructing-Auxiliary-objects","page":"Home","title":"Constructing Auxiliary objects","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SAM.Auxiliary and BAM.Auxiliary are constructed the same two ways.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Immutable auxiliaries are constructed from any bytes-like object which has a MemoryView method. This may be a String, SubString{String}, Memory{UInt8} etc. Auxiliary objects are constructed directly from these:","category":"page"},{"location":"","page":"Home","title":"Home","text":"# Make an IMMUTABLE Auxiliary\njulia> aux = SAM.Auxiliary(\"AB:i:12\\tKN:A:z\")\n2-element XAMAuxData.SAM.Auxiliary{MemoryViews.ImmutableMemoryView{UInt8}}:\n  \"AB\" => 12\n  \"KN\" => 'z'","category":"page"},{"location":"","page":"Home","title":"Home","text":"In order to mutate auxiliary objects, their data need to be able to be resized, and therefore must be backed by a Vector{UInt8}. Mutable auxiliaries can be constructed from a Vector{UInt8} and a starting index. Any data before the starting index is never touched by the Auxiliary object, even if e.g. the object is emptied. This starting index enables Auxiliary objects to share the same Vector that is used to store other data of PAF, SAM or GFA records.","category":"page"},{"location":"","page":"Home","title":"Home","text":"In the example below, the first 22 bytes of the vector (the some not-aux data here part) corresponds to the data before the aux data, and hence the first index of the aux data in the vector is 23.","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> data = collect(codeunits(\"some not-aux data hereAB:i:12\\tKN:A:z\"));\n\njulia> aux = SAM.Auxiliary(data, 23) # Make a MUTABLE Auxiliary\n2-element XAMAuxData.SAM.Auxiliary{Vector{UInt8}}:\n  \"AB\" => 12\n  \"KN\" => 'z'","category":"page"},{"location":"","page":"Home","title":"Home","text":"No matter whether constructed from a memory view or from a Vector, there cannot be any unused bytes at or after the starting index in an Auxiliary. Any trailing bytes will be considered part of the auxiliary data, and may possibly be considered invalid:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> bad_aux = SAM.Auxiliary(\"AB:A:p\\t\\t\"); # trailing tabs\n\njulia> isvalid(bad_aux)\nfalse","category":"page"},{"location":"#Manipulating-Auxiliary-objects","page":"Home","title":"Manipulating Auxiliary objects","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Auxiliarys can be read and written like a normal AbstractDict{AuxTag, Any}:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> aux = SAM.Auxiliary(UInt8[], 1); # empty Auxiliary\n\njulia> # Note: The strings are implicitly `convert`ed to AuxTag\n       aux[\"AX\"] = 'y'; aux[\"AX\"]\n'y': ASCII/Unicode U+0079 (category Ll: Letter, lowercase)\n\njulia> aux[\"cm\"] = 12.1; aux[\"cm\"]\n12.1f0\n\njulia> aux[\"G1\"] = [-1.24, 33.1]; aux[\"G1\"]\n2-element Memory{Float32}:\n -1.24\n 33.1\n\njulia> aux[\"AX\"] = [0x01, 0x02]; aux[\"AX\"] # overwrite AX key\n2-element Memory{UInt8}:\n 0x01\n 0x02","category":"page"},{"location":"","page":"Home","title":"Home","text":"Like Dict, the order of key/value pairs in auxiliaries is arbitrary and cannot be relied on.","category":"page"},{"location":"#More-details-about-element-type","page":"Home","title":"More details about element type","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The SAM and BAM formats support only support a limited set of types, and the supported types differ between SAM and BAM. When writing a value to an Auxiliary object, this package will attempt to convert your value to one of the supported types. Hence, the value written to an Auxiliary may not be the same value when being read back.","category":"page"},{"location":"#Reading:-SAM-element-types","page":"Home","title":"Reading: SAM element types","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SAM type Julia type read\ni Int\nA Char\nf Float32\nZ StringView*\nH Memory{UInt8}\nB:C Memory{UInt8}\nB:c Memory{Int8}\nB:S Memory{UInt16}\nB:s Memory{Int8}\nB:I Memory{UInt32}\nB:i Memory{Int32}\nB:f Memory{Float32}","category":"page"},{"location":"","page":"Home","title":"Home","text":"This is a parametric type. The precise concrete type depends on the type the view looks into, and is an implementation detail.","category":"page"},{"location":"#Writing:-SAM-element-types","page":"Home","title":"Writing: SAM element types","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Input type SAM type stored\nInteger* i\nAbstractChar✝ A\nAbstractFloat‡ f\nAbstractString§ Z\nAbstractVector{UInt8} B:C\nAbstractVector{Int8} B:c\nAbstractVector{UInt16} B:S\nAbstractVector{Int16} B:s\nAbstractVector{<:Signed}¶ B:i\nAbstractVector{<:Unsigned}¶ B:I\nAbstractVector{<:AbstractFloat}‡ B:f\nHex H","category":"page"},{"location":"","page":"Home","title":"Home","text":"* Only values in typemin(Int32):typemax(Int32) are allowed. This is because the SAM specs recommend to limit yourselves to this range of values. However, when reading i fields, an Int is returned, such that values outside the recommended range can still be read on 64-bit systems.\n✝ Permitted Char values are only those in '!':'~'.\n‡ Only values representable by a Float32 are allowed.\n§ Only characters in '!':'~' and spaces (' ') are permitted in strings\n¶ These are stored as Int32 and UInt32 for Signed and Unsigned, respectively.","category":"page"},{"location":"#BAM-element-types","page":"Home","title":"BAM element types","text":"","category":"section"},{"location":"#BAM-integers","page":"Home","title":"BAM integers","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The main difference between SAM and BAM types is that the latter format permits different types of integers. Hence, except the types mentioned below, all the SAM types in the table above are also supported in BAM, with the same Julia <-> BAM type correspondance. Further, reading a value of i will return an Int32 instead of an Int.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Input type BAM type Julia type read\nUInt8 C UInt8\nInt8 c Int8\nUInt16 S UInt16\nInt16 s Int16\nUnsigned I UInt32\nSigned i Int32","category":"page"},{"location":"#BAM-Arrays","page":"Home","title":"BAM Arrays","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Whereas reading a SAM array (type B) will result in a Memory, BAM arrays only promise that they return an AbstractVector of the correct element type:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> aux = BAM.Auxiliary(\"ABBs\\2\\0\\0\\0\\1\\2\\3\\4\");\n\njulia> aux[\"AB\"] isa AbstractVector{Int16}\ntrue","category":"page"},{"location":"#The-Hex-type","page":"Home","title":"The Hex type","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The SAM/BAM type H signifies a hex-encoded byte array. Byte arrays are stored as the B:C type when written, since this type is more legible (in SAM), and more space efficient in BAM. Therefore, the B:C type is generally recommended over the H type when writing byte arrays. However, if you want to write an AbstractVector{UInt8} value explicitly as an H type, wrap it in the Hex type:","category":"page"},{"location":"","page":"Home","title":"Home","text":"aux = SAM.Auxiliary(UInt8[], 1)\naux[\"AB\"] = UInt8[0x01, 0x02]\n\nusing MemoryViews\n# Print the memory content of the aux.\n# The array was written as a value of the type B:c\nprintln(String(MemoryView(aux)))\n\n# Wrap input type in the Hex type\naux[\"AB\"] = Hex(UInt8[0x01, 0x02])\n\n# It is now written as a H instead\nprintln(String(MemoryView(aux)))\n\n# output\nAB:B:C,1,2\nAB:H:0102","category":"page"},{"location":"#Writing-Auxiliarys-to-files","page":"Home","title":"Writing Auxiliarys to files","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Calling MemoryView on an Auxiliary will return a view of the underlying data. This data is guaranteed to be valid SAM/BAM auxiliary data:","category":"page"},{"location":"","page":"Home","title":"Home","text":"(field1, field2) = (\"PG:A:n\", \"ca:i:-241\")\naux = SAM.Auxiliary(field1 * '\\t' * field2)\n\n# Get a view of the data underlying `aux`.\n# This is guaranteed to be valid SAM data (and likewise for BAM)\nusing MemoryViews\nmem = MemoryView(aux)\n\n# We make no guarantees about which order the two fields are,\n# but we DO guarantee the memory is a valid SAM aux data\n# with these two fields\nprintln(\n  in(\n    String(mem),\n    (field1 * '\\t' * field2, field2 * '\\t' * field1)\n  )\n)\n\n# output\ntrue\n","category":"page"},{"location":"#Invalid-data-in-Auxiliaries","page":"Home","title":"Invalid data in Auxiliaries","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The elements of an Auxiliary are lazily loaded, and in the interest of speed, there is no mandatory validation of the data done when constructing an Auxiliary. Hence, they may contain invalid data. This package distinguishes two different kinds of bad data:","category":"page"},{"location":"","page":"Home","title":"Home","text":"If the data is malformed in such a way that it's not possible to identify the keys of the auxiliary, or the data segment of the corresponding values, we say the auxiliary is malformed. Loading keys or values from malformed auxiliaries may throw an exception. The function is_well_formed can be used to check for malformed auxiliaries.\nIf the keys and the data segments corresponding to the values can be identified, but the data itself is corrupt such that the values cannot be loaded, we instead say the auxiliary is invalid. Loading an invalid value return an object of type Error. The validity of an auxiliary can be checked with isvalid. All valid records are also well-formed.","category":"page"},{"location":"","page":"Home","title":"Home","text":"For example:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> # This is completely mangled\n\njulia> aux = SAM.Auxiliary(\"erwlifju093\");\n\njulia> (is_well_formed(aux), isvalid(aux))\n(false, false)\n\njulia> # Keys and values can be identified, but data can't be loaded as an integer\n\njulia> aux = SAM.Auxiliary(\"AB:i:dslkjas\");\n\njulia> (is_well_formed(aux), isvalid(aux))\n(true, false)\n\njulia> only(values(aux)) isa Error\ntrue","category":"page"}]
}