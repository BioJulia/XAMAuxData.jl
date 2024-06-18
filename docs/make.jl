using Documenter, XAMAuxData

meta = quote
    using XAMAuxData: SAM, BAM
    line = "AK:z:some string\ts1:i:2512\tst:A:+\tas:f:211.2\tar:B:c3,-16,21,-100"
end

DocMeta.setdocmeta!(
    XAMAuxData,
    :DocTestSetup,
    meta,
    recursive=true
)

makedocs(
    sitename = "XAMAuxData.jl",
    modules = [XAMAuxData],
    pages = [
        "Home" => "index.md",
    ],
    authors = "Jakob Nybo Nissen",
    checkdocs = :public,
    remotes=nothing, # TODO: Remove
)

# TODO: Call deploydocs
