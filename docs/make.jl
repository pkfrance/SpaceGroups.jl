using Documenter, SpaceGroups

makedocs(
    sitename = "SpaceGroups",
    modules = [SpaceGroups],
    pages = [
        "Home" => "index.md",
        "API Reference" => "api/public.md",
        "Internal API" => "api/private.md"
    ]
)

deploydocs(
    repo = "github.com/pkfrance/SpaceGroups.jl.git"
)