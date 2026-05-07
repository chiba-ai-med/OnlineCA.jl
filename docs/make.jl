using Documenter
using OnlineCA

makedocs(
    sitename = "OnlineCA.jl",
    modules = [OnlineCA],
    format = Documenter.HTML(prettyurls = true),
    pages = [
        "Home" => "index.md",
        "Julia API" => "juliaapi.md",
        "Command line Tool" => "commandline.md"
    ])

deploydocs(
    repo = "github.com/chiba-ai-med/OnlineCA.jl.git",
    devbranch = "main",
    target = "build",
    deps = nothing,
    make = nothing)
