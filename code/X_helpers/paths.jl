# Repository paths. Every helper creates the directory it returns.

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

"Folder of the raw per-tumour files of one cohort: `data/raw/<block>/<cohort>`."
rawdir(block::AbstractString, cohort::AbstractString) = mkpath(joinpath(ROOT, "data", "raw", block, cohort))

"Folder of the processed tables of one block: `data/processed/<block>`."
procdir(block::AbstractString) = mkpath(joinpath(ROOT, "data", "processed", block))

"Folder of the panels of one figure: `figures/<name>`."
figdir(name::AbstractString) = mkpath(joinpath(ROOT, "figures", name))
