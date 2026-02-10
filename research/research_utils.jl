"""
research_utils.jl

Shared helpers for scripts under the `research/` folder.

Goals
- Never write files unless the user explicitly asks for it.
- Keep all generated output under `research/output/...` in a predictable layout.
- Keep scripts easy to run on different machines (paths via environment variables).

Usage pattern (inside a research script)

    include(joinpath(@__DIR__, "..", "research_utils.jl"))

    save_output = get_bool_env("SAVE_OUTPUT", false)
    data_dir = get(ENV, "RESEARCH_DATA_DIR", "")

    out_png = output_path("figure.png"; script_dir=@__DIR__, script_file=@__FILE__)
    if save_output
        ensure_parent_dir(out_png)
        Plots.savefig(p, out_png)
    end

Environment variables
- SAVE_OUTPUT: set to 1/true/yes to enable output saving (default: off)
- RESEARCH_DATA_DIR: optional base directory for input datasets

Notes
- This file intentionally has a friendly tone because it is meant for humans.
- Unicode is allowed in comments/strings in research code.
"""

"""Parse a boolean environment variable in a forgiving way."""
function get_bool_env(name::AbstractString, default::Bool=false)::Bool
    raw = get(ENV, name, "")
    isempty(raw) && return default

    v = lowercase(strip(raw))
    if v in ("1", "true", "t", "yes", "y", "on")
        return true
    elseif v in ("0", "false", "f", "no", "n", "off")
        return false
    end

    @warn "Unrecognized boolean value for ENV[$name]='$raw'. Using default=$default."
    return default
end

"""Absolute path to the `research/` folder (this file lives there)."""
research_root() = @__DIR__

"""Absolute path to `research/output`."""
output_root() = joinpath(research_root(), "output")

"""Return a stable folder name for a script file."""
function _script_stem(script_file::AbstractString)::String
    base = basename(script_file)
    stem, _ = splitext(base)
    return stem
end

"""
Return `script_dir` relative to `research/`.

If `script_dir` is outside `research/`, we fall back to `"misc"`.
"""
function _relative_script_dir(script_dir::AbstractString)::String
    rr = normpath(research_root())
    sd = normpath(script_dir)

    if startswith(sd, rr)
        rel = relpath(sd, rr)
        return rel == "." ? "" : rel
    end

    return "misc"
end

"""Make sure the parent directory of `path` exists."""
function ensure_parent_dir(path::AbstractString)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
    return nothing
end

"""
Compute an output path under `research/output/`.

Layout
- research/output/<relative_script_dir>/<script_stem>/<filename>

We only use the basename of `filename` (so callers can pass "foo/bar.png" safely).
"""
function output_path(filename::AbstractString; script_dir::AbstractString, script_file::AbstractString)::String
    rel_dir = _relative_script_dir(script_dir)
    stem = _script_stem(script_file)

    out_dir = isempty(rel_dir) ? joinpath(output_root(), stem) : joinpath(output_root(), rel_dir, stem)
    return joinpath(out_dir, basename(filename))
end

"""Helpful error if a script expects RESEARCH_DATA_DIR but it is not set."""
function require_data_dir(data_dir::AbstractString; var_name::AbstractString="RESEARCH_DATA_DIR")
    isempty(data_dir) || return nothing

    msg = "This script expects input data files. Please set ENV[$var_name] to the base data directory. " *
          "Example: RESEARCH_DATA_DIR=/path/to/data SAVE_OUTPUT=0 julia <script>.jl"
    throw(ArgumentError(msg))
end
