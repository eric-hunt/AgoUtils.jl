module AgoUtils

export
    # Types
    NucleicAcidGuide, Guide, GuideDNA, GuideRNA,
    # Methods
    compileguides,
    exportguides,
    makeguides,
    # From other package namespaces needed to use AgoUtils.jl
    DNAAlphabet,
    RNAAlphabet,
    @dna_str,
    @rna_str

include("types.jl")
include("utils.jl")
include("printing.jl")

using BioSequences: @dna_str, @rna_str, LongSequence, LongNuc,
    NucleicAcidAlphabet, DNAAlphabet, RNAAlphabet, NucleicAcid
using CSV: CSV
using DataFrames: DataFrames, DataFrame
using Dates: Dates


# TODO Add option to select slicing algorithm (see utils.jl _makeslices()).
"""
    makeguides(
        target::LongNuc, size::Integer, ::Type{A};
        include_overhang::Bool
    ) where {A <: NucleicAcidAlphabet}

Create guides for a sequence
"""
function makeguides(
    target::LongNuc, size::Integer, ::Type{A};
    include_overhang::Bool=true
) where {A<:NucleicAcidAlphabet}
    sliceseqs = _makeslices(target, size; include_overhang=include_overhang)
    guideseqs = convert.(LongSequence{A}, sliceseqs)

    guides = Vector{NucleicAcidGuide}()
    for seq in guideseqs
        push!(guides, NucleicAcidGuide(seq))
    end
    return guides
end # function makeguides


"""
    compileguides(guides::Vector{NucleicAcidGuide})

Compile all guide sequences for a substrate in order.
"""
function compileguides(guides::Vector{NucleicAcidGuide})
    compiled = mapreduce(_fetchseqs, vcat, guides)
    guidedict = [(:wells => _wells_384[1:length(compiled)]), (:seqs => compiled)]
    guidetable = DataFrame(guidedict)
    return guidetable
end # function compileguides


"""
    exportguides(
        guides::DataFrame,
        dir::AbstractString = joinpath(homedir(), "Desktop"),
        suffix::AbstractString = ""
    )

Save guide sequences in CSV format for ordering.
"""
function exportguides(
    guides::DataFrame;
    dir::AbstractString = joinpath(homedir(), "Desktop"),
    suffix::AbstractString = ""
)
    dir = expanduser(normpath(dir))
    if !isdir(dir)
        error("The containing directory path is not valid.")
    end
    suffix = replace(suffix, r"\s" => "")
    timestamp = Dates.format(Dates.now(), "yymmdd_HHMMSS")
    CSV.write(
        joinpath(
            dir,
            string(
                "guides_",
                timestamp,
                # I'm not sure what is considered the best style here..
                # !isempty(suffix) ? "_" : "",
                # !isempty(suffix) && "_$suffix",
                if !isempty(suffix) "_$suffix" end,
                ".csv"
            )
        ),
        guides
    )
end # function exportguides

function exportguides(guides::Vector{NucleicAcidGuide}, varargs...)
    guidetable = compileguides(guides)
    exportguides(guidetable, varargs...)
end # function exportguides

end # module
