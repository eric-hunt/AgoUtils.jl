module AgoUtils

export
    # Types
    NucleicAcidGuide, Guide, GuideDNA, GuideRNA,
    # Methods
    makeguides,
    compileguides,
    # From other package namespaces needed to use AgoUtils.jl
    DNAAlphabet,
    RNAAlphabet,
    @dna_str,
    @rna_str

include("types.jl")
include("utils.jl")
include("printing.jl")

using BioSequences: BioSequences, LongSequence, LongNuc, NucleicAcidAlphabet,
    @dna_str, @rna_str
using BioSymbols: BioSymbols.NucleicAcid, BioSymbols.DNA, BioSymbols.RNA,
    DNA_A, DNA_C, DNA_G, DNA_T, RNA_A, RNA_C, RNA_G, RNA_U


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
    return compiled
end # function compileguides

end # module
