module GuideMaker

export 
    # Types
    NucleicAcidGuide, Guide, GuideDNA, GuideRNA,
    # Methods
    makeguides,
    # From other package namespaces needed to use GuideMaker.jl
    DNAAlphabet,
    RNAAlphabet,
    @dna_str,
    @rna_str

include("types.jl")

using BioSequences: BioSequences, LongSequence, LongNuc, NucleicAcidAlphabet,
    @dna_str, @rna_str
using BioSymbols: BioSymbols.NucleicAcid, BioSymbols.DNA, BioSymbols.RNA,
    DNA_A, DNA_C, DNA_G, DNA_T, RNA_A, RNA_C, RNA_G, RNA_U


const DNA_NUCS = DNA[DNA_A, DNA_C, DNA_G, DNA_T]
const RNA_NUCS = RNA[RNA_A, RNA_C, RNA_G, RNA_U]


"""
    _makeslices(
        sequence::S, window::Integer; 
        include_overhang::Bool = true
    ) where {S <: LongNuc}

Make slices of a sequence using a sliding window.
"""
function _makeslices(
    sequence::S, window::Integer; 
    include_overhang::Bool
) where {S <: LongNuc}
    seqlength = length(sequence)
    ∆ = window - 1 # \increment
    # if sequence::S is (+) strand..
    # create slices from (-) strand in 3' -> 5' orientation
    complement = BioSequences.complement(sequence)
    
    slices = Vector{S}()
    for i = 1:(seqlength-∆)
        push!(slices, complement[range(i, i + ∆)])
    end
    
    if include_overhang
        push!(slices, push!(last(slices)[begin+1:end], 'A'))
    end

    # return (-) strand slices in 5' -> 3' orientation
	reverse!.(slices)
	
    return slices
end # function _makeslices


"""
    makeguides(
        target::LongNuc, size::Integer, ::Type{A};
        include_overhang::Bool
    ) where {A <: NucleicAcidAlphabet}

Create guides for a sequence
"""
function makeguides(
    target::LongNuc, size::Integer, ::Type{A};
    include_overhang::Bool = true
) where {A <: NucleicAcidAlphabet}
    sliceseqs = _makeslices(target, size; include_overhang = include_overhang)
    guideseqs = convert.(LongSequence{A}, sliceseqs)

    guides = Vector{NucleicAcidGuide}()
    for i in eachindex(guideseqs)
        seq = guideseqs[i]
        first = seq[begin]
        altnucs = setdiff(DNA_NUCS, [first])
        altseqs = Vector{LongSequence{A}}(undef, length(altnucs))
        for n in eachindex(altnucs)
            altseqs[n] = pushfirst!(seq[begin+1:end], altnucs[n])
        end
        push!(guides, NucleicAcidGuide{A}(seq, first, altnucs, altseqs))
    end
    return guides
end # function makeguides

end # module
