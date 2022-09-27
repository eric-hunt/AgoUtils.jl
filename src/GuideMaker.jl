module GuideMaker

export Guide, GuideDNA, GuideRNA, makeguides

include("types.jl")

import BioSequences

using BioSequences: @dna_str, @rna_str
using BioSequences: DNA_A, DNA_C, DNA_G, DNA_T, RNA_A, RNA_C, RNA_G, RNA_U

const DNA_NUCS = BioSequences.DNA[DNA_A, DNA_C, DNA_G, DNA_T]
const RNA_NUCS = BioSequences.RNA[RNA_A, RNA_C, RNA_G, RNA_U]

"""
    makeslices(
        sequence::S, window::Integer; 
        include_overhang = true
    ) where {S <: BioSequences.LongSequence}

Make slices of a sequence using a sliding window.
"""
function _makeslices(
    sequence::S, window::Integer; 
    include_overhang::Bool = true
) where {S <: BioSequences.LongSequence}
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
    target::BioSequences.LongSequence, size::Integer;
    GuideType::Type{G} = BioSequences.LongDNA{4}, include_overhang::Bool = true
) where {G <: BioSequences.LongSequence}

Dispatch on traits
"""
function makeguides(
    target::BioSequences.LongSequence, size::Integer;
    GuideType::Type{G} = BioSequences.LongDNA{4}, include_overhang::Bool = true
) where {G <: BioSequences.LongSequence}
    return makeguides(GuideTrait(G), target, size, GuideType; include_overhang)
end # function makeguides

"""
    makeguides(
    ::IsDNA, target::BioSequences.LongSequence, size::Integer, ::Type{G};
    include_overhang::Bool
) where {G <: BioSequences.LongSequence}

Create DNA guides
"""
function makeguides(
    ::IsDNA, target::BioSequences.LongSequence, size::Integer, ::Type{G};
    include_overhang::Bool
) where {G <: BioSequences.LongSequence}
    sliceseqs = _makeslices(target, size; include_overhang = include_overhang)
    guideseqs = convert.(G, sliceseqs)

    guides = Vector{GuideDNA}()
    for i in eachindex(guideseqs)
        seq = guideseqs[i]
        first = seq[begin]
        altnucs = setdiff(DNA_NUCS, [first])
        altseqs = Vector{G}(undef, length(altnucs))
        for n in eachindex(altnucs)
            altseqs[n] = pushfirst!(seq[begin+1:end], altnucs[n])
        end
        push!(guides, GuideDNA(seq, first, altnucs, altseqs))
    end
    return guides
end # function makeguides (IsDNA)

"""
    makeguides(
    ::IsRNA, target::BioSequences.LongSequence, size::Integer, ::Type{G};
    include_overhang::Bool
) where {G <: BioSequences.LongSequence}

Create RNA guides
"""
function makeguides(
    ::IsRNA, target::BioSequences.LongSequence, size::Integer, ::Type{G};
    include_overhang::Bool
) where {G <: BioSequences.LongSequence}
    sliceseqs = _makeslices(target, size; include_overhang = include_overhang)
    guideseqs = convert.(G, sliceseqs)

    guides = Vector{GuideRNA}()
    for i in eachindex(guideseqs)
        seq = guideseqs[i]
        first = seq[begin]
        altnucs = setdiff(RNA_NUCS, [first])
        altseqs = Vector{G}(undef, length(altnucs))
        for n in eachindex(altnucs)
            altseqs[n] = pushfirst!(seq[begin+1:end], altnucs[n])
        end
        push!(guides, GuideRNA(seq, first, altnucs, altseqs))
    end
    return guides
end # function makeguides (IsRNA)

end # module
