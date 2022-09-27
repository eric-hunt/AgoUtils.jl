module GuideMaker

export Guide, GuideDNA, GuideRNA, makeguides

include("types.jl")

import BioSequences

using BioSequences: @dna_str, @rna_str

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

end # module
