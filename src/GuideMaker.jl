module GuideMaker

export Guide, GuideDNA, GuideRNA

include("types.jl")

import BioSequences

using BioSequences: @dna_str, @rna_str

"""
    makeslices(
        sequence::T, window::Integer; 
        include_overhang = true
    ) where {T <: BioSequences.LongSequence}

Make slices of a sequence using a sliding window.
"""
function _makeslices(
    sequence::T, window::Integer; 
    include_overhang = true
) where {T <: BioSequences.LongSequence}
	
    seqlength = length(sequence)
    ∆ = window - 1 # \increment
    # if sequence::T is (+) strand..
    # create slices from (-) strand in 3' -> 5' orientation
    complement = BioSequences.complement(sequence)
    
    slices = Vector{T}()
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

end # module
