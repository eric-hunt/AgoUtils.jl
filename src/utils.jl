"""
Utilities
"""

using BioSequences: BioSequences, LongNuc, DNAAlphabet, RNAAlphabet


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
) where {S<:LongNuc}
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
    _altbases(::A) where {A <: NucleicAcidAlphabet}

Get alternate bases for Alphabet
"""
_altbases(::Type{<:DNAAlphabet}) = BioSequences.symbols(DNAAlphabet{2}())
_altbases(::Type{<:RNAAlphabet}) = BioSequences.symbols(RNAAlphabet{2}())


"""
    _fetchseqs(guide::NucleicAcidGuide)

Get all sequences for a single guide.
"""
function _fetchseqs(guide::NucleicAcidGuide)
    seqs = [guide.seq, guide.altseqs...]
    return seqs
end # function _fetchseqs


const WELL_LETTERS = 'A':'P'
const WELL_NUMBERS = 1:24

_wells_96 = reduce(
    vcat,
    map(n -> map(l -> l * string(n), WELL_LETTERS[1:8]), WELL_NUMBERS[1:12])
)

_wells_384 = reduce(
    vcat,
    vcat(
        map(n -> map(l -> l * string(n), WELL_LETTERS[1:8]), WELL_NUMBERS[1:12]),
        map(n -> map(l -> l * string(n), WELL_LETTERS[9:16]), WELL_NUMBERS[1:12]),
        map(n -> map(l -> l * string(n), WELL_LETTERS[1:8]), WELL_NUMBERS[13:24]),
        map(n -> map(l -> l * string(n), WELL_LETTERS[9:16]), WELL_NUMBERS[13:24]),
    )
)
