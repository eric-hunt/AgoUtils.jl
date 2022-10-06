"""
Types
"""

using BioSequences: LongSequence, LongNuc, LongDNA, LongRNA,
    NucleicAcidAlphabet, DNAAlphabet, RNAAlphabet
using BioSymbols: BioSymbols.NucleicAcid, BioSymbols.DNA, BioSymbols.RNA


"""
Concrete guide types
"""
struct NucleicAcidGuide{A<:NucleicAcidAlphabet}
    alphabet::Type{A}
    seq::LongSequence{A}
    length::Integer
    firstbase::NucleicAcid
    altbases::Vector{NucleicAcid}
    altseqs::Vector{LongSequence{A}}
end # struct Guide

"An alias for NucleicAcidGuide{A <: NucleicAcidAlphabet{N}}"
const Guide{A} = NucleicAcidGuide{A}

"An alias for NucleicAcidGuide{DNAAlphabet{N}}"
const GuideDNA{N} = NucleicAcidGuide{DNAAlphabet{N}}

"An alias for NucleicAcidGuide{RNAAlphabet{N}}"
const GuideRNA{N} = NucleicAcidGuide{RNAAlphabet{N}}

"""
    NucleicAcidGuide(seq::LongSequence{A}) where {A<:NucleicAcidAlphabet}

Construct a new NucleicAcidGuide from a sequence.
"""
function NucleicAcidGuide(seq::LongSequence{A}) where {A<:NucleicAcidAlphabet}
    first = seq[begin]
    altnucs = setdiff(_altbases(A), [first])
    altseqs = Vector{LongSequence{A}}(undef, length(altnucs))
    for n in eachindex(altnucs)
        altseqs[n] = pushfirst!(seq[begin+1:end], altnucs[n])
    end
    NucleicAcidGuide{A}(
        A,
        seq,
        seq.len,
        first,
        altnucs,
        altseqs
    )
end

Guide(seq) = NucleicAcidGuide(seq)
