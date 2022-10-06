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
