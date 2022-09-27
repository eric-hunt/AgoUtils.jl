"""
Types
"""

import BioSequences

"""
Abstract guide type
"""
abstract type Guide end # abstract Guide

"""
Concrete guide types
"""
struct GuideDNA <: Guide
    seq::Union{BioSequences.LongDNA, Nothing}
    firstbase::Union{BioSequences.DNA, Nothing}
    altbases::Union{Vector{BioSequences.DNA}, Nothing}
    altseqs::Union{Vector{BioSequences.LongDNA}, Nothing}
end # struct GuideDNA

struct GuideRNA <: Guide
    seq::Union{BioSequences.LongRNA, Nothing}
    firstbase::Union{BioSequences.RNA, Nothing}
    altbases::Union{Vector{BioSequences.RNA}, Nothing}
    altseqs::Union{Vector{BioSequences.LongRNA}, Nothing}
end # struct GuideRNA

"""
Traits
"""
abstract type GuideTrait end
struct IsDNA <: GuideTrait end
struct IsRNA <: GuideTrait end

GuideTrait(::Type) = IsDNA() # default
GuideTrait(::Type{<: BioSequences.LongDNA}) = IsDNA()
GuideTrait(::Type{<: BioSequences.LongRNA}) = IsRNA()
