"""
Printing
"""


function Base.show(io::IO, ::MIME"text/plain", g::NucleicAcidGuide)
    println(io, summary(g.seq), ':')
    BioSequences.showcompact(io, g.seq)
end # function Base.show for NucleicAcidGuide

function Base.show(io::IO, ::MIME"text/plain", s::Vector{NucleicAcidGuide})
    println(io, "$(length(s)) x $(summary(s[1].seq))s:")
    foreach(s) do g
        println(io, g.seq)
    end
end # function Base.show for Vector{NucleicAcidGuide}