"""
Printing
"""


function Base.show(io::IO, ::MIME"text/plain", g::NucleicAcidGuide)
    println(io, summary(g.seq), ':')
    BioSequences.showcompact(io, g.seq)
    println(io, "\n", g.gc * 100, "% GC")
end # function Base.show for NucleicAcidGuide

function Base.show(io::IO, ::MIME"text/plain", s::Vector{NucleicAcidGuide})
    println(io, "$(length(s)) x $(summary(s[1].seq))s:")
    foreach(s) do g
        println(io, g.seq)
    end
end # function Base.show for Vector{NucleicAcidGuide}