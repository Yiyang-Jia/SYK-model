using Plots
include("/Users/yiyang/Desktop/SYKcodes/JuliaCodes/RMTfunctions.jl")

t=0:0.1:12
sffplots = plot(t,[sffGUE, sffGOE, sffGSE], labels = ["GUE" "GOE" "GSE"], title="Spectral form factors")
savefig(sffplots,"sffplots.png")
