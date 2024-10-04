# SYK-model
Codes for studying the SYK model.

MATLABcodes folder contains codes for constructing the SYK Hamiltonian and doing the exact diagonalization.

spectral-form-factor.nb: 

    1. import data from exact diagonalization of SYK Hamiltonian (raw data generated from the ''MATLABcodes'' folder in this  repository),
    
    2. calcualte average adjacent level spacing ratio,
    
    3. unfold the spectra (slow, use Julia codes)
    
    4. calculate and plot the unfolded and connected spectral form factor(slow, use Julia codes). 
    
JuliaCodes folder: do all the spectral form factor computation in spectral-form-factor.nb (and more) but much faster.

SYKmoments.nb: calculates the single-trace moments (up to Tr H^8) and double-trace moments (up to TrH^6 TrH^6) of SYK. Algorithms follow arXiv:1801.02696 and arXiv:1912.11923


double-trace-chord-diagram-classification: 

    1. find equivalent chord diagrams due to cyclic permutation and reflection of ordering of gamaa matrices inside the two traces,
    
    2. different classes chord diagrams can still have the same intersection structures after this classification, that I worked out by hands.
