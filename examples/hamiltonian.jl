using PerturbedLattices

h = StraussHamiltonian(1.0, 1.3, Grid(20, 2))

hc = HardCoreHamiltonian(0.5, Grid(20, 2))

local_energy(h, 1)