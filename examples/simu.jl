using PerturbedLattices
#using Plots


grid =  Grid(20, 2)
h = StraussHamiltonian(1.0, 1.3,grid)
move = GaussianMoveModel([0.5 0.; 0. 0.5], 2)

# Create the lattice
pl = PerturbedLatticeModel(h, move, grid)

# Warmup phase
println("Starting warmup ...")
@time begin
    rand(pl)
end
println("Warmup completed!\n")

# Plot the point grid connection
#p = plot_point_grid_connection(pl, Window)

display(p)


println("All simulations completed!")