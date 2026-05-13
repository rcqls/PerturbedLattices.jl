""" 
Simulation of a perturbed lattice model with Strauss potential and Gaussian perturbations.

"""

using PerturbedLattices

beta_values = [-log(0.9)]
sigma_values = [0.5, 1.0]
RS_values = [0.5, 1.5]
N_lattice = 30
warmup_iterations = 100
Window = [-5.0 5.0; -5.0 5.0]

for beta in beta_values
    for sigma in sigma_values
        for rs in RS_values
            println("="^60)
            println("Simulating lattice with beta = $beta, sigma = $sigma")
            println("="^60)

            # Create the lattice
            pl = PerturbedLattice(N_lattice, RS=rs, beta=beta, sigma=sigma, seed=1234)

            # Warmup phase
            println("Starting warmup ($warmup_iterations iterations)...")
            @time begin
                for iter in 1:warmup_iterations
                    iterate!(pl)
                end
            end
            println("Warmup completed!\n")

            # Plot the point grid connection
            p = plot_point_grid_connection(pl, Window)

            # Save the image
            filename = "lattice_beta_$(round(beta, digits=3))_sigma_$(sigma)_RS_$(rs).png"
            savefig(p, filename)
            println("Saved plot to $filename\n")
        end
    end
end

println("All simulations completed!")