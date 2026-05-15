""" 
Estimation procedure for Strauss interaction and gaussian interaction parameters based on Takacs-Fiksel method.
"""

using PerturbedLattices
using Statistics

liste_RS = [1.5, 0.5]
liste_beta = [-log(0.9),-log(0.5), -log(0.1)]
liste_sigma = [0.5, 1.0]
liste_n = [8., 12., 16.]

beta_init = 0.5
theta_init = 0.5

N_est = 30
B = [-5. 5.; -5. 5.]
base_seed = 1234

# Construire toutes les combinaisons de paramètres
params = [(r, b, s, n) for r in liste_RS for b in liste_beta for s in liste_sigma for n in liste_n]
n_params = length(params)

# Pré-allouer les résultats (chaque indice est accédé par un seul thread)
res_array = Vector{Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64, Float64}}(undef, n_params)
full_beta_local = [Vector{Tuple{Float64, Float64, Float64, Float64, Float64}}() for _ in 1:n_params]
full_theta_local = [Vector{Tuple{Float64, Float64, Float64, Float64, Float64}}() for _ in 1:n_params]

println("Lancement de $n_params tâches sur $(Threads.nthreads()) threads")

Threads.@threads for idx in 1:n_params
    r, b, s, n = params[idx]
    println("[Thread $(Threads.threadid())] Début: RS=$r, beta=$b, sigma=$s, n=$n")

    plS = PerturbedLattice(30, RS=r, beta=b, sigma=s, seed=base_seed + idx)

    # Warm-up
    for iter in 1:1_000_000
        iterate!(plS)
    end

    m_n = log(n)
    W_n = [-n n; -n n]

    beta_hat = Float64[]
    theta_hat = Float64[]
    for n_simu in 1:10
        for iter in 1:10_000
            iterate!(plS)
        end
        points_D_n, grid_D_n, points_W_n_boundary = points_in_window(plS, W_n, m_n)
        beta_opt, theta_opt, dlr_opt, vec_loc_en_W, config_adjacency = fit(plS, beta_init, theta_init, points_D_n, grid_D_n, points_W_n_boundary, W_n, N_est, B)
        push!(beta_hat, beta_opt)
        push!(theta_hat, theta_opt)
        push!(full_beta_local[idx], (beta_opt, r, b, s, n))
        push!(full_theta_local[idx], (theta_opt, r, b, s, n))
    end

    res_array[idx] = (mean(theta_hat), std(theta_hat), (mean(theta_hat) - 1/(s^2))^2 + var(theta_hat), mean(beta_hat), std(beta_hat), (mean(beta_hat) - b)^2 + var(beta_hat), r, b, s, n)
    println("[Thread $(Threads.threadid())] Terminé: RS=$r, beta=$b, sigma=$s, n=$n")
end

# Rassembler les résultats
res = collect(res_array)
full_beta = vcat(full_beta_local...)
full_theta = vcat(full_theta_local...)


using CSV, DataFrames
df = DataFrame(theta_bar_1 = [r[1] for r in res],
                   sigma_theta_1 = [r[2] for r in res],
                   mse_1 = [r[3] for r in res],
                   theta_bar_2 = [r[4] for r in res],
                   sigma_theta_2 = [r[5] for r in res],
                   mse_2 = [r[6] for r in res],
                   RS = [r[7] for r in res],
                   beta = [r[8] for r in res],
                   sigma = [r[9] for r in res],
                   n = [r[10] for r in res])
CSV.write("output_final.csv", df, delim=';')

df1 = DataFrame(beta_opt = [r[1] for r in full_beta],
                    RS = [r[2] for r in full_beta],
                    beta_true = [r[3] for r in full_beta],
                    sigma = [r[4] for r in full_beta],
                    n = [r[5] for r in full_beta])
CSV.write("full_beta_final.csv", df1, delim=';')

df2 = DataFrame(theta_opt = [r[1] for r in full_theta],
                    RS = [r[2] for r in full_theta],
                    theta_true = [r[3] for r in full_theta],
                    sigma = [r[4] for r in full_theta],
                    n = [r[5] for r in full_theta])
CSV.write("full_theta_final.csv", df2, delim=';')