"""
Visualization functions for perturbed lattice models.
"""

"""
    plot_points(pl::PerturbedLatticeModel; title_str="")

Visualize the current point cloud in 2D or 3D.
"""
function plot_points(pl::PerturbedLatticeModel; title_str="")
    if pl.d == 2
        x_coords = [p[1] for p in pl.points]
        y_coords = [p[2] for p in pl.points]
        scatter(x_coords, y_coords,
                markersize=2,
                markerstrokewidth=0,
                xlabel="x",
                ylabel="y",
                title=title_str,
                aspect_ratio=:equal)
    elseif pl.d == 3
        x_coords = [p[1] for p in pl.points]
        y_coords = [p[2] for p in pl.points]
        z_coords = [p[3] for p in pl.points]
        scatter(x_coords, y_coords, z_coords,
                markersize=2,
                markerstrokewidth=0,
                xlabel="x",
                ylabel="y",
                zlabel="z",
                title=title_str)
    end
end

"""
    plot_point_grid_connection(pl::PerturbedLatticeModel, Window::Matrix{Float64})

Visualize connections between grid points and their perturbed positions.

Shows:
- Red points: Original grid positions
- Blue points: Perturbed positions
- Gray lines/arrows: Connections showing the perturbation
"""
function plot_point_grid_connection(pl::PerturbedLatticeModel, Window::Matrix{Float64})
    if pl.d == 2
        grid_x = [p[1] for p in pl.grid]
        grid_y = [p[2] for p in pl.grid]
        points_x = [p[1] for p in pl.points]
        points_y = [p[2] for p in pl.points]

        p = scatter(grid_x, grid_y, color=:red, markersize=3, label="Grid")

        # Add arrows from grid to perturbed points
        for i in 1:length(pl.points)
            quiver!([pl.grid[i][1]], [pl.grid[i][2]],
                   quiver=([pl.points[i][1] - pl.grid[i][1]],
                          [pl.points[i][2] - pl.grid[i][2]]),
                   color=:gray, alpha=0.3)
        end

        scatter!(points_x, points_y,
                color=:blue, markersize=3, label="Perturbed")
        xlims!(Window[1,1], Window[1,2])
        ylims!(Window[2,1], Window[2,2])
        xlabel!("x")
        ylabel!("y")
        plot!(aspect_ratio=:equal)

        return p
    elseif pl.d == 3
        grid_x = [p[1] for p in pl.grid]
        grid_y = [p[2] for p in pl.grid]
        grid_z = [p[3] for p in pl.grid]
        points_x = [p[1] for p in pl.points]
        points_y = [p[2] for p in pl.points]
        points_z = [p[3] for p in pl.points]

        p = plot3d()

        # Connection lines
        for i in 1:length(pl.points)
            plot3d!([pl.grid[i][1], pl.points[i][1]],
                   [pl.grid[i][2], pl.points[i][2]],
                   [pl.grid[i][3], pl.points[i][3]],
                   color=:gray, alpha=0.3, label="", linewidth=0.5)
        end

        # Grid points
        scatter3d!(grid_x, grid_y, grid_z,
                  color=:red, label="Grid", markersize=2)

        # Perturbed points
        scatter3d!(points_x, points_y, points_z,
                  color=:blue, label="Perturbed", markersize=2)
        xlims!(Window[1,1], Window[1,2])
        ylims!(Window[2,1], Window[2,2])
        zlims!(Window[3,1], Window[3,2])
        xlabel!("x")
        ylabel!("y")
        zlabel!("z")

        return p
    end
end


