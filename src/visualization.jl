function RecipesBase.plot(pl::PerturbedLatticeModel, Window; arrow=true)
     if pl.grid.d == 2
        grid_x = [p[1] for p in pl.grid.points]
        grid_y = [p[2] for p in pl.grid.points]
        points_x = [p[1] for p in pl.points]
        points_y = [p[2] for p in pl.points]

        p = scatter(grid_x, grid_y, color=:red, markersize=3, label="Grid")

        # Add arrows from grid to perturbed points
        if arrow
            for i in 1:length(pl.points)
                quiver!([pl.grid.points[i][1]], [pl.grid.points[i][2]],
                       quiver=([pl.points[i][1] - pl.grid.points[i][1]],
                              [pl.points[i][2] - pl.grid.points[i][2]]),
                       color=:gray, alpha=0.3)
            end
        end

        scatter!(points_x, points_y,
                color=:blue, markersize=3, label="Perturbed")
        xlims!(Window[1,1], Window[1,2])
        ylims!(Window[2,1], Window[2,2])
        xlabel!("x")
        ylabel!("y")
        plot!(aspect_ratio=:equal)

        return p
    elseif pl.grid.d == 3
        grid_x = [p[1] for p in pl.grid.points]
        grid_y = [p[2] for p in pl.grid.points]
        grid_z = [p[3] for p in pl.grid.points]
        points_x = [p[1] for p in pl.points]
        points_y = [p[2] for p in pl.points]
        points_z = [p[3] for p in pl.points]

        p = plot3d()

        # Connection lines
        for i in 1:length(pl.points)
            plot3d!([pl.grid.points[i][1], pl.points[i][1]],
                   [pl.grid.points[i][2], pl.points[i][2]],
                   [pl.grid.points[i][3], pl.points[i][3]],
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
        plot!(aspect_ratio=:equal)

        return p
    end
end