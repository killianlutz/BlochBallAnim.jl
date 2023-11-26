function sphereplot!(; fig=current_figure(), axis=current_axis(), radius=1/sqrt(2), npointsθ=20, npointsϕ=15)
    θ = range(0, 2π - 1e-6, npointsθ)
    ϕ = range(0, π - 1e-6, npointsϕ)

    x = radius * cos.(θ)       * sin.(ϕ)'
    y = radius * sin.(θ)       * sin.(ϕ)'
    z = radius * ones(npointsθ) * cos.(ϕ)'

    latitude_points = (x, y, z) .|> vec
    longitude_points = (x, y, z) .|> transpose .|> vec

    latitude_points, longitude_points = map((latitude_points, longitude_points)) do points
        [Point3(X, Y, Z) for (X, Y, Z) in zip(points...)]
    end
    equator = [Point3f(radius*cos(t), radius*sin(t), 0) for t in range(0, 2π - 1e-6, 2*npointsθ)]
    sphere_xaxis = ([Point3(0, 0, 0)], [Point3(1.1radius, 0, 0)])
    sphere_yaxis = ([Point3(0, 0, 0)], [Point3(0, 1.1radius, 0)])
    sphere_zaxis = ([Point3(0, 0, 0)], [Point3(0, 0, 1.1radius)])

    lines!(axis, equator, color=:black)
    arrows!(axis, sphere_xaxis..., color=:black, arrowsize=0.05, linewidth=0.01, alpha=0.5)
    arrows!(axis, sphere_yaxis..., color=:black, arrowsize=0.05, linewidth=0.01, alpha=0.5)
    arrows!(axis, sphere_zaxis..., color=:black, arrowsize=0.05, linewidth=0.01, alpha=0.5)

    scat_origin = scatter!(axis, Point3(0.0, 0.0, 0.0), color=:black, markersize=10, marker=:diamond)
    scatlines_latitude = scatterlines!(axis, latitude_points, color=(:purple, 0.2), markersize=0.1, marker=:circle)
    scatlines_longitude = lines!(axis, longitude_points, color=(:purple, 0.2))

    return (scatlines_latitude, scatlines_longitude, scat_origin)
end

function ball_mesh(; radius=1/sqrt(2), npointsθ=8, npointsϕ=8, npoints_radius=5)
    r = range(0.1, radius, npoints_radius)
    ## more uniform distribution of points wrt. radius
    θs_given_r = map(r) do ρ
        nθ = Int( round( 1/(1.1 - ρ)*npointsθ ) )
        range(0, 2π - 1e-6, nθ)
    end
    ϕs_given_r = map(r) do ρ
        nϕ = Int( round( 0.5/(1.1 - ρ)*npointsϕ ) )
        range(0, π - 1e-6, nϕ)
    end
    x = [ρ * cos.(θ) * sin.(ϕ)' for (ρ, θ, ϕ) in zip(r, θs_given_r, ϕs_given_r)]
    y = [ρ * sin.(θ) * sin.(ϕ)' for (ρ, θ, ϕ) in zip(r, θs_given_r, ϕs_given_r)]
    z = [ρ * ones(length(θ)) * cos.(ϕ)' for (ρ, θ, ϕ) in zip(r, θs_given_r, ϕs_given_r)]

    points = map(x, y, z) do X, Y, Z # loop over radius size
        (X, Y, Z) .|> vec
    end

    points = map(points) do fixed_radius_points # again
        [Point3(u, v, w) for (u, v, w) in zip(fixed_radius_points...)]
    end

    return cat(points..., dims=1) # cat over all radii
end