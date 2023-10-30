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

function bloch_velocity(points, bloch_velocity_parameters)
    scalar_type = eltype(first(points))
    temp_vector = zeros(scalar_type, 3)

    directions = similar(points)
    inward_components = similar(directions, scalar_type)
    for (i, base) in enumerate(points)
        bloch_velocity!(temp_vector, base, bloch_velocity_parameters)
        inward_components[i] = norm(temp_vector)
        directions[i] = temp_vector/inward_components[i]
    end

    return (; directions, inward_components)
end

function bloch_velocity!(dx, x, p)
    velocity_field, control, time_horizon = p

    mul!(dx, velocity_field.free, x, time_horizon, 0)
    for (control_bloch_lindbladian, control_component) in zip(velocity_field.controls, control)
        mul!(dx, control_bloch_lindbladian, x, time_horizon*control_component, 1)
    end
    axpy!(time_horizon, velocity_field.drift, dx)

    return Point3(dx)
end

function to_bloch_matrix(A::AbstractMatrix, basis_vectors::AbstractVector)
    return [real(dot(y, A, x)) for y in basis_vectors, x in basis_vectors]
end

function to_bloch_drift(A::AbstractMatrix, basis_vectors::AbstractVector)
    return [real(dot(y, A, vec(I(2))))/2 for y in basis_vectors]
end

function to_bloch_matrix(lindbladian_matrices::NamedTuple, npauli::AbstractVector)
    vec_npauli = map(vec, npauli)

    bloch_free = to_bloch_matrix(lindbladian_matrices.free, vec_npauli)
    bloch_controls = map(lindbladian_matrices.controls) do A
        to_bloch_matrix(A, vec_npauli)
    end
    drift = to_bloch_drift(lindbladian_matrices.free, vec_npauli)

    return (; free=bloch_free, controls=bloch_controls, drift)
end