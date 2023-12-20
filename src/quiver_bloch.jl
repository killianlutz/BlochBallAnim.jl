function bloch_velocity(points, bloch_velocity_parameters)
    scalar_type = eltype(first(points))
    temp_vector = zeros(scalar_type, 3)

    directions = similar(points)
    velocity_magnitudes = similar(directions, scalar_type)
    for (i, base) in enumerate(points)
        bloch_velocity!(temp_vector, base, bloch_velocity_parameters)
        velocity_magnitudes[i] = norm(temp_vector)
        directions[i] = temp_vector/velocity_magnitudes[i]
    end

    return (; directions, velocity_magnitudes)
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