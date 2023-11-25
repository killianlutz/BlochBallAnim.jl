function normalized_pauli(C::DataType)
    npauli = [
        [zero(C) 1; 1 0], 
        [zero(C) -im; im 0],
        [one(C) 0; 0 -1]
    ] ./ C(sqrt(2))

    return npauli
end

function blochvector_to_vectorized_density(X::Real, Y::Real, Z::Real, vectorized_npauli::AbstractVector)
    ϱ0 = zero(first(vectorized_npauli))
    components = [X, Y, Z]
    length = norm(components)

    if length > 1 # project bloch ball radius 1
        components ./= length
    end
    for i in 1:3
        ϱ0 .+= components[i] .* vectorized_npauli[i]
    end

    return ϱ0
end

function vecgate(U::AbstractArray{<:Number})
    # maps a quantum gate U acting on qudits to the corresponding 
    # gate vecU = (conj(U) ⊗ U) acting on density operators ϱ
    kron(conj(U), U)
end

function gate1b_hadamard(C::DataType, M::Union{DataType,UnionAll})
    return M([one(C) 1; 1 -1] / C(sqrt(2)))
end

function gate1b_rx(θ::Real, C::DataType, M::Union{DataType,UnionAll})
    t = C(0.5 * θ)
    return M(
        [cos(t) -im*sin(t);
         -im*sin(t) cos(t)]
         )
end

function gate1b_rz(θ::Real,  C::DataType, M::Union{DataType,UnionAll})
    t = C(0.5 * θ)
    return M([exp(-im * t) 0; 0 exp(im * t)])
end

function gates(C::DataType)
    rx = vecgate(gate1b_rx(π/2, C, Matrix{C}))
    rz = vecgate(gate1b_rz(3π/2, C, Matrix{C}))
    hadamard = vecgate(gate1b_hadamard(C, Matrix{C}))

    return (; rx, hadamard, rz)
end


" Builds Lindbladians, i.e. matrices (super-operators) associated to GKSL equation parameters "
function lindbladians(H, Hc, h, γ)
    dim = size(H, 1)
    dim2 = dim^2
    Lf  = oftype(H, spzeros(eltype(H), dim2, dim2)) # free component of GKSL generator
    tmp = zero(H)
    eye = one(H)

    for (hk, γk) in zip(h, γ)
        tmp .+= -γk * hk' * hk             # \sum_k    γk hk† * hk
        Lf .+= kron(2 * γk * conj(hk), hk) # \sum_k 2 *γk hk†^T ⊗ hk
    end
    Lf .+= kron(eye, -im * H + tmp) .+ kron(transpose(im * H) + transpose(tmp), eye)

    Lc = map(Hc) do A
        kron(eye, -im * A) .+ kron(transpose(im * A), eye)
    end # controlled component of GKSL generator

    return (; free=Lf, controls=Lc)
end

function gksl!(dx, x, p, t)
    lindbladian_matrices, control, time_horizon = p

    mul!(dx, lindbladian_matrices.free, x, time_horizon, 0)
    for (control_lindbladian, control_component) in zip(lindbladian_matrices.controls, control)
        mul!(dx, control_lindbladian, x, time_horizon*control_component, 1)
    end

    nothing
end

function simulate_gksl(control::AbstractVector, lindbladian_matrices::NamedTuple, x0::AbstractVecOrMat, time_horizon::Real; tspan=(0.0, 1.0), kwargs...)
    parameters = (lindbladian_matrices, control, time_horizon)
    prob = ODEProblem(gksl!, x0, tspan, parameters)
    solve(prob, Tsit5(); kwargs...)
end

function to_blochball(density_operator::AbstractVecOrMat, npauli::AbstractVector)
    map(density_operator) do ϱ
        map(npauli) do P
            real(tr(reshape(ϱ, 2, 2) * P'))
        end |> Point3
    end
end