""" Runs the Bloch ball (in fact: COHERENCE ball) animation in an external window using GLMakie. \n
Available kwargs: \n 
``H0<:AbstractMatrix{<:Number}'': free Hamiltonian \n
``h<:Vector{<:AbstractMatrix{<:Number}}'': jump_operators."""
function bloch_animation(; kwargs...)
    fig, axis_bloch = setup_figure(; location=[1, 2]) 

    control_spoint = setup_selectpoint!(; fig, location=[1, 3], sublocation=[1, 1])
    slider_grid = setup_slidergrid!(; fig, location=[1, 3], sublocation=[2, 1])   
    toggles = setup_toggles!(; fig, location=[1, 3], sublocation=[3, 1])
    menu = setup_menu!(; fig, location=[1, 3], sublocation=[4, 1])
    
    orbit_obs, quiver_obs = liftobservables(toggles, control_spoint, slider_grid, menu; kwargs...)
    
    vectorplots = blochvectorplot!(orbit_obs; fig, axis=axis_bloch)
    quiverplots = blochquiverplot!(quiver_obs, vectorplots; fig, axis=axis_bloch)
    hide_and_show!(toggles.quiver, vectorplots, quiverplots)

    colgap!(fig.layout, 1, Relative(1/20))
    resize_to_layout!(fig)

    display(fig)
    return fig
end


#### SOURCE 
function setup_figure(; location=[1, 1])
    radius = 1/sqrt(2)
    fig = Figure(size=(1000, 770))
    axis_bloch = Axis3(
        fig[location...], 
        title="Bloch ball",
        limits=(-radius, radius, -radius, radius, -radius, radius),
        xzpanelcolor = (:black, 0.05), 
        yzpanelcolor = (:black, 0.05),
        xypanelcolor = (:black, 0.05),
        xlabel=L"X",
        ylabel=L"Y",
        zlabel=L"Z",
        aspect=:equal
    )
    hidedecorations!(axis_bloch, label=false)
    axis_bloch.azimuth = π/4

    sphereplot!(; fig, axis=axis_bloch, radius)
    return (fig, axis_bloch)
end

function blochquiverplot!(quiver_observables, vectorplots; fig=current_figure(), axis=current_axis(), isvisible=false)
    points, directions, velocity_magnitudes, colorrange = quiver_observables
    
    arrow_object = arrows!(axis, points, directions, 
        color=velocity_magnitudes, 
        lengthscale=0.03f0, 
        arrowsize=0.03f0*Vec3f(1.f0, 1, 1),
        linewidth=0.015f0, 
        colormap=:viridis,
        colorrange=colorrange,
        visible=isvisible
    )
    
    Colorbar(fig[1, 1][1, 1],
        label="purity",
        vectorplots.lines_orbitbloch, 
        tellwidth=true, 
        tellheight=true,
        ticks=[0.5, 1.0],
        labelpadding=-25.0
    )
    Colorbar(fig[1, 1][2, 1],
        label="velocity magnitude",
        arrow_object, 
        tellwidth=true, 
        tellheight=true,
        ticklabelsvisible=false,
        ticksvisible=false,
        labelpadding=0.0
    )

    return (; arrow_object, )
end

function blochvectorplot!(orbit_observables; fig=current_figure(), axis=current_axis())
    orbit, orbit_color, orbit_begin, orbit_end, gate_endstate, control_inbloch = orbit_observables
    arrow_base = [Point3f(0.f0, 0.f0, 0.f0)]

    lines_orbitbloch = lines!(axis, orbit, color=orbit_color, colormap=c=cgrad(:jet, rev=true), colorrange=(0.5, 1.0), linewidth=6)
    scatter_beginbloch = scatter!(axis, orbit_begin, marker=:star5, markersize=30, color=(:orange, 0.95))
    scatter_endbloch = scatter!(axis, orbit_end, marker=:xcross, markersize=30, color=(:blue, 0.95))
    scatter_endgate = scatter!(axis, gate_endstate, marker=:cross, markersize=30, color=(:green, 0.95))
    arrows!(axis, arrow_base, control_inbloch; color=:red, arrowsize=0.08, linewidth=0.02)
    
    Legend(fig[1, 1][3, 1], 
        [lines_orbitbloch, scatter_beginbloch, scatter_endbloch, scatter_endgate], 
        [L"ϱ(t)", L"ϱ(0)", L"ϱ(T)", L"\mathcal{Q}ϱ(0)"],
        orientation=:vertical,
        labelsize=16
    )

    return (; lines_orbitbloch, scatter_beginbloch, scatter_endbloch, scatter_endgate)
end

" When quiver plot: orbit is off and quiver is on; else: vice-versa "
function hide_and_show!(quiver_toggle, vectorplots, quiverplots)
    on(quiver_toggle.active) do isvisible
        for object in quiverplots
            object.visible = isvisible
        end
        for object in vectorplots
            object.visible = !(isvisible)
        end
    end
end

function liftobservables(toggles, control_spoint, slider_grid, menu; kwargs...)
    # pre-allocations
    C = ComplexF32 # faster since Makie converts to Float32 anyways
    saveat = range(0.f0, 1.f0, 500)
    qgates = gates(C)
    npauli = normalized_pauli(C)
    vectorized_npauli = vec.(npauli)# ./ C(sqrt(2)))
    pairings = map(menu) do m 
        options = to_value(m.options)
        Dict(s => i for (s, i) in zip(options, 1:3))
    end

    ###### ORBIT OBSERVABLES
    odep = gksl_odeparameters(toggles.gksl; kwargs...)
    lindbladian_matrices = lift(odep.free_hamiltonian, odep.control_hamiltonians, odep.damping_rates) do x, y, z
        lindbladians(x, y, odep.jump_operators, z)        
    end

    control_effective = lift(control_spoint, slider_grid.control_amplitude) do direction, amplitude
        isapprox(direction, zero(direction)) ? zero(direction) : amplitude.*(direction./norm(direction))
    end

    # ϱ0 density matrix
    initial_density = lift(slider_grid.Xcomponent_ϱ0, slider_grid.Ycomponent_ϱ0, slider_grid.Zcomponent_ϱ0) do X, Y, Z
        blochvector_to_vectorized_density(X, Y, Z, vectorized_npauli)
    end
    # ϱ(t) GKSL solution
    orbit_density = lift(control_effective, lindbladian_matrices, initial_density, slider_grid.time_horizon) do control_eff, lindbladian_mat, ϱ0, time_horizon
        simulate_gksl(control_eff, lindbladian_mat, ϱ0, time_horizon; saveat).u
    end
    # x(t) coherence vectors associated to ϱ(t)
    orbit_bloch = lift(orbit_density) do orbit
        to_blochball(orbit, npauli)
    end

    orbit_color = lift(orbit_bloch) do orbit; 0.5f0 .+ norm.(orbit).^2; end # purity
    orbit_begin = lift(first, orbit_bloch)
    orbit_end = lift(last, orbit_bloch)

    # target quantum state
    gate_endstate = lift(initial_density, menu.gate.selection) do ϱ0, gate_selec
        gates_pairing = first(pairings)
        to_blochball(
            [reshape(qgates[gates_pairing[gate_selec]] * ϱ0, 2, 2)],
            npauli
        )
    end

    # equatorial angular velocity Ω
    control_inbloch = lift(control_spoint) do control
        if isapprox(control, zero(control))
            return [Point3f(first(control), last(control), 0.f0)]
        else
            c = maximum(abs.(control)) * 3
            return [Point3f(first(control)/c, last(control)/c, 0.f0)]
        end
    end

    ###### QUIVER OBSERVABLES
    quiver_points = ball_mesh() # location of the tip of velocity vectors
    scalar_type = eltype(first(quiver_points))

    directions = Observable(zero(quiver_points))
    velocity_magnitudes = Observable( zeros(scalar_type, length(quiver_points)) )
    colorrange = Observable( (zero(scalar_type), one(scalar_type)) )
    
    onany(toggles.quiver.active, lindbladian_matrices, control_effective, slider_grid.time_horizon) do isquiver_active, L, u, T
        if isquiver_active
            velocity_field = to_bloch_matrix(L, npauli)
            bvelocity = bloch_velocity(quiver_points, (velocity_field, u, T))
            minmax = extrema(bvelocity.velocity_magnitudes)

            directions[] = bvelocity.directions
            velocity_magnitudes[] = bvelocity.velocity_magnitudes # FOR NOW ITS THE LENGTH
            if isapprox(first(minmax), last(minmax))
                colorrange[] = (zero(scalar_type), one(scalar_type))
            else
                colorrange[] = minmax
            end
        else
            nothing
        end
    end

    ######## USEFUL OBSERVABLES
    orbit_observables = (;
        orbit_bloch, 
        orbit_color, 
        orbit_begin, 
        orbit_end, 
        gate_endstate, 
        control_inbloch
    )
    quiver_observables = (;
        quiver_points, 
        directions,
        velocity_magnitudes,
        colorrange
    )

    return (orbit_observables, quiver_observables)
end

function gksl_odeparameters(toggles; H0::M=[-1 2; 2 1], h::V=[[0 1; 0 0], [0 0; 1 0]]) where {M<:AbstractMatrix{<:Number}, V<:Vector{<:AbstractMatrix{<:Number}}}
    C = ComplexF32
    if !(ishermitian(H0))
        throw(ArgumentError("Matrix H0 must be HERMITIAN, i.e. H0 = adjoint(H0)"))
    end
    # DIVIDE CONTROL MATRICES BY 2 => INTERPRET CONTROL AS ANGULAR VELOCITY LEGITIMATE!
    H_c = Matrix{C}.([[0 1; 1 0]/2, [0 -im; im  0]/2])
    H_0 = Matrix{C}(H0)
    jump_operators = Matrix{C}.(h)
    γ = ones(real(C), length(h))

    free_hamiltonian = lift(x -> x ? H_0 : zero(H_0), toggles.free.active)
    control_hamiltonians = lift(x -> x ? H_c : zero.(H_c), toggles.control.active)
    damping_rates = lift(x -> x ? γ : zero(γ), toggles.noisy.active)

    return (; free_hamiltonian, control_hamiltonians, jump_operators, damping_rates)
end