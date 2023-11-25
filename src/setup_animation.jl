""" Runs the Bloch ball (in fact: COHERENCE ball) animation in an external window using GLMakie. \n
Available kwargs: \n 
``H0<:AbstractMatrix{<:Number}'': free Hamiltonian \n
``h<:Vector{<:AbstractMatrix{<:Number}}'': jump_operators."""
function bloch_animation(; figure_theme=Theme(), kwargs...)
    fig, axis_bloch = setup_figure(; figure_theme, location=[1, 2]) 

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
end


#### SOURCE 
function setup_figure(; figure_theme, location=[1, 1])
    with_theme(figure_theme) do
        radius = 1/sqrt(2)
        fig = Figure(resolution=(1000, 770))
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
end

function blochquiverplot!(quiver_observables, vectorplots; fig=current_figure(), axis=current_axis(), isvisible=false)
    points, directions, inward_components, colorrange = quiver_observables
    
    arrow_object = arrows!(axis, points, directions, 
        color=inward_components, 
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
        label="quiver strength",
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

    initial_density = lift(slider_grid.Xcomponent_ϱ0, slider_grid.Ycomponent_ϱ0, slider_grid.Zcomponent_ϱ0) do X, Y, Z
        blochvector_to_vectorized_density(X, Y, Z, vectorized_npauli)
    end
    orbit_density = lift(control_effective, lindbladian_matrices, initial_density, slider_grid.time_horizon) do control_eff, lindbladian_mat, ϱ0, time_horizon
        simulate_gksl(control_eff, lindbladian_mat, ϱ0, time_horizon; saveat).u
    end
    orbit_bloch = lift(orbit_density) do orbit
        to_blochball(orbit, npauli)
    end

    orbit_color = lift(orbit_bloch) do orbit; 0.5f0 .+ norm.(orbit).^2; end
    orbit_begin = lift(first, orbit_bloch)
    orbit_end = lift(last, orbit_bloch)

    gate_endstate = lift(initial_density, menu.gate.selection) do ϱ0, gate_selec
        gates_pairing = first(pairings)
        to_blochball(
            [reshape(qgates[gates_pairing[gate_selec]] * ϱ0, 2, 2)],
            npauli
        )
    end

    control_inbloch = lift(control_spoint) do control
        if isapprox(control, zero(control))
            return [Point3f(first(control), last(control), 0.f0)]
        else
            c = maximum(abs.(control)) * 3
            return [Point3f(first(control)/c, last(control)/c, 0.f0)]
        end
    end

    ###### QUIVER OBSERVABLES
    quiver_points = ball_mesh()
    scalar_type = eltype(first(quiver_points))

    directions = Observable(zero(quiver_points))
    inward_components = Observable( zeros(scalar_type, length(quiver_points)) )
    colorrange = Observable( (zero(scalar_type), one(scalar_type)) )
    
    onany(toggles.quiver.active, lindbladian_matrices, control_effective, slider_grid.time_horizon) do isquiver_active, L, u, T
        if isquiver_active
            velocity_field = to_bloch_matrix(L, npauli)
            bvelocity = bloch_velocity(quiver_points, (velocity_field, u, T))
            minmax = extrema(bvelocity.inward_components)

            directions[] = bvelocity.directions
            inward_components[] = bvelocity.inward_components # FOR NOW ITS THE LENGTH
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
        inward_components,
        colorrange
    )

    return (orbit_observables, quiver_observables)
end

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

function setup_toggles!(; fig=current_figure(), location=[1, 2], sublocation=[1, 1])
    height = 15; width = 30; fontsize = 15;
    free_toggle = Toggle(fig, active=true, halign=:left, height=height, width=width)
    free_label = Label(fig, lift(x -> x ? "Free Hamilt." : "Trivial", free_toggle.active), fontsize=fontsize)

    control_toggle = Toggle(fig, active=true, halign=:left, height=height, width=width)
    control_label = Label(fig, lift(x -> x ? "Controlled" : "No control", control_toggle.active), fontsize=fontsize)

    noisy_toggle = Toggle(fig, active=true, halign=:left, height=height, width=width)
    noisy_label = Label(fig, lift(x -> x ? "Noisy" : "Isolated", noisy_toggle.active), fontsize=fontsize)

    quiver_toggle = Toggle(fig, active=false, halign=:left, height=height, width=width)
    quiver_label = Label(fig, lift(x -> x ? "Quiver" : "Orbit", quiver_toggle.active), fontsize=fontsize)

    content = [free_toggle free_label;
            control_toggle control_label;
            noisy_toggle noisy_label;
            quiver_toggle quiver_label]
    fig[location...][sublocation...] = grid!(
        content, 
        tellheight=false, 
        tellwidth=false
    )

    return (;
        gksl=(;
            free=free_toggle, 
            control=control_toggle, 
            noisy=noisy_toggle,
        ),
        quiver=quiver_toggle
    )
end

function setup_slidergrid!(; fig=current_figure(), location=[1, 2], sublocation=[1, 3])
    sg = SliderGrid(
        fig[location...][sublocation...],
        (label = L"\omega_{\mathrm{max}}", range = 0.f0:0.1f0:20.f0, format = "{:.1f}", startvalue = 1.f0),
        (label = L"T", range = 0.f0:0.1f0:10.f0, format = "{:.1f}", startvalue = 1.f0),
        (label = L"⟨ϱ(0),X⟩", range = -1.f0:0.01f0:1.f0, format = "{:.2f}", startvalue = 1.f0),
        (label = L"⟨ϱ(0),Y⟩", range = -1.f0:0.01f0:1.f0, format = "{:.2f}", startvalue = 0.f0),
        (label = L"⟨ϱ(0),Z⟩", range = -1.f0:0.01f0:1.f0, format = "{:.2f}", startvalue = 0.f0),
        width = 200,
        tellheight = false
    )

    return (; 
        control_amplitude=sg.sliders[1].value, 
        time_horizon=sg.sliders[2].value,
        Xcomponent_ϱ0=sg.sliders[3].value,
        Ycomponent_ϱ0=sg.sliders[4].value,
        Zcomponent_ϱ0=sg.sliders[5].value
    )
end

function setup_menu!(; fig=current_figure(), location=[1, 2], sublocation=[1, 2])
    gate_menu = Menu(
        fig, 
        options = [
            "Gate RX(π/2)",
            "Hadamard",
            "RZ(3π/2)"
        ],
        default = "Gate RX(π/2)",
    )
    fig[location...][sublocation...] = vgrid!(
        gate_menu, 
        tellheight = true, 
        width = 200
    )

    return (; gate=gate_menu)
end


function setup_selectpoint!(; fig=current_figure(), location=[1, 2], sublocation=[4, 1])
    axis_control = Axis(
        fig[location...][sublocation...],
        title=L"\mathrm{Control}~ \omega",
        limits=(-1.0, 1.0, -1.0, 1.0),
        aspect=AxisAspect(1),
        xlabel=L"ω_x",
        ylabel=L"ω_y"
    )
    hidexdecorations!(axis_control, label=false, grid=false)
    hideydecorations!(axis_control, label=false, grid=false)

    Makie.deactivate_interaction!(axis_control, :rectanglezoom)
    spoint = select_point(axis_control.scene, marker=:circle)
    
    arrow_head = lift(spoint) do p
        if isapprox(p, zero(p))
            [(8 .* p) ./ 10]
        else
            [(8 .* p ./ norm(p)) ./ 10]
        end
     end
    arrow_base = [Point2(0.f0, 0.f0)]
    arrows!(axis_control, arrow_base, arrow_head; color=:red, linewidth=5, arrowsize=20)
    
    return spoint
end