function setup_toggles!(; fig=current_figure(), location=[1, 2], sublocation=[1, 1])
    height = 15; width = 30; fontsize = 15;
    free_toggle = Toggle(fig, active=true, halign=:left, height=height, width=width)
    free_label = Label(fig, lift(x -> x ? L"H_0 \neq 0" : L"H_0 = 0", free_toggle.active), fontsize=fontsize)

    control_toggle = Toggle(fig, active=true, halign=:left, height=height, width=width)
    control_label = Label(fig, lift(x -> x ? "Controlled" : "No control", control_toggle.active), fontsize=fontsize)

    noisy_toggle = Toggle(fig, active=true, halign=:left, height=height, width=width)
    noisy_label = Label(fig, lift(x -> x ? "Noisy" : "Isolated", noisy_toggle.active), fontsize=fontsize)

    quiver_toggle = Toggle(fig, active=false, halign=:left, height=height, width=width)
    quiver_label = Label(fig, lift(x -> x ? "Velocity" : "Orbit", quiver_toggle.active), fontsize=fontsize)

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
        (label = L"Ω_{\mathrm{max}}", range = 0.f0:0.1f0:20.f0, format = "{:.1f}", startvalue = 1.f0),
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
        title=L"\mathrm{Angular~ Control}~ Ω",
        limits=(-1.0, 1.0, -1.0, 1.0),
        aspect=AxisAspect(1),
        xlabel=L"Ω_x",
        ylabel=L"Ω_y"
    )
    hidexdecorations!(axis_control, label=false, grid=false)
    hideydecorations!(axis_control, label=false, grid=false)

    Makie.deactivate_interaction!(axis_control, :rectanglezoom)
    spoint = select_point(axis_control.scene, marker=:circle, color=:red)
    
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