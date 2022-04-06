function ssfm_solve(sim::Simulation, coeffs::Coefficients, state::InitialState)
    pyplot()
    time_steps = Int(floor(sim.T/sim.dt))
    space_steps = Int(floor(sim.S/sim.ds))
    time = LinRange(0, sim.T, time_steps)
    space = LinRange(-sim.S/2, sim.S/2, space_steps)

    print("Simulation diagnostics:\n")
    @printf("Time  : from %6.3i ms, to  %6.3i ms; %6d steps \n", time[1]*1e3, time[time_steps]*1e3, time_steps)
    @printf("Space : from %6.3f mm, to  %6.3f mm; %6d steps \n", space[1]*1e3, space[space_steps]*1e3, space_steps)

    # Spatial frequency range computation
    k = 2*pi* LinRange(-1/(2*sim.ds), 1/(2*sim.ds), space_steps)
    k = fftshift(k)
    #fig = Plots.plot(k, show=true)
    ψ = zeros(ComplexF64, space_steps, time_steps)
    ψ_spect = zeros(ComplexF64, space_steps, time_steps)
    
    #Implement external initial state 
    wave1(s::Float64) = sqrt(1/ (sqrt(2*pi)* state.width)) * exp.(-(s).^2 / (4*state.width^2))
    wave2(s::Float64) = sqrt(1/(2 * state.width)) * 2 ./(exp.(-(s/(state.width))) .+ exp.(s/(state.width)))
  
    if state.sech_flag == 0 # Gaussian Pulse
      wave = wave1
    else # Sech Pulse
      wave = wave2
    end
    fig = plot(space, wave1.(space), show=true)
    plot!(space, wave2.(space), show=true)
    ## check waveform integral normalization
    integral, error = quadgk(s -> abs(wave(s))^2, -sim.S/2, +sim.S/2)
    display(integral)
    display(error)
    waveform = wave.(space)
    

    ## [SPACE INDEX, TIME INDEX]
    ψ[:, 1] = waveform

    fwd_disp = exp.(sim.dt/2 .* coeffs.α * k.^2)
    fwd_curvature = exp.(sim.dt/2 .* coeffs.β.(space))
    # display("Gamma at peak")
    # display(coeffs.γ.(ψ[Int(floor(space_steps/2)),1]))
    @showprogress "Propagating the field... " for n = 1:time_steps-1
      ψ_spect[:, n] = fft(ψ[:, n])
      ψ_spect[:, n] = ψ_spect[:, n] .* fwd_disp .* fwd_curvature
      ψ[:, n] = ifft(ψ_spect[:, n])
      ψ[:, n+1] = ψ[:, n] .* exp.(sim.dt/2 .* coeffs.γ.(ψ[:, n]))  ## this is an Euler step
    end
    ψ_spect[:, time_steps] = fft(ψ[:, time_steps])

    tmargin_l = Int(floor(time_steps/2))
    tmargin_r = Int(ceil(time_steps/2))
    t_points = 10000
    z_points = 1000
    t_skip = Int(ceil(20 * state.width/sim.dt /t_points))
    z_skip = Int(ceil(space_steps/z_points))

    # Natural orientation choice for plot 
    fig2 = heatmap(time*1e3, 
                   space*1e3, 
                   abs.(ψ[:, :]).^2, 
                   show=true, 
                   title = "wavefunction", 
                   ylabel="space [mm]", 
                   xlabel="time [ms]", 
                   colorrange=(0, 1),
                   reuse=false)

    fig3 = plot(space*1e3,
                abs.(ψ[:, time_steps]).^2,
                title = "evolved wavefunction",
                xlabel="space [mm]", 
                reuse=false, 
                label="t=t_max")
    plot!(space*1e3, 
          abs.(ψ[:, 1]).^2, 
          show=true, 
          reuse=false, 
          label="t=0")
  end