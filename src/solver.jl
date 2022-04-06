function ssfm_solve(sim::Simulation, coeffs::Coefficients, state::InitialState)
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
    if state.sech_flag == 0 # Gaussian Pulse
      waveform = 1/ (sqrt(2*pi)* state.width) * exp.(-(space).^2 / (2*state.width^2))
    else # Sech Pulse
      waveform = 1/(sqrt(2*pi)* state.width) * 2 ./(exp.(-(space/(2*state.width))) .+ exp.(space/(2*state.width)))
    end
    #fig = Plots.plot(space, waveform, show=true)

    ## [SPACE INDEX, TIME INDEX]
    ψ[:, 1] = waveform

    fwd_disp = exp.(sim.dt/2 * coeffs.α * k.^2)
    fwd_curvature = exp.(sim.dt/2 .* coeffs.β.(space))

    @showprogress "Propagating the field... " for n = 1:time_steps-1
      ψ_spect[:, n] = fft(ψ[:, n])
      ψ_spect[:, n] = ψ_spect[:, n] .* fwd_disp .* fwd_curvature
      ψ[:, n] = ifft(ψ_spect[:, n])
      ψ[:, n+1] = ψ[:, n] .* exp.(sim.dt/2 * coeffs.γ.(ψ[:, n]))  ## this is an Euler step
    end
    ψ_spect[:, time_steps] = fft(ψ[:, time_steps])

    tmargin_l = Int(floor(time_steps/2))
    tmargin_r = Int(ceil(time_steps/2))
    t_points = 10000
    z_points = 1000
    t_skip = Int(ceil(20 * state.width/sim.dt /t_points))
    z_skip = Int(ceil(space_steps/z_points))

    # Natural orientation choice for plot 
    fig2 = Plots.heatmap!(time[1:t_skip:time_steps]*1e3, space[1:z_skip:space_steps]*1e3, abs.(ψ[1:z_skip:space_steps, 1:t_skip:time_steps]), show=true, title = "wavefunction", ylabel="space [mm]", xlabel="time [ms]", size = (1920, 1080))
    #fig3 = Plots.surface(time*1e3, space*1e3, abs.(ψ[:, :]), show=true, title = "wavefunction", ylabel="space [mm]", xlabel="time [ms]")
  end