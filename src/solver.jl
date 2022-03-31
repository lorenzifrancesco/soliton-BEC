using FFTW
using Plots
import Plotly
using ProgressMeter

function solve(sim::Simulation, coeffs::Coefficients)
    time_steps = Int(floor(sim.T/sim.dt))
    space_steps = Int(floor(sim.S/sim.ds))
    
    # disp_ratio = sim.S / LD
    # nonl_ratio = sim.S / LNL
    print("Simulation diagnositics:\n")
    # @printf("\tT/T0  = %6.4f  window form factor\n", sim.T/state.S0)
    # @printf("\tL/LD  = %6.4f\n", disp_ratio)
    # @printf("\tL/LNL = %6.4f\n", nonl_ratio)
    # @printf("\tN     = %6.4f \n", sqrt(LD / LNL))
    @printf("time steps : %6d \n", time_steps)
    @printf("space steps: %6d \n", space_steps)
  
    time = LinRange(-sim.T/2, sim.T/2, time_steps)
    space = LinRange(0, sim.S, space_steps)
    
    # Spatial frequency range computation
    k = 2*pi* LinRange(-1/(2*sim.ds), 1/(2*sim.ds), space_steps)
    k = fftshift(k)
    ψ = zeros(ComplexF64, space_steps, time_steps)
    ψ_spect = zeros(ComplexF64, space_steps, time_steps)
    
    #Implement external initial state 
    waveform = 1 .* 2 ./(exp.(space/10e-12).+exp.(-space/10e-12))
  
    ## [SPACE INDEX, TIME INDEX]
    ψ[:, 1] = waveform
    #  fig = Plots.plot(time, abs.(ψ[:, 1]), show=true)
    ψ_spect[:, 1] = fft(ψ[:, 1])
  
    fwd_disp = exp.(sim.ds * coeffs.α)
    curvature = coeffs.β(time)
    print(typeof(coeffs.γ(1.0+im)))
    @showprogress "Propagating the field... " for n = 1:time_steps-1
      ψ_spect[:, n] = fft(ψ[:, n])
      ψ_spect[:, n+1] = ψ_spect[:, n] .* fwd_disp .* exp(curvature)
      ψ[:, n+1] = ifft(ψ_spect[:, n+1])
      ψ[:, n+1] = ψ[:, n+1] .* exp.(im * sim.dt * coeffs.γ.(ψ[:, n]))  ## this is an Euler step
    end
    ψ_spect[:, time_steps] = fft(ψ[:, time_steps])
  
    tmargin_l = Int(floor(time_steps/2))
    tmargin_r = Int(ceil(time_steps/2))
    t_points = 10000
    z_points = 1000
    t_skip = Int(ceil(20 * state.S0/sim.dt /t_points))
    z_skip = Int(ceil(space_steps/z_points))
    fig2 = Plots.heatmap(abs.(ψ[1:z_skip:space_steps, tmargin_l:t_skip:tmargin_r, ]), show=true)
    fig3 = Plots.surface(abs.(ψ[1:z_skip:space_steps, tmargin_l:t_skip:tmargin_r]), show=true)
  end
  