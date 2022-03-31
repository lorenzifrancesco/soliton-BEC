using FFTW
using Plots
import Plotly
using ProgressMeter

function solve(sim::Simulation, coeffs)
    time_steps = Int(floor(sim.T/sim.dt))
    space_steps = Int(floor(sim.S/sim.ds))
    
    disp_ratio = apparatus.L / LD
    nonl_ratio = apparatus.L / LNL
    print("Simulation diagnositics:\n")
    # @printf("\tT/T0  = %6.4f  window form factor\n", apparatus.T/state.T0)
    # @printf("\tL/LD  = %6.4f\n", disp_ratio)
    # @printf("\tL/LNL = %6.4f\n", nonl_ratio)
    # @printf("\tN     = %6.4f \n", sqrt(LD / LNL))
    @printf("time steps : %6d \n", time_steps)
    @printf("space steps: %6d \n", space_steps)
  
    time = LinRange(-apparatus.T/2, apparatus.T/2, time_steps)
    space = LinRange(0, apparatus.L, space_steps)
    
    # Frequency range computation
    w = 2*pi* LinRange(-1/(2*apparatus.dt), 1/(2*apparatus.dt), time_steps)
    w = fftshift(w)
    ψ = zeros(ComplexF64, time_steps, space_steps)
    A_spect = zeros(ComplexF64, time_steps, space_steps)
    
    if (state.sech_flag == 1)
      waveform = sqrt(state.P0) .* 2 ./(exp.(time/state.T0).+exp.(-time/state.T0))
      if collision == 1
        waveform = sqrt(state.P0) .* ( 2 ./ (exp.((time.-10*state.T0)./state.T0).+exp.(-(time.-10*state.T0)./state.T0)) + 2 ./ (exp.((time.+10*state.T0)./state.T0).+exp.(-(time.+10*state.T0)./state.T0)))
      end
    else
      waveform = sqrt(state.P0) .* exp.(-(time/state.T0) .^ 2 ./ 2)
    end
    
    A[:, 1] = waveform
    #  fig = Plots.plot(time, abs.(A[:, 1]), show=true)
    A_spect[:, 1] = fft(A[:, 1])
  

    # Computation for alpha independent of time and beta dependent
    fwd_disp = exp.(apparatus.dz * (im * apparatus.beta_2/2 * w .^ 2 .- alpha/2))
    @showprogress "Propagating the field... " for n = 1:space_steps-1
      ψ_spect[:, n] = fft(ψ[:, n])
      ψ_spect[:, n+1] = ψ_spect[:, n] .* fwd_disp
      ψ[:, n+1] = ifft(ψ_spect[:, n+1])
      ψ[:, n+1] = ψ[:, n+1] .* exp.(im * sim.ds * apparatus.gamma .* (abs.(ψ[:, n+1])).^2)  ## this is an euler step
    end
    A_spect[:, space_steps] = fft(A[:, space_steps])
  
    tmargin_l = Int(floor(time_steps/2 - 10 * state.T0/apparatus.dt))
    tmargin_r = Int(ceil(time_steps/2 + 10 * state.T0/apparatus.dt))
    t_points = 10000
    z_points = 1000
    t_skip = Int(ceil(20 * state.T0/apparatus.dt /t_points))
    z_skip = Int(ceil(space_steps/z_points))
    fig2 = Plots.heatmap(abs.(ψ[tmargin_l:t_skip:tmargin_r, 1:z_skip:space_steps]), show=true)
    fig3 = Plots.surface(abs.(ψ[tmargin_l:t_skip:tmargin_r, 1:z_skip:space_steps]), show=true)
  end

  function mem_estimate(apparatus::apparatus)
    return 2* apparatus.T/apparatus.dt * apparatus.L/apparatus.dz * 128
  end

  function 
    
  end

  