
function ssfm_solve(num::Numerics, coeffs::Coefficients)
  time_steps = Int(floor(num.T / num.dt))
  space_steps = Int(floor(num.S / num.ds))
  time = LinRange(0, num.T, time_steps)
  space = LinRange(-num.S / 2, num.S / 2, space_steps)

  # print("Simulation diagnostics:\n")
  # @printf("Time  : from %6.3i ms, to  %6.3i ms; %6d steps \n", time[1] * 1e3, time[time_steps] * 1e3, time_steps)
  # @printf("Space : from %6.3f mm, to  %6.3f mm; %6d steps \n", space[1] * 1e3, space[space_steps] * 1e3, space_steps)

  # Spatial frequency range computation
  k = 2 * pi * LinRange(-1 / (2 * num.ds), 1 / (2 * num.ds), space_steps)
  k = fftshift(k)
  ψ = zeros(ComplexF64, space_steps, time_steps)
  ψ_spect = zeros(ComplexF64, space_steps, time_steps)

  waveform = coeffs.initial.(space)
  # fig0 = plot(space * 1e3,
  #   abs.(waveform) .^ 2,
  #   title="initial configuration",
  #   xlabel="space [mm]",
  #   show=true,
  #   reuse=false,
  #   label="initial state")

  # fig1 = plot(space * 1e3,
  #   abs.(coeffs.β.(space)),
  #   show=true,
  #   reuse=false,
  #   label="potential")

  ## [SPACE INDEX, TIME INDEX]
  ψ[:, 1] = waveform

  fwd_disp = exp.(num.dt / 2 .* coeffs.α * k .^ 2)
  fwd_curvature = exp.(num.dt / 2 .* coeffs.β.(space))
  # display("Gamma at peak")
  # display(coeffs.γ.(ψ[Int(floor(space_steps/2)),1]))
  #@showprogress "Propagating the field... " 
  for n = 1:time_steps-1
    # if 1+2*app.as*(app.N - 1) * abs.(ψ[Int(floor(space_steps/2)), n])^2 > 1
    #   display("Collapse detected")
    #   return
    # end
    ψ_spect[:, n] = fft(ψ[:, n])
    ψ_spect[:, n] = ψ_spect[:, n] .* fwd_disp
    ψ[:, n] = ifft(ψ_spect[:, n])
    ψ[:, n+1] = ψ[:, n] .* exp.(num.dt / 2 .* coeffs.γ.(ψ[:, n])) .* fwd_curvature  ## this is an Euler step
  end
  ψ_spect[:, time_steps] = fft(ψ[:, time_steps])

  return time, space, ψ, ψ_spect
end

function ssfm_propagate(num::Numerics, coeffs::Coefficients)
  time_steps = Int(floor(num.T / num.dt))
  space_steps = Int(floor(num.S / num.ds))
  time = LinRange(0, num.T, time_steps)
  space = LinRange(-num.S / 2, num.S / 2, space_steps)

  # Spatial frequency range computation
  k = 2 * pi * LinRange(-1 / (2 * num.ds), 1 / (2 * num.ds), space_steps)
  k = fftshift(k)
  ψ = zeros(ComplexF64, space_steps)
  ψ_spect = zeros(ComplexF64, space_steps)

  waveform = coeffs.initial.(space)

  ## [SPACE INDEX, TIME INDEX]
  ψ = waveform

  fwd_disp = exp.(num.dt / 2 .* coeffs.α * k .^ 2)
  fwd_curvature = exp.(num.dt / 2 .* coeffs.β.(space))
  max_amplitude = maximum(abs.(ψ).^2)
  for n = 1:time_steps-1

    ψ_spect = fft(ψ)
    ψ_spect = ψ_spect .* fwd_disp
    ψ = ifft(ψ_spect)
    ψ= ψ .* exp.(num.dt / 2 .* coeffs.γ.(ψ)) .* fwd_curvature  ## this is an Euler step
    if max_amplitude < maximum(abs.(ψ).^2)
      max_amplitude = maximum(abs.(ψ).^2)
    end
  end
  ψ_spect = fft(ψ)

  return time, space, ψ, ψ_spect, max_amplitude
end

function ground_state_solve(num::Numerics, coeffs::Coefficients)
  time_steps = Int(floor(num.T / num.dt))
  space_steps = Int(floor(num.S / num.ds))
  time = LinRange(0, num.T, time_steps)
  space = LinRange(-num.S / 2, num.S / 2, space_steps)

  # print("Simulation diagnostics:\n")
  # @printf("Time  : from %6.3i ms, to  %6.3i ms; %6d steps \n", time[1] * 1e3, time[time_steps] * 1e3, time_steps)
  # @printf("Space : from %6.3f mm, to  %6.3f mm; %6d steps \n", space[1] * 1e3, space[space_steps] * 1e3, space_steps)

  # Spatial frequency range computation
  k = 2 * pi * LinRange(-1 / (2 * num.ds), 1 / (2 * num.ds), space_steps)
  k = fftshift(k)
  #fig = Plots.plot(k, show=true)
  ψ = zeros(ComplexF64, space_steps)
  ψ_spect = zeros(ComplexF64, space_steps)

  #Implement gaussian initial state 
  wave(s::Float64) = sqrt(1 / (sqrt(2 * pi) * num.S / 3)) * exp.(-(s) .^ 2 / (4 * num.S / 3^2))
  integral, error = quadgk(s -> abs(wave(s))^2, -num.S / 2, +num.S / 2)

  waveform = wave.(space) / sqrt(integral)

  ## [SPACE INDEX, TIME INDEX]
  ψ[:] = waveform
  fwd_disp = exp.(-im * num.dt / 2 .* coeffs.α * k .^ 2)
  fwd_curvature = exp.(-im * num.dt / 2 .* coeffs.β.(space))
  # display("Gamma at peak")
  # display(coeffs.γ.(ψ[Int(floor(space_steps/2)),1]))
  target_Δ = 1e-10
  Δ = 1
  let prev, current
    prev = 0
    while Δ > target_Δ
      # if 1+2*app.as*(app.N - 1) * abs.(ψ[Int(floor(space_steps/2)), n])^2 > 1
      #   display("Collapse detected")
      #   return
      # end
      ψ_spect = fft(ψ)
      ψ_spect = ψ_spect .* fwd_disp
      ψ = ifft(ψ_spect)
      ψ = ψ .* exp.(im * num.dt / 2 .* coeffs.γ.(ψ)) .* fwd_curvature  ## this is an Euler step
      # renormalization
      #isplay(sum(abs.(ψ) .^ 2 * num.ds))
      current = sum(abs.(ψ) .^ 2 * num.ds)
      Δ = abs(prev - current)
      ψ = ψ / sqrt(current)
      prev = current
    end
  end
  ψ_spect = fft(ψ)
  return time, space, ψ, ψ_spect
end


function plot_dynamics(time, space, ψ, ψ_spect)
  pyplot()
  # Natural orientation choice for plot 
  fig2 = heatmap(time * 1e3,
    space * 1e3,
    abs.(ψ[:, :]) .^ 2,
    show=true,
    title="|ψ|^2",
    ylabel="space [mm]",
    xlabel="time [ms]",
    colorrange=(0, 1),
    reuse=false,
    size=(800, 600))

  # fig2_bis = heatmap(time * 1e3,
  #   space * 1e3,
  #   angle.(ψ[:, :]),
  #   show=true,
  #   title="|ψ|^2",
  #   ylabel="space [mm]",
  #   xlabel="time [ms]",
  #   colorrange=(0, 1),
  #   reuse=false,
  #   size=(800, 600))

  fig3 = plot(space * 1e3,
    abs.(ψ[:, end]) .^ 2,
    title="evolved |ψ|^2",
    xlabel="space [mm]",
    reuse=false,
    label="t=t_max")
  plot!(space * 1e3,
    abs.(ψ[:, 1]) .^ 2,
    show=true,
    reuse=false,
    label="t=0",
    size=(800, 600))
end


function plot_ground_state(time, space, ψ, ψ_spect)
  pyplot()
  # Natural orientation choice for plot 
  fig3 = plot(space * 1e3,
    abs.(ψ) .^ 2,
    title="stationary |ψ|^2",
    xlabel="space [mm]",
    reuse=false,
    show=true,
    label="t=t_max")
end