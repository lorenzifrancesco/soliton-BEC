
function ssfm_solve(num::Numerics, coeffs::Coefficients)
  time_steps = Int(floor(num.T / num.dt))
  space_steps = Int(floor(num.S / num.ds))
  time = LinRange(0, num.T, time_steps)
  space = LinRange(-num.S / 2, num.S / 2, space_steps)

  # Spatial frequency range computation
  k = 2 * pi * LinRange(-1 / (2 * num.ds), 1 / (2 * num.ds), space_steps)
  k = fftshift(k)
  ψ = zeros(ComplexF64, space_steps, time_steps)
  ψ_spect = zeros(ComplexF64, space_steps, time_steps)

  waveform = coeffs.initial.(space)

  ## [SPACE INDEX, TIME INDEX]
  ψ[:, 1] = waveform

  fwd_disp = exp.(num.dt / 2 .* coeffs.α * k .^ 2)
  fwd_curvature = exp.(num.dt / 2 .* coeffs.β.(space))

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

# impossible to keep in memory the full time evolution
function ssfm_propagate_3d(num3D::Numerics_3D, coeffs3d::Coefficients_3D)
  time_steps = Int(floor(num3D.T / num3D.dt))
  axial_steps = Int(floor(num3D.S / num3D.ds))
  transverse_steps = Int(floor(num3D.Transverse / num3D.dtr))

  time = LinRange(0, num3D.T, time_steps)
  axial = LinRange(-num3D.S / 2, num3D.S / 2, axial_steps)
  x_axis =  LinRange(-num3D.Transverse / 2, num3D.Transverse / 2, transverse_steps)
  y_axis =  LinRange(-num3D.Transverse / 2, num3D.Transverse / 2, transverse_steps)

  # Spatial frequency range computation
  ks = 2 * pi * LinRange(-1 / (2 * num3D.ds), 1 / (2 * num3D.ds), axial_steps)
  ks = fftshift(ks)
  kx = 2 * pi * LinRange(-1 / (2 * num3D.dtr), 1 / (2 * num3D.dtr), transverse_steps)
  kx = fftshift(kx)
  ky = 2 * pi * LinRange(-1 / (2 * num3D.dtr), 1 / (2 * num3D.dtr), transverse_steps)
  ky = fftshift(ky)
  ψ = zeros(ComplexF64, transverse_steps, transverse_steps, axial_steps)
  ψ_spect = zeros(ComplexF64, transverse_steps, transverse_steps, axial_steps)
  # this is wrong
  waveform = Array{ComplexF64, 3}(undef, (transverse_steps, transverse_steps, axial_steps))
  axial_waveforms = coeffs3d.initial_axial.(axial)
  print("AODFNBAOSNBA")
  display(axial)
  idx = 1
  for x in x_axis
    idy = 1
    for y in y_axis
      waveform[idx, idy, :] .= coeffs3d.initial_radial.((x^2 + y^2).^(1/2))' .* axial_waveforms
      idy +=1
    end
    idx+=1
  end

  # Natural orientation choice for plot 
  # fig2 = heatmap(x_axis * 1e3,
  # y_axis * 1e3,
  # abs.(waveform[:, :, 1]),
  # show=true,
  # title="|ψ|^2",
  # ylabel="space [mm]",
  # xlabel="time [ms]",
  # colorrange=(0, 1),
  # reuse=false,
  # size=(800, 600))
  # display(fig2)

  ## [axial INDEX, TIME INDEX]
  ψ = waveform
  gr()
  fig_axis = plot(axial, abs2.(waveform[3, 3, :]), title="initial axial distribution")
  display(fig_axis)

  fwd_disp_s = exp.(num3D.dt / 2 .* coeffs3d.α * ks .^ 2)
  fwd_disp_x = exp.(num3D.dt / 2 .* coeffs3d.α * kx .^ 2)
  fwd_disp_y = exp.(num3D.dt / 2 .* coeffs3d.α * ky .^ 2)
  transverse_disp = Array{ComplexF64, 2}(undef, (transverse_steps, transverse_steps))

  transverse_disp .= fwd_disp_x .* fwd_disp_y'
  disp = Array{ComplexF64 , 3}(undef, (transverse_steps, transverse_steps, axial_steps))
  idk = 1
  for k in ks
    disp[:, :, idk] .= exp(num3D.dt / 2 * coeffs3d.α * k ^ 2) * transverse_disp
    idk+=1
  end

  fwd_s = exp.(num3D.dt / 2 .* coeffs3d.β.(axial))
  fwd_x = exp.(num3D.dt / 2 .* coeffs3d.β.(x_axis))
  fwd_y = exp.(num3D.dt / 2 .* coeffs3d.β.(y_axis))
  transverse_beta = Array{ComplexF64, 2}(undef, (transverse_steps, transverse_steps))
  transverse_beta .= fwd_x .* fwd_y'
  fwd_beta = Array{ComplexF64, 3}(undef, (transverse_steps, transverse_steps, axial_steps))
  ids = 1
  for s in axial
    fwd_beta[:, :, ids] = exp(num3D.dt / 2 * coeffs3d.α * s ^ 2) * transverse_beta
    idk+=1
  end
  display(abs2.(ψ[3,3,:]))

  max_amplitude = maximum(abs.(ψ).^2)
  for n = 1:time_steps-1
    ψ_spect = fft(ψ)
    ψ_spect .= ψ_spect .* disp
    ψ = ifft(ψ_spect)
    ψ .= ψ .* exp.(num3D.dt / 2 .* coeffs3d.γ.(ψ))  ## this is an Euler step
    if max_amplitude < maximum(abs.(ψ).^2)
      max_amplitude = maximum(abs.(ψ).^2)
    end
    #display(ψ[3, 3, :])
  end
  ψ_spect = fft(ψ)

  print("\naxial distribution:\n ")
  display(abs2.(ψ[3,3,:]))

  return time, axial, ψ, ψ_spect, max_amplitude
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


function pseudospectral_solve(num::Numerics, coeffs::Coefficients)
  time_steps = Int(floor(num.T / num.dt))
  time = LinRange(0, num.T, time_steps)

  space_steps = Int(floor(num.S / num.ds))
  space = LinRange(-num.S / 2, num.S / 2, space_steps)
  
  waveform = coeffs.initial.(space) |> complex

  # Scale factor in fourier transform are inessential, do not include dμx / dμk
  u0 = fft(waveform)
  
  tspan = (time[1], time[end])
  p = Para(num, coeffs)

  print("\n\nDefine the ODEProblem\n")
  problem = ODEProblem(propagation_function, u0, tspan, p)
  print("\n\n=========================== Launch Solve ============================\n")
  @time time, ψ_k = solve(problem, Tsit5(), dense=false,maxiters=1e5,progress=true)

  #print(abs2.(ψ_k))
  ψ_x = ψ_k

  for i in 1:length(time)
    ψ_x[:, i] = ifft(ψ_k[:, i])
  end

  return time, space, ψ_x, ψ_k
end