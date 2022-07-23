using CUDA.CUFFT
hbar = 6.62607015e-34 / (2 * pi)

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
  gr()
  fig_axis = plot(space, abs2.(ψ), title="axial distribution", reuse=false)
  display(fig_axis)
  return time, space, ψ, ψ_spect, max_amplitude
end

function ssfm_solve_3d(num3D::Numerics_3D, coeffs3d::Coefficients_3D)
  time_steps = Int(floor(num3D.T / num3D.dt))  
  axial_steps = Int(floor(num3D.S / num3D.ds))
  transverse_steps = Int(floor(num3D.Transverse / num3D.dtr))

  time = CuArray(LinRange(0, num3D.T, time_steps))
  axial = LinRange(-num3D.S / 2, num3D.S / 2, axial_steps)
  x_axis =  LinRange(-num3D.Transverse / 2, num3D.Transverse / 2, transverse_steps)
  y_axis = LinRange(-num3D.Transverse / 2, num3D.Transverse / 2, transverse_steps)

  # Spatial frequency range computation
  ks = 2 * pi * Array(LinRange(-1 / (2 * num3D.ds), 1 / (2 * num3D.ds), axial_steps))
  ks = fftshift(ks)
  kx = 2 * pi * Array(LinRange(-1 / (2 * num3D.dtr), 1 / (2 * num3D.dtr), transverse_steps))
  ky = 2 * pi * Array(LinRange(-1 / (2 * num3D.dtr), 1 / (2 * num3D.dtr), transverse_steps))

  ψ = CuArray{ComplexF64, 3}(undef, (transverse_steps, transverse_steps, axial_steps))
  ψ_spect = CuArray{ComplexF64, 3}(undef, (transverse_steps, transverse_steps, axial_steps))
  # this is wrong
  waveform = Array{ComplexF64, 3}(undef, (transverse_steps, transverse_steps, axial_steps))

  axial_waveforms = coeffs3d.initial_axial.(axial)
  idx = 1
  for x in x_axis
    idy = 1
    for y in y_axis
      waveform[idx, idy, :] .= coeffs3d.initial_radial.((x^2 + y^2).^(1/2))' .* axial_waveforms
      idy +=1
    end
    idx+=1
  end

  ## [axial INDEX, TIME INDEX]
  ψ = CuArray(waveform)

  fwd_disp_s = exp.(num3D.dt / 2 .* coeffs3d.α * ks .^ 2)
  fwd_disp_x = exp.(num3D.dt / 2 .* coeffs3d.α * -kx .^ 2)
  fwd_disp_y = exp.(num3D.dt / 2 .* coeffs3d.α * ky .^ 2)
  print("\n\nsize of x: ", size(fwd_disp_x))
  transverse_disp = Array{ComplexF64, 2}(undef, (transverse_steps, transverse_steps))

  transverse_disp .= fftshift(fwd_disp_x .* (fwd_disp_y)')
  print("\n\nTransverse_disp", typeof(transverse_disp),"\n\n")

  disp = Array{ComplexF64 , 3}(undef, (transverse_steps, transverse_steps, axial_steps))
  idk = 1
  for k in ks
    disp[:, :, idk] .= exp(num3D.dt / 2 * coeffs3d.α * k ^ 2) * transverse_disp
    idk+=1
  end
  disp = CuArray(disp)

  fwd_beta = Array{ComplexF64, 3}(undef, (transverse_steps, transverse_steps, axial_steps))
  idx = 1
  for x in x_axis
    idy = 1
    for y in y_axis
      fwd_beta[idx, idy, :] .= exp.(-im / hbar * num3D.dt / 2 * coeffs3d.confinment.((x^2 + y^2).^(1/2))) * exp.(num3D.dt / 2 * coeffs3d.β.(axial))
      idy +=1
    end
    idx+=1
  end

  fwd_beta = CuArray(fwd_beta)
  ###############
  fig = heatmap((Array(angle.(fwd_beta[:, :, 1]))))
  #display(fig)
  ###############
  ψ_abs2_result = Array{ComplexF64, 2}(undef, (axial_steps, time_steps))
  cross_section = Array{ComplexF64, 2}(undef, (axial_steps, time_steps))
  square_distance_mask = Array{ComplexF64, 2}(undef, (transverse_steps, transverse_steps))
  idx = 1
  for x in x_axis
    idy = 1
    for y in y_axis
      square_distance_mask[idx, idy] = (x^2 + y^2)
      idy +=1
    end
    idx+=1
  end

  square_distance_mask = CuArray(square_distance_mask)


  ψ_abs2_result[:, 1] = sum(abs2.(ψ), dims=(1, 2))
  print("\n\nsize(SUM waveform): ", size(sum(ψ, dims=(1, 2))))
  print("\n\nsize(SUM waveform .* square_distance_mask): ", size(sum(ψ.*square_distance_mask, dims=(1, 2))[1, 1, :]))
  print("\n\nsize( ψ_abs2_result[:, 1]) ", size(ψ_abs2_result[:, 1]))
  ψ_abs2_result = CuArray(ψ_abs2_result)
  cross_section = CuArray(cross_section)
  
  cross_section[:, 1] .= sum(abs2.(ψ) .* square_distance_mask , dims=(1, 2))[1, 1, :] ./ ψ_abs2_result[:, 1]

  max_amplitude = maximum(abs.(ψ).^2)
  Δ_radial=x_axis[2]-x_axis[1]
  Δ_axial=axial[2]-axial[1]
  #ψ .= ψ / (sqrt(sum(abs2.(ψ)* Δ_radial^2 * Δ_axial) ))
  print("\n\nintegral of modulus squared: " , sqrt(sum(abs2.(ψ) * Δ_axial * Δ_radial^2)))
  integral = sqrt(sum(abs2.(ψ)))

  #print("\n\nSEEK THE DIFFERENCES fwd_beta :" , fwd_beta[:, :, 1].-fwd_beta[:, :, 1]')

  for n = 1:time_steps-1
    ψ_spect = fft(ψ)
    ψ_spect .= ψ_spect .* disp
    
    ψ = ifft(ψ_spect)
    ψ .= ψ .* fwd_beta  .* exp.(num3D.dt / 2 .* coeffs3d.γ(ψ))

    # Renormalize
    print("\nrenorm factor: ", sqrt(sum(abs2.(ψ)))/integral)
    ψ .= ψ / sqrt(sum(abs2.(ψ)))* integral 
    ψ_abs2_result[:, n+1] = sum(abs2.(ψ), dims=(1, 2))
    cross_section[:, n+1] .= sum(abs2.(ψ) .* square_distance_mask , dims=(1, 2))[1, 1, :] ./ ψ_abs2_result[:, n+1]
  end

  ###############
  fig = heatmap(abs.(1e12*Array(abs2.(ψ[:, :, Int(floor(axial_steps/4))]))))
  display(fig)
  ###############


  ψ_spect = fft(ψ)
  print("\n\nComputation completed!\n\n")
  return time, axial, ψ_abs2_result, cross_section
end


function ssfm_propagate_3d(num3D::Numerics_3D, coeffs3d::Coefficients_3D)
  time_steps = Int(floor(num3D.T / num3D.dt))  
  axial_steps = Int(floor(num3D.S / num3D.ds))
  transverse_steps = Int(floor(num3D.Transverse / num3D.dtr))

  time = CuArray(LinRange(0, num3D.T, time_steps))
  axial = LinRange(-num3D.S / 2, num3D.S / 2, axial_steps)
  x_axis =  LinRange(-num3D.Transverse / 2, num3D.Transverse / 2, transverse_steps)
  y_axis = LinRange(-num3D.Transverse / 2, num3D.Transverse / 2, transverse_steps)

  # Spatial frequency range computation
  ks = 2 * pi * Array(LinRange(-1 / (2 * num3D.ds), 1 / (2 * num3D.ds), axial_steps))
  ks = fftshift(ks)
  kx = 2 * pi * Array(LinRange(-1 / (2 * num3D.dtr), 1 / (2 * num3D.dtr), transverse_steps))
  ky = 2 * pi * Array(LinRange(-1 / (2 * num3D.dtr), 1 / (2 * num3D.dtr), transverse_steps))

  ψ = CuArray{ComplexF64, 3}(undef, (transverse_steps, transverse_steps, axial_steps))
  ψ_spect = CuArray{ComplexF64, 3}(undef, (transverse_steps, transverse_steps, axial_steps))
  avg_cross_section = Array{ComplexF64, 2}(undef, (axial_steps, time_steps))

  # this is wrong
  waveform = Array{ComplexF64, 3}(undef, (transverse_steps, transverse_steps, axial_steps))

  axial_waveforms = coeffs3d.initial_axial.(axial)
  idx = 1
  for x in x_axis
    idy = 1
    for y in y_axis
      if idx<=2 || idx>=transverse_steps-2 || idy<=2 || idy>=transverse_steps-2
        waveform[idx, idy, :] .= 0 * axial_waveforms
      else waveform[idx, idy, :] .= coeffs3d.initial_radial.((x^2 + y^2).^(1/2))' .* axial_waveforms
      end
      idy +=1
    end
    idx+=1
  end

  ## [axial INDEX, TIME INDEX]
  ψ = CuArray(waveform)

  fwd_disp_s = exp.(num3D.dt / 2 .* coeffs3d.α * ks .^ 2)
  fwd_disp_x = exp.(num3D.dt / 2 .* coeffs3d.α * -kx .^ 2)
  fwd_disp_y = exp.(num3D.dt / 2 .* coeffs3d.α * ky .^ 2)
  transverse_disp = Array{ComplexF64, 2}(undef, (transverse_steps, transverse_steps))

  transverse_disp .= fftshift(fwd_disp_x .* (fwd_disp_y)')
  disp = Array{ComplexF64 , 3}(undef, (transverse_steps, transverse_steps, axial_steps))
  idk = 1
  for k in ks
    disp[:, :, idk] .= exp(num3D.dt / 2 * coeffs3d.α * k ^ 2) * transverse_disp
    idk+=1
  end
  disp = CuArray(disp)

  fwd_beta = Array{ComplexF64, 3}(undef, (transverse_steps, transverse_steps, axial_steps))
  idx = 1
  for x in x_axis
    idy = 1
    for y in y_axis
      if idx<=2 || idx>=transverse_steps-2 || idy<=2 || idy>=transverse_steps-2
        fwd_beta[idx, idy, :] .= 0 * exp.(-im / hbar * num3D.dt / 2 * coeffs3d.confinment.((x^2 + y^2).^(1/2))) * exp.(num3D.dt / 2 * coeffs3d.β.(axial))
      else
        fwd_beta[idx, idy, :] .= exp.(-im / hbar * num3D.dt / 2 * coeffs3d.confinment.((x^2 + y^2).^(1/2))) * exp.(num3D.dt / 2 * coeffs3d.β.(axial))
      end
      idy +=1
    end
    idx+=1
  end

  fwd_beta = CuArray(fwd_beta)
  integral = sqrt(sum(abs2.(ψ)))
  max_amplitude = maximum(abs.(ψ).^2)
  for n = 1:time_steps-1
    ψ_spect = fft(ψ)
    ψ_spect .= ψ_spect .* disp
    ψ = ifft(ψ_spect)
    ψ .= ψ .* exp.(num3D.dt / 2 .* coeffs3d.γ(ψ))  .* fwd_beta  ## this is an Euler step
    # Renormalize
    ψ .= ψ / sqrt(sum(abs2.(ψ)))* integral 
    # if max_amplitude < maximum(abs.(ψ).^2)
    #   max_amplitude = maximum(abs.(ψ).^2)
    # end
    #display(ψ[3, 3, :])
    
    # Average cross section computed as mean variance
  end
  ψ_spect = fft(ψ)

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