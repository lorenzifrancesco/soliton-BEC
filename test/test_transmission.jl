### TEST FOR THE TRANSMISSION COEFFICIENT FOR DIFFERENT SOLITON VELOCITIES AND BARRIER STRENGTHS
## PLOTTING RESULTS IN NORMALIZED COORDINATES AS GARDNER et al
using SolitonBEC
using Printf
using Plots
using Elliptic
using Distributed
#using FileIO
hbar = 6.62607015e-34 / (2 * pi)
width = 1.42e-6
mass = 6.941 * 1.660539e-27
omega_perp = 2 * pi * 710
l_perp = sqrt(hbar / (mass * omega_perp))
N = 4e3
as = -0.21e-9

interaction_g = abs(2*hbar^2 * as / mass / l_perp^2)
ggg = 2*hbar*omega_perp*abs(as)

print("\ninteraction_g: ", interaction_g)

# mass as total mass or mass per atom??
time_unit = hbar^2/(mass * ggg * N)
space_unit = hbar^3/(mass * ggg^2 * N^2)
velocity_unit = ggg * N / hbar
print("\n\ttime_unit: ", time_unit)
print("\n\tspace_unit: ", space_unit)
print("\n\tvelocity unit: ", velocity_unit, " m/s")

# --------- Simulation ---------
khaykovich_gpe = Simulation(
  "GPE",
  "barrier",
)
# --------- Simulation ---------
khaykovich_npse = Simulation(
  "NPSE",
  "barrier",
)

γ = 0 * 1.77e-11
# --------- Apparata ---------
std_apparatus = Apparatus(
  mass, #m (conversion AMU -> kg)
  as, # as
  omega_perp, # ω_perp
  N, #N
  γ, #γ
  2 * pi * 4,# ω_z
)

# Energy = 1 / 2 * (std_apparatus.as^2 / (width^2) + std_apparatus.as^2 / (2 * l_perp^2) + width^2 / std_apparatus.as^2 - interaction_g * std_apparatus.as^3 / (l_perp * width^2)) * hbar * omega_perp
L = 2 * space_unit

v0 = 1*velocity_unit
x0 =-L/4
# --------- InitialStates ---------   
function initial_state_velocity(velocity::Float64)
  init = InitialState1 = InitialState(
    "sech", #type
    width*1.8, # width
    velocity, #v0
    x0,
  )
  return init
end

# --------- Numerics ---------
function adaptive_numerics(velocity::Float64, L, x0, velocity_unit)
  if (velocity==0)
    T = abs(x0)/ velocity_unit * 2
  else
    T = abs(x0)/ velocity * 2
  end
  num = Numerics(
    T, #
    T * 1e-3, #dt
    L, #S
    L * 1e-3, #ds
  )
  return num
end

#--------- Potential ---------

function barrier_height(energy::Float64)
  r = Potential(
    "barrier", # type 
    L * 1e-3, #width
    0, #position
    energy, #energy
    0, #ϵ
    10e-6 #a
  )
  return r
end

energy_unit = mass * N^2 * ggg^2 / hbar^2
print("\n\tenergy unit: ", energy_unit, " J")

GSEnergy = energy_unit * (-N/24)
#print("\nground state energy: ", GSEnergy/hbar, " hbar\n")
phys_vel = 4.918e-3
normd_vel = phys_vel*hbar/ggg/N
Energy = GSEnergy + energy_unit * normd_vel^2*N/2
#print("\ntraveling state energy: ", Energy/hbar, " hbar\n")


## ==================== Transmission grid configuration
configs_GPE = []
configs_NPSE = []
num_barr = 50
barrier_list = LinRange(0, 0.3, num_barr)
velocity_list = LinRange(0, 0.3, num_barr)

for vel in velocity_list
  for barrier_energy in barrier_list
    potential = barrier_height(barrier_energy * energy_unit /500000)
    state = initial_state_velocity(vel * velocity_unit)
    numerics = adaptive_numerics(vel * velocity_unit, L, x0, velocity_unit)
    push!(configs_GPE, (numerics, khaykovich_gpe, std_apparatus, potential, state))
    push!(configs_NPSE, (numerics, khaykovich_npse, std_apparatus, potential, state))

  end
end

mem_limit = 15000000000 #byte
pyplot()
plt_width = 800
plt_height = 600

## Soliton - barrier collision --------------------------------------------
#run_dynamics(configs[20*20 - 5]...)


Tg = zeros(Float64, length(velocity_list), length(barrier_list))
Tn = zeros(Float64, length(velocity_list), length(barrier_list))
ΔT = zeros(Float64, length(velocity_list), length(barrier_list))
max_amplitude_g = zeros(Float64, length(velocity_list), length(barrier_list))
max_amplitude_n = zeros(Float64, length(velocity_list), length(barrier_list))
nth = Threads.nthreads() #print number of threads
print("\n--number of threads: ", nth)

 for iv in axes(velocity_list, 1)
  
  Threads.@threads for ib in axes(barrier_list, 1)

    (numerics_g, sim_g, app_g, pot_g, state_g) = configs_GPE[(iv-1) * num_barr + ib]
    (numerics_n, sim_n, app_n, pot_n, state_n) = configs_NPSE[(iv-1) * num_barr + ib]

    # potential space index
    coeffs_g = get_coefficients(sim_g, app_g, pot_g, state_g)
    coeffs_n = get_coefficients(sim_n, app_n, pot_n, state_n)

    time, space, ψg, ψ_spectg = ssfm_solve(numerics_g, coeffs_g)
    time, space, ψn, ψ_spectn = ssfm_solve(numerics_n, coeffs_n)

    Tg[iv, ib] = sum(abs.(ψg[Int(floor(length(space)/2)):end, end]) .^ 2 * numerics_g.ds)
    Tn[iv, ib] = sum(abs.(ψn[Int(floor(length(space)/2)):end, end]) .^ 2 * numerics_n.ds)
    ΔT[iv, ib] = Tg[iv, ib] - Tn[iv, ib]
    max_amplitude_g[iv, ib] = maximum(abs.(ψg).^2)
    max_amplitude_n[iv, ib] = maximum(abs.(ψn).^2)

    print("\nTg[", iv, ", ", ib ,"] = ", Tg[iv, ib], "<--> Tn[", iv, ", ", ib ,"] = ", Tn[iv, ib])
  end
end

write("zoom_well_T_50_n.bin", Tn)
write("zoom_well_T_50_g.bin", Tg)
write("zoom_max_amplitude_g.bin", max_amplitude_g)
write("zoom_max_amplitude_n.bin", max_amplitude_n)