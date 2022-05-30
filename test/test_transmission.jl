### TEST FOR THE TRANSMISSION COEFFICIENT FOR DIFFERENT SOLITON VELOCITIES AND BARRIER STRENGTHS
## PLOTTING RESULTS IN NORMALIZED COORDINATES AS GARDNER et al
using Pkg
Pkg.instantiate()
using SolitonBEC
using Printf
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
    T, #T
    T * 5e-4, #dt
    L, #S
    L * 5e-4, #ds
  )
  return num
end

#--------- Potential ---------

function barrier_height(energy::Float64)
  r = Potential(
    "barrier", # type 
    1e-6, #width
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
configs = []
num_barr = 100
barrier_list = LinRange(0, 0.2 * 15036/8, num_barr)
velocity_list = LinRange(0, 0.2, num_barr)

for vel in velocity_list
  for barrier_energy in barrier_list
    potential = barrier_height(barrier_energy * hbar)
    state = initial_state_velocity(vel * velocity_unit)
    numerics = adaptive_numerics(vel * velocity_unit, L, x0, velocity_unit)
    push!(configs, (numerics, khaykovich_gpe, std_apparatus, potential, state))
  end
end

mem_limit = 15000000000 #byte

## Soliton - barrier collision --------------------------------------------
#run_dynamics(configs[20*20 - 5]...)


T = zeros(Float64, length(velocity_list), length(barrier_list))
max = zeros(Float64, length(velocity_list), length(barrier_list))

nth = Threads.nthreads() #print number of threads
print("\nNumber of threads: ", nth, "\n")


for iv in axes(velocity_list, 1)
  print("\nVel. ", iv)
  Threads.@threads for ib in axes(barrier_list, 1)
    print("\n\t Barr. ", ib, " in vel. ", iv)

    (numerics, sim, app, pot, state) = configs[(iv-1) * num_barr + ib]
    
    # potential space index
    coeffs = get_coefficients(sim, app, pot, state)
    time, space, ψ, ψ_spect, max_amplitude = ssfm_propagate(numerics, coeffs)
    #potential_idx = Int64(floor((pot.position+numerics.S/2) / numerics.ds))

    T[iv, ib] = sum(abs.(ψ[Int(floor(length(space)/2)):end]) .^ 2 * numerics.ds)
    max[iv, ib] = max_amplitude
  end
end

write("T_discontinuity_better.bin", T)
#write(".bin", max)