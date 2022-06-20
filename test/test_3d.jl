### TEST FOR THE TRANSMISSION COEFFICIENT FOR DIFFERENT SOLITON VELOCITIES AND BARRIER STRENGTHS
using SolitonBEC
using Printf
using Elliptic
using Distributed
using Plots
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

γ = 0 * 1.77e-11 # three body losses, phenomenological
# --------- Apparata ---------
std_apparatus = Apparatus(
  mass, #m (conversion AMU -> kg)
  as * 2*pi, # as
  omega_perp, # ω_perp
  N, #N
  γ, #γ
  2 * pi * 4,# ω_z
)

L = 2 * space_unit

v0 = 1*velocity_unit
x0 =-L/4
# --------- InitialStates ---------   
function initial_state_vel(velocity::Float64)
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
  num = Numerics_3D(
    T, #T
    T * 1e-3, #dt
    L, #S
    L * 1e-3, #ds
    l_perp, 
    l_perp * 5e-2, 
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

## Configurations
vel = 0.007
energy = 3000 *hbar
configs = [
  (adaptive_numerics(vel/2, L, x0, velocity_unit), khaykovich_gpe, barrier_height(energy), std_apparatus, initial_state_vel(vel)),
  #(num, khaykovich_gpe, barrier, std_apparatus, InitialState1), 
]
mem_limit = 15000000000 #byte
cnt = 1

angle = LinRange(0, 2 * pi, 100)

for (num, sim, pot, app, state) in configs

  global cnt

  coeffs = get_coefficients_3d(sim, app, pot, state)
  print("\nLaunch pseudospectral solver")
  time, axial, ψ_abs2_result = @time ssfm_solve_3d(num, coeffs)
  gr()
  print(size(ψ_abs2_result))
  fig = heatmap(abs.(ψ_abs2_result))
  display(fig)
  print("Maximum time: ", time[end])
  print("\nMaximum space: ", axial[end])
  cnt += 1

end