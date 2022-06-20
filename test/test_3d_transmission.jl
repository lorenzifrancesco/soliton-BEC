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
  "NPSE",
  "barrier",
)

γ = 0 * 1.77e-11
# --------- Apparata ---------
std_apparatus = Apparatus(
  mass, #m (conversion AMU -> kg)
  as *2*pi, # as
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
    l_perp * 1e-1, 
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

## ==================== Transmission grid configuration
configs = []
num_barr = 15
barrier_list = LinRange(0, 15036/16, num_barr)
velocity_list = LinRange(0, 1, num_barr)

for vel in velocity_list
  for barrier_energy in barrier_list
    potential = barrier_height(barrier_energy * hbar)
    state = initial_state_vel(vel * velocity_unit)
    numerics = adaptive_numerics(vel * velocity_unit, L, x0, velocity_unit)
    push!(configs, (numerics, khaykovich_gpe, std_apparatus, potential, state))
  end
end

mem_limit = 15000000000 #byte
plt_width = 800
plt_height = 600

## Soliton - barrier collision --------------------------------------------
#run_dynamics(configs[20*20 - 5]...)


T = zeros(Float64, length(velocity_list), length(barrier_list))

nth = Threads.nthreads() #print number of threads
print("\n-->Number of threads: ", nth)
#Threads.@threads
 for iv in axes(velocity_list, 1)
  
    for (ib, barrier_energy) in enumerate(barrier_list)

    (numerics, sim, app, pot, state) = configs[(iv-1) * num_barr + ib]
    # potential space index
    coeffs = get_coefficients_3d(sim, app, pot, state)
    time, space, ψ, ψ_spect = ssfm_propagate_3d(numerics, coeffs)
    axis_center = Int(floor((numerics.Transverse / 2/numerics.dtr)))
    #display(ψ - ψ_old)
    #potential_idx = Int64(floor((pot.position+numerics.S/2) / numerics.ds))
    # fig1 = plot(title="|ψ|^2 after collision with barrier",
    # xlabel="space [mm]",
    # ylabel="|ψ|^2",
    # reuse=false,
    # size=(800, 400),
    # legend=:topleft)
    # heatmap!(abs2.(ψ[:, :, :]))
    # display(fig1)    
    
    T[iv, ib] = sum(abs.(ψ[:, :, Int(floor(length(space)/2)):end]) .^ 2 * numerics.ds)
    print("\nT[", iv, ", ", ib ,"] = ", T[iv, ib])
  
  end
end

write("T_3D_GPU.bin", T)
