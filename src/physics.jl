struct Simulation
    T::Float64
    dt::Float64
    S::Float64
    ds::Float64
end

struct Apparatus
    m::Float64 # dB/km
    as::Float64 # ps^2/km
    omega_perp::Float64 # km^-1
    N::Int64
    l_perp::Float64
  end
  
  struct InitialState
    P0::Float64
    T0::Float64
    sech_flag::Float64
  end

function get_coefficients(a::Apparatus, s:InitialState)


    return α, β, γ
end

function run(a::Apparatus, s:InitialState)
    coeffs = get_coefficients(a, s)
    solve(coeffs)

end