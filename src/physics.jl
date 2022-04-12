struct Simulation
    T::Float64
    dt::Float64
    S::Float64
    ds::Float64
    equation::String #("GPE", "NPSE")
    potential::String #("barrier", "ellipse")
end


struct Apparatus
    m::Float64 # kg ω
    as::Float64 # m
    ω_perp::Float64 # rad/s
    N::Int64
    γ::Float64
    ω_z::Float64
    ϵ::Float64
    a::Float64
end


struct InitialState
    width::Float64 # m
    sech_flag::Float64 # 1 or 0
end


struct Coefficients
    α::ComplexF64
    β::Function
    γ::Function
end

function β(sim::Simulation, app::Apparatus, state::InitialState, s::Float64)
    hbar = 6.62607015e-34 / (2 * pi)
    l_perp = sqrt(hbar / (app.m * app.ω_perp))
    l_z = sqrt(hbar / (app.m * app.ω_z))
    if (sim.equation == "NPSE")
        value = im * hbar / (8 * app.m)
    elseif (sim.equation == "GPE")
        value = -im * app.ω_perp
    end

    if sim.potential == "barrier"
        value += -im * hbar * exp(-(s - 15e-6)^2 / (5e-6)^2) * 9e12 / (8 * app.m)
    elseif sim.potential == "ellipse"
        target_error = 0.01 * app.a
        lower_ϕ = -pi
        upper_ϕ = 3 * pi
        let error
            let midpoint
                error = target_error+1
                while error > target_error
                    midpoint = (lower_ϕ + upper_ϕ) / 2
                    display(midpoint)
                    if s - app.a * E(midpoint, app.ϵ) > 0
                        lower_ϕ = midpoint
                    else
                        upper_ϕ = midpoint
                    end
                    error = abs(s - app.a * E(midpoint, app.ϵ))
                end
                display(midpoint)
                value +=  im *  hbar * (1 / app.a * (sqrt(1 - app.ϵ^2)) / (sin(midpoint)^2 + sqrt(1 - app.ϵ^2) * cos(midpoint)^2)^(3 / 2))^2 / (8*app.m)  * 150
            end
        end
    end
    return value
end

function γ(sim::Simulation, app::Apparatus, state::InitialState, ψ::ComplexF64)
    hbar = 6.62607015e-34 / (2 * pi)
    l_perp = sqrt(hbar / (app.m * app.ω_perp))
    l_z = sqrt(hbar / (app.m * app.ω_z))
    npse_gamma = app.γ * hbar * app.ω_perp / l_z^6 * 1e-4
    σ(ψ::ComplexF64) = l_perp * sqrt(sqrt(1 + 2 * app.as * (app.N - 1) * abs.(ψ)^2))
    if sim.equation == "NPSE"
        value = -im * hbar / (2 * app.m * σ(ψ)^2) - im * app.m * app.ω_perp^2 / (2 * hbar) * σ(ψ)^2 - im * 2 * hbar * app.as * (app.N - 1) / (app.m * σ(ψ)^2) * abs.(ψ)^2 - npse_gamma * app.N^2 * abs.(ψ)^4
    elseif sim.equation == "GPE"
        value = -im * 2 * hbar * app.as * (app.N - 1) / (app.m * l_perp^2) * abs.(ψ)^2
    end
    return value
end


function get_coefficients(sim::Simulation, app::Apparatus, state::InitialState)
    # implement with PhysicalConstants.CODATA2018.h
    hbar = 6.62607015e-34 / (2 * pi)
    l_perp = sqrt(hbar / (app.m * app.ω_perp))
    l_z = sqrt(hbar / (app.m * app.ω_z))
    npse_gamma = app.γ * hbar * app.ω_perp / l_z^6 * 1e-4
    @assert(sim.equation in ["NPSE", "GPE"])
    @assert(sim.potential in ["barrier", "ellipse", ""])

    α = -im * hbar / (2 * app.m)
    beta(s::Float64) = β(sim::Simulation, app::Apparatus, state::InitialState, s)
    gamma(ψ::ComplexF64) = γ(sim::Simulation, app::Apparatus, state::InitialState, ψ)

    return Coefficients(α, beta, gamma)
end


function run_simulation(sim::Simulation, app::Apparatus, state::InitialState)
    coeffs = get_coefficients(sim, app, state)
    describe_simulation(sim, app, state)
    ssfm_solve(sim, coeffs, state, app)
end