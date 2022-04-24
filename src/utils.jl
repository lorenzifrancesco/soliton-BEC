function graphics_config()
    Plots.plotly()
    Plots.default(size=(1300, 800),
        guidefont=("times", 10),
        tickfont=("times", 10))
end

function mem_estimate(num::Numerics)
    return 2 * num.T / num.dt * num.S / num.ds * 128
end

function describe_simulation(sim::Simulation, app::Apparatus, state::InitialState)

end