##
using MicrobeAgents, DataFrames, Statistics, Random
using Plots; Plots.gr()
include("../plot_defaults.jl")
include("../utils.jl")
##

##
function concentration_time_derivative(pos, model)
    dt = model.timestep
    t = model.t * dt
    t_on = model.t_on
    C₀, C₁ = model.C₀, model.C₁
    if t_on-dt/2 ≤ t < t_on+dt/2
        return (C₁-C₀)/dt
    else
        return 0.0
    end
end
function concentration_field(pos, model)
    dt = model.timestep
    t = model.t * dt
    t_on = model.t_on
    C₀, C₁ = model.C₀, model.C₁
    concentration_field(t, t_on, C₀, C₁)
end
concentration_field(t, t_on, C₀, C₁) = t < t_on ? C₀ : C₁
##

##
L = 1000.0
extent = (L,)
dt = 0.1
C₀, C₁ = 0.0, 3.0 # μM
t_on = 10.0
properties = Dict(
    :C₀ => C₀, :C₁ => C₁, :t_on => t_on,
    :concentration_field => concentration_field,
    :concentration_time_derivative => concentration_time_derivative
)
model = UnremovableABM(Celani{1}, extent, dt; properties)
add_agent!(model; gain=1.3)
##

##
S₀(a) = a.state[1]
S₁(a) = a.state[2]
S₂(a) = a.state[3]
S(a) = -a.state[4]+1
ν(a) = max(0, a.turn_rate * a.state[4])
adata = [S, ν, S₀, S₁, S₂]
adf, mdf = run!(model, round(Int, 2t_on/dt); adata)
##

##
ts = adf.step .* dt
C = @. concentration_field(ts, t_on, C₀, C₁)
plot(
    plot(ts.-t_on, C, ylab="C (μM)", xformatter=_->""),
    #plot(ts.-t_on, [adf.S₀ adf.S₁ adf.S₂], ylab="Markovian variables", xformatter=_->""),
    plot(ts.-t_on, adf.S, ylab="S", xlab="time (s)"),
    plot(ts.-t_on, adf.ν, ylab="ν (1/s)", xlab="time (s)"),
    leg=false, lw=4, layout=(2,2), size=(900,600),
    left_margin=-5Plots.mm, bottom_margin=-4Plots.mm
)
##