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
C₀, C₁ = 1.0, 3.0 # μM
t_on = 10.0
properties = Dict(
    :C₀ => C₀, :C₁ => C₁, :t_on => t_on,
    :concentration_field => concentration_field,
    :concentration_time_derivative => concentration_time_derivative
)
model = UnremovableABM(BrownBerg{1}, extent, dt; properties)
add_agent!(model; adaptation_time=1)
##

##
S(a) = a.state
ν(a) = min(a.turn_rate * exp(-a.motor_gain * a.state), 1/dt)
adata = [S, ν]
adf, mdf = run!(model, round(Int, 2t_on/dt); adata)
##

##
ts = adf.step .* dt
C = @. concentration_field(ts, t_on, C₀, C₁)
KD = model[1].receptor_binding_constant
Pb = @. 1 / (1 + KD/C)
plot(
    plot(ts.-t_on, C, ylab="C (μM)", xformatter=_->""),
    plot(ts.-t_on, Pb.*100, ylab="Pb (%)", xformatter=_->""),
    plot(ts.-t_on, adf.S.*100, ylab="S (%/s)", xlab="time (s)"),
    plot(ts.-t_on, adf.ν, ylab="ν (1/s)", xlab="time (s)"),
    leg=false, lw=4, layout=(2,2), size=(900,600),
    left_margin=-5Plots.mm, bottom_margin=-4Plots.mm
)
##