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
    t_off = model.t_off
    C₀, C₁ = model.C₀, model.C₁
    if t_on-dt/2 ≤ t < t_on+dt/2
        return (C₁-C₀)/dt
    elseif t_off-dt/2 ≤ t < t_off+dt/2
        return (C₀-C₁)/dt
    else
        return 0.0
    end
end
function concentration_field(pos, model)
    dt = model.timestep
    t = model.t * dt
    t_on = model.t_on
    t_off = model.t_off
    C₀, C₁ = model.C₀, model.C₁
    concentration_field(t, t_on, t_off, C₀, C₁)
end
function concentration_field(t, t_on, t_off, C₀, C₁)
    if t < t_on || t > t_off
        return C₀
    else
        return C₁
    end
end
##

##
L = 1000.0
extent = (L,)
dt = 0.1
C₀, C₁ = 0.1, 1.0 # μM
t_on = 10.0
t_off = Inf
properties = Dict(
    :C₀ => C₀, :C₁ => C₁, :t_on => t_on, :t_off => t_off,
    :concentration_field => concentration_field,
    :concentration_time_derivative => concentration_time_derivative
)
model = UnremovableABM(Brumley{1}, extent, dt; properties)
add_agent!(model; chemotactic_precision=1, receptor_gain=1, motor_gain=1)
add_agent!(model; chemotactic_precision=6, receptor_gain=1, motor_gain=1)
##

##
S(a) = a.state
ν(a) = min(a.turn_rate * (1 + exp(-a.motor_gain*a.state))/2, 1/dt)
adata = [S, ν]
adf, mdf = run!(model, round(Int, 20/dt); adata)
##

##
ts = unique(adf.step) .* dt
C = @. concentration_field(ts, t_on, t_off, C₀, C₁)
_S = vectorize_adf_measurement(adf, :S)
_ν = vectorize_adf_measurement(adf, :ν)
plot(
    plot(ts.-t_on, C, ylab="C (μM)", xformatter=_->""),
    plot(formatter=_->"",axis=false),
    plot(ts.-t_on, _S, ylab="S", xlab="time (s)"),
    plot(ts.-t_on, _ν, ylab="ν (1/s)", xlab="time (s)"),
    leg=false, lw=4, layout=(2,2), size=(900,600),
    left_margin=-5Plots.mm, bottom_margin=-4Plots.mm
)
##