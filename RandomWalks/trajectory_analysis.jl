##
using MicrobeAgents, DataFrames, Distributions, Random
using Plots
include("../plot_defaults.jl")
# initialize a random number generator for reproducibility
rng = MersenneTwister(45)
##

## simulation parameters

##
d = 3
L = 300
dt = 0.01
extent = ntuple(_ -> L, d)
T = 60
nsteps = round(Int, T/dt)
nbacteria = 100
##

## model setup
model = UnremovableABM(Microbe{d}, extent, dt; rng)
for i in 1:nbacteria
	add_agent!(model)
end
##

##
adata = [:pos, :vel]
adf, mdf = run!(model, nsteps; adata)
##

## analysis
lags = 0:15:size(x,1)-1
M = msd(adf; L)
ϕ = acf(adf, :vel, lags)
##

## plotting
begin
	ts = (1:nsteps-1) .* dt
	plot(ts, M, lw=5, scale=:log10, leg=false)
	plot!(ts, t -> (30*t)^2, ls=:dash)
	plot!(ts, t -> 2d * (30^2/3)*t, ls=:dash)
	plot!(xlab="time (s)", ylab="MSD (μm²)")
end
begin
	plot((0:nsteps).*dt, t->exp(-t), lw=5)
	scatter!(lags.*dt, ϕ, xlims=(0,6), m=:x, msw=0, ms=5)
	plot!(xlab="time (s)", ylab="VACF")
end
##


## individual trajectories
gdf = groupby(adf, :id)
##
## individual acf
ϕ = acf.(collect(gdf), :vel)
plot(gdf[1].step.*dt, ϕ, xlims=(0,6), lab=false, lw=0.3, ylims=(-0.15,1.05))
plot!(gdf[1].step.*dt, mean(ϕ), lw=4, lc=:white, lab=false)
##
## individual msd
M = msd.(collect(gdf); L)
plot(ts, M, lab=false, lw=0.3)
plot!(ts, mean(M), lw=4, lc=:white, lab=false)
##