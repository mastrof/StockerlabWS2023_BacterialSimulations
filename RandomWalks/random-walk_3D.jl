##
using MicrobeAgents, DataFrames, Distributions, Random
using Plots
include("../plot_defaults.jl")
# initialize a random numbre generator for reproducibility
rng = MersenneTwister(45)
##

## simulation parameters
L = 1000
extent = (L,L,L) # domain size (μm)
dt = 0.01 # integration timestep (s)
T = 120 # simulation time (s)
nsteps = round(Int, T/dt)
nbacteria = 3
##

## model setup
# initialize the model with appropriate boundary conditions
periodic = true
model = UnremovableABM(Microbe{3}, extent, dt; rng, periodic)
# add bacteria to the model, all with the same starting position
for i in 1:nbacteria
	pos = (L/2,L/2,L/2) # initial position
	turn_rate = 1.0 # average turn rate (Hz)
	add_agent!(pos, model; turn_rate)
end
##

## simulation
# select what properties to measure during simulation
adata = [:pos]
# run the model
adf, mdf = run!(model, nsteps; adata)
##

## post processing
# extract trajectories of each bacterium
trajectories = if periodic
	MicrobeAgents.unfold(vectorize_adf_measurement(adf, :pos), L)
else
	vectorize_adf_measurement(adf, :pos)
end
x = first.(trajectories) .- L/2
y = second.(trajectories) .- L/2
z = last.(trajectories) .- L/2
##

## plotting
t = (0:nsteps) .* dt
plot(x, y, z,
	lab=false,
	xlab="x (μm)", ylab="y (μm)", zlab="z (μm)",
	ratio=1
)
##
