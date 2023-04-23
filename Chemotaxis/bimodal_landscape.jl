##
using MicrobeAgents, DataFrames, Statistics, Random
using Plots; Plots.gr()
include("../plot_defaults.jl")
include("../utils.jl")
rng = MersenneTwister(33)
##

## concentration field
gausspeak(x,A,x₀,σ) = A*exp(-(x-x₀)^2/(2σ^2))
∇gausspeak(x,A,x₀,σ) = -(x-x₀)/σ^2 * gausspeak(x,A,x₀,σ)

function concentration_field(pos, model)
    A₁ = model.A₁
    A₂ = model.A₂
    x₁ = model.x₁
    x₂ = model.x₂
    σ₁ = model.σ₁
    σ₂ = model.σ₂
    x = first(pos)
    concentration_field(x,A₁,x₁,σ₁,A₂,x₂,σ₂)
end
concentration_field(x,A₁,x₁,σ₁,A₂,x₂,σ₂) = gausspeak(x,A₁,x₁,σ₁)+gausspeak(x,A₂,x₂,σ₂)

function concentration_gradient(pos, model)
    A₁ = model.A₁
    A₂ = model.A₂
    x₁ = model.x₁
    x₂ = model.x₂
    σ₁ = model.σ₁
    σ₂ = model.σ₂
    x = first(pos)
    (concentration_gradient(x,A₁,x₁,σ₁,A₂,x₂,σ₂),)
end
concentration_gradient(x,A₁,x₁,σ₁,A₂,x₂,σ₂) = ∇gausspeak(x,A₁,x₁,σ₁)+∇gausspeak(x,A₂,x₂,σ₂)
##

## simulation parameters
L = 400
extent = (L,) # domain size (μm)
Δt = 0.1 # timestep (s)
T = 300 # simulation time (s)
nsteps = round(Int, T/Δt)
n = 1000
##

## model setup
periodic = false
A₁, A₂ = 1.0, 5.0 # μM
x₁, x₂ = 100, 300 # μm
σ₁, σ₂ = 40, 40 # μm
properties = Dict(
    :A₁ => A₁, :x₁ => x₁, :σ₁ => σ₁,
    :A₂ => A₂, :x₂ => x₂, :σ₂ => σ₂,
    :concentration_field => concentration_field,
    :concentration_gradient => concentration_gradient
)
model = UnremovableABM(Brumley{1}, extent, Δt; rng, periodic, properties)
# add n bacteria with perfect sensing
for i in 1:n
    pos = (0.0,)
    chemotactic_precision = 0
    add_agent!(pos, model; chemotactic_precision)
end
# n bacteria with noisy sensing
for i in 1:n
    pos = (0.0,)
    chemotactic_precision = 20
    add_agent!(pos, model; chemotactic_precision)
end
##

## plot landscape
xmesh = range(0,L,length=500)
c = concentration_field.(xmesh,A₁,x₁,σ₁,A₂,x₂,σ₂)
plot(xmesh, c, lc=:white, lw=4, lab=false, xlab="x (μm)", ylab="C (μM)")
##

## run
adata = [:pos]
when(model,s) = s==nsteps
adf, mdf = run!(model, nsteps; adata, when)
##

## postprocessing
traj = vectorize_adf_measurement(adf, :pos)
x = first.(traj)
isitnoisy(microbe) = microbe.chemotactic_precision > 0
noisy_idx = findall(isitnoisy.(model.agents))
perfect_idx = setdiff(eachindex(model.agents), noisy_idx)
##

## histograms
transform!(adf, :pos => ByRow(p -> p[1]) => :pos)
xmesh = range(0,L,length=500)
c = concentration_field.(xmesh,A₁,x₁,σ₁,A₂,x₂,σ₂)
plot(xmesh, c./20A₂, lc=:white, lw=4, lab=false, xlab="x (μm)", yticks=false)
histogram!(adf[adf.step.==nsteps, :].pos[perfect_idx],
    lw=0, normalize=:pdf, bins=range(0,L,length=50), alpha=0.75,
    lab="Perfect"
)
histogram!(adf[adf.step.==nsteps, :].pos[noisy_idx],
    lw=0, normalize=:pdf, bins=range(0,L,length=50), alpha=0.75,
    lab="Noisy"
)
##