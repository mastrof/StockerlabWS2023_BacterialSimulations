
##
using MicrobeAgents
using DelimitedFiles, DataFrames, Distributions, StatsBase, LinearAlgebra
using Plots; Plots.gr()
include("../plot_defaults.jl")
include("../utils.jl")
##

## import and normalize angular distribution from Berg&Brown 1972
θraw, wraw = eachcol(readdlm("ecoli_angularpdf.txt"; comment_char='#', comments=true)[:,1:2])
θ = deg2rad.(θraw)
A = sum(diff(θ).*midpoints(wraw))
w = wraw./A
scatter(θ, w)
##

## extend pdf to [-π,π] by mirroring
ϕ = [-reverse(θ); θ]
p = [reverse(w); w]
A = sum(diff(ϕ).*(midpoints(p)))
p = p./A
scatter(θ, w, lab="Original")
scatter!(ϕ, p, lab="Extended")
##

## find a mixture model that matches the empirical pdf
scatter(θ,w)
ξ, ω, α = 0.45, 0.93, 3.4
#plot!(0:π/50:π, x -> pdf(SkewNormal(ξ, ω, α), x), lw=3)

scatter(ϕ, p, lab="Empirical")
MixMod = MixtureModel([SkewNormal(-ξ,ω,-α), SkewNormal(ξ,ω,α)])
plot!(-π:π/50:π, x -> pdf(MixMod, x), lw=3, lab="Model")

# try to sample from the model pdf
ϕ_sampled = rand(MixMod, 10_000)
histogram!(ϕ_sampled, normalize=:pdf, alpha=0.4, lab="Samples", lw=0.2)
##

## model setup
L = 1000.0
extent = (L,L)
dt = 0.1
model = UnremovableABM(Microbe{2}, extent, dt)
n = 500
for _ in 1:n
    motility = RunTumble(polar=MixMod)
    add_agent!(model; motility)
    add_agent!(model)
end
adata = [:pos, :vel]
##

## run
nsteps = 3000
adf, mdf = run!(model, nsteps; adata)
##

## postprocessing
traj = MicrobeAgents.unfold(vectorize_adf_measurement(adf, :pos), L)
x = first.(traj)
y = last.(traj)
motidx = findall(m -> m.motility.polar == MixMod, model.agents)
##

## plot trajectories
plotidx = rand(eachindex(model.agents), 8)
lc = [i ∈ motidx ? 1 : 2 for i in plotidx]
xn = @view x[:,plotidx]
yn = @view y[:,plotidx]
plot(xn, yn, lab=false, lw=1, axis=false, ratio=1, lc=lc')
##

## angular distributions from simulation
vel = vectorize_adf_measurement(adf, :vel)
ψ = [safe_acos(dot(vel[i,j],vel[i-1,j])) for i in axes(vel,1)[2:end], j in axes(vel,2)]
ψ₁ = filter(s -> s≠0, ψ[:,motidx])
ψ₂ = filter(s -> s≠0, ψ[:,setdiff(1:2n,motidx)])
histogram(ψ₁, normalize=:pdf, bins=50, alpha=0.8)
histogram!(ψ₂, normalize=:pdf, bins=50, alpha=0.5)
scatter!(θ, w)
vline!([mean(ψ₁)], ls=:dash, lc=1, lab=false ,lw=3)
vline!([mean(ψ₂)], ls=:dash, lc=2, lab=false ,lw=3)
##

## compare msd between populations
transform!(adf, :id => ByRow(i -> i∈motidx) => :motidx)
gdf = groupby(adf, :motidx)
M₁ = MicrobeAgents.msd(gdf[(true,)]; L)
M₂ = MicrobeAgents.msd(gdf[(false,)]; L)
ts = eachindex(M₁) .* dt
plot(ts, [M₁, M₂],
    lw=4, lab=["Custom PDF" "Uniform"],
    legend=:topleft, scale=:log10,
    xlab="time (s)", ylab="MSD (μm²)",
)
##

## ratio between the two MSDs
plot(ts, M₁./M₂, lab=false, lw=4, xlab="time (s)", ylab="M₁ / M₂")
##