##
using MicrobeAgents, DataFrames, Statistics, Random
using Plots; Plots.gr()
include("../plot_defaults.jl")
include("../utils.jl")
rng = MersenneTwister(33)
##

## concentration field
@inline function concentration_field(pos, model)
    C₀ = model.C₀
    C₁ = model.C₁
    Lx = first(model.space.extent)
    concentration_field(pos,Lx,C₀,C₁)
end
@inline concentration_field(pos,Lx,C₀,C₁) = C₀ + (C₁-C₀)*pos[1]/Lx

@inline function concentration_gradient(pos, model)
    C₀ = model.C₀
    C₁ = model.C₁
    Lx = first(model.space.extent)
    ntuple(i -> i==1 ? concentration_gradient(pos[i],Lx,C₀,C₁) : 0.0, length(pos))
end
@inline concentration_gradient(x,Lx,C₀,C₁) = (C₁-C₀)/Lx
##

## simulation parameters
Lx, Ly = 3000, 1500
extent = (Lx, Ly) # domain size (μm)
Δt = 0.1 # timestep (s)
T = 120 # simulation time (s)
nsteps = round(Int, T/Δt)
n = 100
##

## model setup
periodic = false
C₀, C₁ = 0.0, 20.0 # μM
properties = Dict(
    :C₀ => C₀,
    :C₁ => C₁,
    :concentration_field => concentration_field,
    :concentration_gradient => concentration_gradient
)
model = UnremovableABM(Brumley{2}, extent, Δt; rng, periodic, properties)
# add n bacteria with perfect sensing
for i in 1:n
    pos = (0.0, rand(model.rng)*Ly)
    chemotactic_precision = 0
    add_agent!(pos, model; chemotactic_precision)
end
# n bacteria with noisy sensing
for i in 1:n
    pos = (0.0, rand(model.rng)*Ly)
    chemotactic_precision = 60
    add_agent!(pos, model; chemotactic_precision)
end
##

## run
adata = [:pos, :vel]
adf, mdf = run!(model, nsteps; adata)
##

## postprocessing
traj = vectorize_adf_measurement(adf, :pos)
x = first.(traj)
y = second.(traj)
isitnoisy(microbe) = microbe.chemotactic_precision > 0
noisy_idx = findall(isitnoisy.(model.agents))
perfect_idx = setdiff(eachindex(model.agents), noisy_idx)
##

## plotting
ts = unique(adf.step) .* Δt
lw = eachindex(ts) ./ length(ts) .* 3
xmesh = range(0,Lx,length=100)
ymesh = range(0,Ly,length=100)
c = concentration_field.(Iterators.product(xmesh,ymesh),Lx,C₀,C₁)
idxplot = rand(rng,axes(x,2),10)
xn = @view x[:,idxplot]
yn = @view y[:,idxplot]
col = [i ∈ noisy_idx ? 1 : 2 for i in idxplot]
heatmap(xmesh, ymesh, c', cbar=false, ratio=1, axis=false, c=:bone)
plot!(xn, yn, lab=false, lw=lw, lc=col')
scatter!(xn[end,:], yn[end,:], lab=false, m=:c, mc=col, msw=0.5, ms=6)
##

## drift
gdf = groupby(adf, :id)
vd = map(df -> driftvelocity_direction(df, (1,0)), collect(gdf))
vd_perfect = vd[perfect_idx]
vd_noisy = vd[noisy_idx]
#plot(ts, mean(vd), lw=1, lab="All", leg=:topright)
plot(ts, [mean(vd_perfect) mean(vd_noisy)], lab=["Perfect" "Noisy"], lw=1)
##

## lateral transport
adfx = transform(adf, :pos => ByRow(p -> (p[1],)) => :pos)
Mx_perfect = map(msd, collect(groupby(adfx, :id)[perfect_idx]))
Mx_noisy = map(msd, collect(groupby(adfx, :id)[noisy_idx]))
begin
    plot(
        plot(ts[2:end-1], Mx_perfect, lw=0.25, title="Perfect", ylab="MSDₓ (μm²)"),
        plot(ts[2:end-1], Mx_noisy, lw=0.25, yformatter=_->"", title="Noisy"),
        scale=:log10, legend=false,
        xlims=(0.5,1e2), ylims=(1e2,2e6), xticks=exp10.(0:2),
        xlab="time (s)"
    )
    plot!(subplot=1, ts[2:end-1], mean(Mx_perfect), lw=3, lc=:white)
    plot!(subplot=2, ts[2:end-1], mean(Mx_noisy), lw=3, lc=:white)
    plot!(subplot=1, 0.6:20, t -> 1e3*t^2, ls=:dash, lc=:cyan)
    plot!(subplot=2, 0.6:20, t -> 1e3*t^2, ls=:dash, lc=:cyan)
    plot!(subplot=1, 1:50, t -> 1e4*t, ls=:dash, lc=:orange)
    plot!(subplot=2, 1:50, t -> 1e4*t, ls=:dash, lc=:orange)
end
##