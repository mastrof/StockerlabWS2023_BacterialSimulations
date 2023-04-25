##
using MicrobeAgents, DataFrames, Statistics, BubbleBath, Distributions
using Plots; Plots.gr()
include("../plot_defaults.jl")
include("../utils.jl")
##

## function to draw circles
function circle(x,y,r,n=500)
    θ = LinRange(0, 2π, n)
    x .+ r.*sin.(θ), y .+ r.*cos.(θ)
end # function
##

## create an environment packed with spheres
Lx, Ly = 1000.0, 500.0
extent = (Lx, Ly)
Rmin, R₀, Rmax = 20.0, 50.0, 80.0
radius_pdf = Truncated(Exponential(R₀), Rmin, Rmax)
ϕ = 0.5 # packing fraction
min_distance = 10
spheres = bubblebath(radius_pdf, ϕ, extent; min_distance)
##

## plot spheres
plot(
    xlims=(0,Lx), ylims=(0,Ly),
    ratio=1, legend=false, grid=false,
)
for s in spheres
    plot!(circle(s.pos..., s.radius), seriestype=:shape, fc=1)
end
plot!()
##

## create walkmap
wm = walkmap(spheres, extent, 0.5)
Nx, Ny = size(wm)
contourf(wm', ratio=1, xlims=(0,Lx), ylims=(0,Ly), levels=1)
##

## setup model with pathfinder
dt = 0.1
model = UnremovableABM(Microbe{2}, extent, dt; periodic=false)
pathfinder!(model, wm)
n = 50
accessible_idx = findall(wm) # list of all accessible positions
for i in 1:n
    ix, iy = Tuple(rand(accessible_idx))
    pos = (Lx*(ix-1)/Nx, Ly*(iy-1)/Ny)
    add_agent!(pos, model; rotational_diffusivity=0.2)
end
##

## run
adata = [:pos]
nsteps = 500
adf, mdf = run!(model, microbe_pathfinder_step!, model.update!, nsteps; adata)
#adf, mdf = run!(model, nsteps; adata)
##

## get trajectories
traj = vectorize_adf_measurement(adf, :pos)
x = first.(traj)
y = last.(traj)
##

## plot trajectories
xmesh = range(0,Lx,length=Nx)
ymesh = range(0,Ly,length=Ny)
contourf(xmesh, ymesh, wm',
    xlims=(0,Lx), ylims=(0,Ly), ticks=false,
    cbar=false, c=:bone, ratio=1
)
plot!(x, y, lab=false, lw=1)
##