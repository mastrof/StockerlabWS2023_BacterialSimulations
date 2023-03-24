### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ c4e34e92-ca6c-11ed-0b0b-8b5939ad16c4
import Pkg; Pkg.activate("/home/riccardo/Science/StockerlabWS2023_BacterialSimulations")

# ╔═╡ 37a93738-3379-4ecd-9f1f-21772c529aa6
using MicrobeAgents, Plots, Random

# ╔═╡ 2600d3f9-907a-4d9e-be7d-ff838f4082b9
# initialize a random numbre generator for reproducibility
rng = MersenneTwister(rand(UInt8))

# ╔═╡ b573b317-4a36-47bb-a5f8-b856778265de
# define all the parameters for the simulation
begin
	L = 1000
	extent = (L,) # domain size (μm)
	dt = 0.01 # integration timestep (s)
	T = 600 # simulation time (s)
	nsteps = round(Int, T/dt)
	nbacteria = 5
end

# ╔═╡ 2c037456-f769-42b3-86eb-2bbdbde92b35
begin
	# initialize the model with appropriate boundary conditions
	periodic = true
	model = ABM(Microbe{1}, extent, dt; rng, periodic)
	# add bacteria to the model, all with the same starting position
	for i in 1:nbacteria
		pos = (L/2,) # initial position
		turn_rate = 1.0 # average turn rate (Hz)
		add_agent!(pos, model; turn_rate)
	end
	# select what properties to measure during simulation
	adata = [:pos]
end

# ╔═╡ 9d46bee2-21a4-42ad-a0c8-dd6824e7bea2
# run the model
adf, mdf = run!(model, nsteps; adata)

# ╔═╡ ac8f4cb8-296f-4416-a95e-84f63104a2af
# extract trajectories of each bacterium
begin
	trajectories = if periodic
		MicrobeAgents.unfold(vectorize_adf_measurement(adf, :pos), L)
	else
		vectorize_adf_measurement(adf, :pos)
	end
	x = first.(trajectories) .- L/2
end

# ╔═╡ 14227ed0-e0c0-4d59-9b0b-a3cdb3cedf0f
begin
	t = (0:nsteps) .* dt
	plot(t, x, lab=false, xlab="time (s)", ylab="displacement (μm)")
end

# ╔═╡ 23f188ae-67b3-4929-8276-daccfd73a74b
md"""
# Tasks
- modify boundary conditions
- modify the average reorientation rate (`turn_rate`)
- increase the simulation time and the number of bacteria
"""

# ╔═╡ Cell order:
# ╠═c4e34e92-ca6c-11ed-0b0b-8b5939ad16c4
# ╠═37a93738-3379-4ecd-9f1f-21772c529aa6
# ╠═2600d3f9-907a-4d9e-be7d-ff838f4082b9
# ╠═b573b317-4a36-47bb-a5f8-b856778265de
# ╠═2c037456-f769-42b3-86eb-2bbdbde92b35
# ╠═9d46bee2-21a4-42ad-a0c8-dd6824e7bea2
# ╠═ac8f4cb8-296f-4416-a95e-84f63104a2af
# ╠═14227ed0-e0c0-4d59-9b0b-a3cdb3cedf0f
# ╟─23f188ae-67b3-4929-8276-daccfd73a74b
