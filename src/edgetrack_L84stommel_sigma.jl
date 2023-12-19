include("edgetrack.jl")
include("gottwald_model_noice.jl")

using Plots
using Statistics
using ProgressMeter
using Random
using Printf
using HDF5
using Dates

write_output = true

function output_matrix(dataset; dtype=Float32)
    out_matrix = zeros(dtype, size(dataset))
    for (i, state) in enumerate(dataset)
        out_matrix[i,:] = convert(Array{dtype}, collect(state))
    end
    return out_matrix
end

# G21 model with scaling of (x,y,z)
sys2(u,p,t) = gottwald_noice_scaled(u,p,t;stochsys=false)
scale = 5000.
# Parameters
sigma_0 = tryparse(Float64,ARGS[1])
if isnothing(sigma_0) == true
    throw(ArgumentError(
        @sprintf("Value for sigma_0 provided to script (%s) cannot be converted to Float64", ARGS[1])
    ))
end
p_default = param_gwn_default()
p_scaled = zeros(length(p_default)+1)
p_scaled[1:end-1]=p_default[:]
p_scaled[7]=0.1 # F1
p_scaled[5]=3e-4 # epsilon_f
p_scaled[12]=sigma_0  # freshwater flux, cmd line argument
p_scaled[end]=scale
# Initialize ds
u1=[2.,1.,0.,1.5,0.5]./[scale, scale, scale, 1., 1.]
u2=[1., 0., 2., 1.5, 1.125]./[scale, scale, scale, 1., 1.]
ds = ContinuousDynamicalSystem(sys2, u1, p_scaled)

# Timestepping arguments
diffeq_args = (alg = DP5(), adaptive=false, dt=7.5e-6)
t_tot = 100. # total time in time units (for attractor trajectories)
spinup = 150. # spinup time in time units (-"-)
Δt = diffeq_args.dt*10
tsol = 0:Δt:t_tot

# Find attractors by forward simulation
sol_u1 = trajectory(ds, t_tot, u1, Ttr=spinup, Δt=Δt, diffeq = diffeq_args);
sol_u2 = trajectory(ds, t_tot, u2, Ttr=spinup, Δt=Δt, diffeq = diffeq_args);
attrs = Dict(1 =>Dataset([mean(sol_u1.data)]), 2 => Dataset([mean(sol_u2.data)]))

lya_a1 = lyapunov(ds, t_tot; u0 = u1, Ttr=spinup, upper_threshold=1e-7, lower_threshold=1e-11, Δt=Δt, diffeq = diffeq_args)
lya_a2 = lyapunov(ds, t_tot; u0 = u2, Ttr=spinup, upper_threshold=1e-7, lower_threshold=1e-11, Δt=Δt, diffeq = diffeq_args)

println("Attractor 1: ", Matrix(attrs[1]))
println("(Lyapunov = ", lya_a1, ")")
println("Attractor 2: ", Matrix(attrs[2]))
println("(Lyapunov = ", lya_a2, ")")

t1 = now()

Niter = 2000 # iterations for edge tracking algorithm
edge, track1, track2 = edgetracking(
                        ds, u1, u2, attrs,
                        eps1 = 0.0025, eps2 = 0.004, dt=diffeq_args.dt, ϵ_mapper=0.01,
                        diffeq=diffeq_args, maxit=Niter, mx_chk_lost=round(Int, 10000/diffeq_args.dt),
                        output_level=2
                        )

tedge = 0:diffeq_args.dt:(size(edge)[1]-1)*diffeq_args.dt

# Write output to HDF5 file
if write_output
    h5open(@sprintf("edgetrack_L84stommel_sig%04d.h5", sigma_0*1000), "w") do file
        write(file, "t", convert(Array{Float32}, tedge))
        write(file, "u_edge", output_matrix(edge))
        write(file, "u_track1", output_matrix(track1))
        write(file, "u_track2", output_matrix(track2))
    end
end

t2 = now()
print("Edge tracking finished successfully. Wall-clock time: ")
println(canonicalize(t2 - t1))
