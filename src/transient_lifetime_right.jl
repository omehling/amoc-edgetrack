include("model_definitions.jl")

using DynamicalSystemsBase
using ChaosTools
using OrdinaryDiffEq
using Statistics
using Printf
using Dates
using HDF5
using Distributions: quantile, Gamma

println(@sprintf("Running with %d threads", Threads.nthreads()))
t1 = now()

# From CriticalTransitions.jl
function intervals_to_box(bmin::Vector, bmax::Vector)
    # Generates a box from specifying the interval limits
    intervals = []
    dim = length(bmin)
    for i in 1:dim
        push!(intervals, bmin[i]..bmax[i])
    end
    box = intervals[1]
    for i in 2:dim
        box = box × intervals[i]
    end
    box
end

# G21 model (no ice) with scaling of (x,y,z)
sys2(u,p,t) = gottwald_noice_scaled(u,p,t;stochsys=false)
scale = 5000.

# Parameters (freshwater default = 0.9)
p_default = param_gwn_default()
p_scaled = zeros(length(p_default)+1)
p_scaled[1:end-1]=p_default[:]
p_scaled[7]=0.1 # F1
p_scaled[5]=3e-4 # eps_f
p_scaled[end]=scale
# Initialize ds
u1=[2.,1.,0.,1.5,0.5]./[scale, scale, scale, 1., 1.]
u2=[1., 0., 2., 1.5, 1.125]./[scale, scale, scale, 1., 1.]
ds = ContinuousDynamicalSystem(sys2, u1, p_scaled)

# Timestepping arguments
diffeq_args = (alg = DP5(), adaptive=false, dt=7.5e-6)
t_tot = 200. # total time in time units (for attractor trajectories)
spinup = 20. # spinup time in time units (-"-)
Δt = diffeq_args.dt*100
tsol = 0:Δt:t_tot

function attractors_fw(fw; p_orig = p_scaled, t_tot = 200., spinup = 70., diffeq = NamedTuple(), Δt=1e-2)
    u1=[2.,1.,0.,1.5,0.5]./[scale, scale, scale, 1, 1]
    u2=[1., 0., 2., 1.5, 1.125]./[scale, scale, scale, 1, 1]
    p_mod = deepcopy(p_orig)
    p_mod[12] = fw
    ds_mod = ContinuousDynamicalSystem(sys2, u1, p_mod)
    
    sol_u1_mod = trajectory(ds_mod, t_tot, u1, Ttr=spinup, Δt=Δt, diffeq = diffeq) # q > 0
    sol_u2_mod = trajectory(ds_mod, t_tot, u2, Ttr=spinup, Δt=Δt, diffeq = diffeq) # q < 0
    
    return sol_u1_mod, sol_u2_mod
end

function ds_sol_bbox(ds_sol; eps = 0.)
    mins = vec(minimum(Matrix(ds_sol); dims=1))
    maxs = vec(maximum(Matrix(ds_sol); dims=1))
    return intervals_to_box(mins.-eps, maxs.+eps)
end

function calc_tau_transient(
        ds_mod, tr_init, Δt_sample, Tmax, transient_bbox, stable_bbox;
        diffeq_args = NamedTuple(), τ_break=1000*Δt_sample
    )
    integ_transient = integrator(ds_mod, tr_init; diffeq = diffeq_args)
    τ_stable = τ_transient = 0.
    t0 = integ_transient.t
    while integ_transient.t < Tmax
        step!(integ_transient, Δt_sample)
        # Check whether state is on transient "attractor"
        if (get_state(integ_transient) in transient_bbox)
            τ_transient += (integ_transient.t - t0)
        end
        # Check whether state is on stable attractor (for break criterion)
        if (get_state(integ_transient) in stable_bbox)
            τ_stable += (integ_transient.t - t0)
        end
        if τ_stable >= τ_break
            return τ_transient
        end
        # Update for next step
        t0 = integ_transient.t
    end
    return τ_transient
end

function tau_transient_ensemble(
        fw, transient_inits;
        transient_bbox = nothing, Tmax = 10000, Δt_sample = 1., t_tot_stable = 200., spinup_stable = 70.,
        p_orig = p_scaled, diffeq_args = NamedTuple(), box_eps = 0.005
    )
    if isnothing(transient_bbox)
        transient_bbox = ds_sol_bbox(transient_inits, eps = box_eps)
    end
    # Stable solution
    stable_init = [2.,1.,0.,1.5,1.5]./[scale, scale, scale, 1., 1.]
    p_mod = deepcopy(p_orig)
    p_mod[12] = fw
    ds_mod = ContinuousDynamicalSystem(sys2, u1, p_mod)

    sol_stable = trajectory(ds_mod, t_tot_stable, stable_init, Ttr=spinup_stable, Δt=Δt_sample, diffeq = diffeq_args)
    stable_bbox = ds_sol_bbox(sol_stable, eps = box_eps)
    # Transient solution
    τ_transients = zeros(size(transient_inits)[1])
    Threads.@threads for i = 1:size(transient_inits)[1]
        # Integrate
        tr_init = transient_inits[i,:]
        # Lifetime of transient state
        τ_transient = calc_tau_transient(ds_mod, tr_init, Δt_sample, Tmax, transient_bbox, stable_bbox, diffeq_args = diffeq_args, τ_break=100*Δt_sample)
        
        τ_transients[i] = τ_transient
    end
    return τ_transients
end

fw_crit = 0.926
u1_fw_crit, u2_fw_crit = attractors_fw(fw_crit, spinup=50, t_tot = 1000, diffeq = diffeq_args)

transient_inits = u1_fw_crit[1:1000:end-1,:] # 100 samples
println("initializing ", size(transient_inits)[1], " trajectories on ghost attractor...")

fw_range = .9304:0.0002:.934
τ_transients_fw = zeros(length(fw_range), size(transient_inits)[1])

for (i, fw_transient) in enumerate(fw_range)
    Δt_sample = 0.01
    τ_transients = tau_transient_ensemble(
        fw_transient, transient_inits,
        transient_bbox = ds_sol_bbox(u1_fw_crit, eps = 0.015), Tmax = 20000, Δt_sample = Δt_sample,
        diffeq_args = diffeq_args
    );
    τ_transients_fw[i, :] = τ_transients[:]

    # Calculate Statistics
    tau_est = mean(τ_transients)
    tau_CIs = mean(τ_transients)*length(τ_transients) ./ quantile.(Gamma(length(τ_transients), 1), [0.025, 0.975])
    println("fw = ", fw_transient, ", Lifetime: ", tau_est, " ", tau_CIs)
end

### Write output to HDF5 file
h5open("transients_lifetimes_right.h5", "w") do file
    write(file, "fw", convert(Array{Float32}, fw_range))
    write(file, "tau", convert(Array{Float32}, τ_transients_fw))
end

t2 = now()
println("Sampling ICs finished successfully. Wall-clock time:")
println(canonicalize(t2 - t1))
