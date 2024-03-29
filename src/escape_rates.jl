include("model_definitions.jl")

using DynamicalSystemsBase
using ChaosTools
using OrdinaryDiffEq
using Statistics
using Printf
using Dates
using HDF5
using Random
using LinearAlgebra: norm
using DataFrames, CSV

function load_edge(fw; subdir = "/home/omehling/work/reading/edge_tracking/")
    filename = @sprintf("edgetrack_L84stommel_sig%04d.h5", fw*1000)
    
    # Load computed edge state from hdf5 file -> generated by edgetrack_L84stommel_sigma.jl
    u_edge = h5open(subdir * filename, "r") do file
        read(file, "u_edge")
    end
    t_edge = h5open(subdir * filename, "r") do file
        read(file, "t")
    end
    return t_edge, u_edge
end

function random_Q0(ds, e; dtype=nothing)
    D = dimension(ds)
    T = (isnothing(dtype)) ? eltype(ds) : dtype
    Q0 = randn(Random.GLOBAL_RNG, SVector{D, T})
    Q0 = e * Q0 / norm(Q0)
end

function idx_ranges(ts, lim1, lim2; descending=false)
    tfirst=length(ts)
    tlast=length(ts)
    for t=1:length(ts)
        if descending
            if (ts[t]<lim1) & (t<tfirst)
                tfirst = t
            end
            if (ts[t]<lim2) & (t<tlast)
                tlast = t
            end
        else
            if (ts[t]>lim1) & (t<tfirst)
                tfirst = t
            end
            if (ts[t]>lim2) & (t<tlast)
                tlast = t
            end
        end
    end
    tfirst, tlast
end

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

function escape_rates(ds::DynamicalSystem, box::IntervalBox, start_samples, S, eps_init, Tsim, pmin, Δt;
                     diffeq=NamedTuple(), output="aggregated")
    tperturb = 0:Δt:Tsim
    is_inbox = zeros(Bool, size(start_samples)[1]*S, length(tperturb))
    edge_inits = []
    tperturb = []
    p_inbox = []
    # Generate initial conditions (start_samples + perturbations)
    for i = 1:size(start_samples)[1]
        start_sample = start_samples[i,:]
        for j = 1:S
            append!(edge_inits, [start_sample .+ random_Q0(ds, eps_init, dtype=eltype(start_sample))])
        end
    end
    
    pinteg = parallel_integrator(ds, edge_inits; diffeq = diffeq)
    pinteg_dim = dimension(pinteg)
    println("Initializing ", pinteg_dim, " perturbed states for escape rate calculation...")
    append!(p_inbox, sum([get_state(pinteg, i) in box for i=1:pinteg_dim])/pinteg_dim)
    append!(tperturb, pinteg.t)
    t1 = now()
    while (pinteg.t < Tsim) && (p_inbox[end] > pmin)
        step!(pinteg, Δt)
        append!(p_inbox, sum([get_state(pinteg, i) in box for i=1:pinteg_dim])/pinteg_dim)
        append!(tperturb, pinteg.t)
    end
    t2 = now()
    println("Escape rates calculated successfully (", canonicalize(t2 - t1), ")")
    return tperturb, p_inbox
end

# G21 model (no ice) with scaling of (x,y,z)
sys2(u,p,t) = gottwald_noice_scaled(u,p,t;stochsys=false)
scale = 5000.
# Parameters
sigma_0 = 0.8
p_default = param_gwn_default()
p_scaled = zeros(length(p_default)+1)
p_scaled[1:end-1]=p_default[:]
p_scaled[7]=0.1 # F1
p_scaled[5]=3e-4 # epsilon_f
p_scaled[12]=sigma_0  # default freshwater flux
p_scaled[end]=scale
# Initialize ds
u1=[2.,1.,0.,1.5,0.5]./[scale, scale, scale, 1., 1.]
u2=[1., 0., 2., 1.5, 1.125]./[scale, scale, scale, 1., 1.]
ds_default = ContinuousDynamicalSystem(sys2, u1, p_scaled)

# Timestepping arguments
diffeq_args = (alg = DP5(), adaptive=false, dt=7.5e-6)
spinup_steps = 300000

function ds_fw(fw)
    p_fw = deepcopy(p_scaled)
    p_fw[12] = fw
    return ContinuousDynamicalSystem(sys2, u1, p_fw)
end

is_inbox_dict = Dict()
for fw_sel in [0.8:0.01:0.92; 0.921:0.001:0.926]
    println(fw_sel)
    t_edge_sel, u_edge_sel = load_edge(fw_sel)
    ds_mod = ds_fw(fw_sel)

    starts_mod = spinup_steps:10000:size(u_edge_sel)[1]
    edgestate_box_mod = intervals_to_box(
        vec(minimum(Matrix(u_edge_sel[spinup_steps:end,:]); dims=1)),
        vec(maximum(Matrix(u_edge_sel[spinup_steps:end,:]); dims=1))
    )

    t_perturb, is_inbox = escape_rates(ds_mod,
        edgestate_box_mod, u_edge_sel[starts_mod,:],
        5, 1e-9, 20., 0.007, diffeq_args.dt*100, diffeq=diffeq_args
    )
    println("Total steps: ", length(is_inbox), ", time: ", t_perturb[end])
    is_inbox_dict[@sprintf("%.03f", fw_sel)] = is_inbox
end

# pad with missing values and export to CSV via dataframe
maxlength = maximum([length(is_inbox_dict[key]) for key in keys(is_inbox_dict)])
is_inbox_dict_pad = Dict()
for key in keys(is_inbox_dict)
    is_inbox_dict_pad[key] = [is_inbox_dict[key]; [missing for _ in 1:(maxlength - length(is_inbox_dict[key]))]]
end
CSV.write("escape_rates.csv", DataFrame(is_inbox_dict_pad))
