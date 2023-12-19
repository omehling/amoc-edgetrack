# Minimal test of the Stommel-L84 model *without* scaling of atmospheric variables
#
# This script defines the Stommel-L84 model as used in Mehling et al. (2023),
# integrates a long trajectory of each attractor, and plots the attractors
# in T-S subspace
#
include("../src/model_definitions.jl")

using DynamicalSystemsBase
using OrdinaryDiffEq
using Plots

sys(u,p,t) = gottwald_noice(u,p,t;stochsys=false)
# Parameters
sigma_0 = 0.8 # freshwater flux (control parameter)
p_default = param_gwn_default()
p_default[7]=0.1 # F_1
p_default[5]=3e-4 # epsilon_f (timescale separation)
p_default[12]=sigma_0
# Initialize ds
u1=[2.,1.,0.,1.5,0.5] # for q > 0 ("on")
u2=[1., 0., 2., 1.5, 1.125] # for q < 0 ("off")
ds_default = ContinuousDynamicalSystem(sys, u1, p_default)

# Timestepping arguments
diffeq_args = (alg = DP5(), adaptive=false, dt=7.5e-6);
t_tot = 200.
spinup = 70.
Δt_output = diffeq_args.dt*100
#tsol = 0:Δt_output:t_tot;

sol_u1 = trajectory(ds_default, t_tot, u1, Ttr=spinup, Δt=Δt_output, diffeq = diffeq_args) # q > 0
sol_u2 = trajectory(ds_default, t_tot, u2, Ttr=spinup, Δt=Δt_output, diffeq = diffeq_args) # q < 0

p = plot(xlabel="T", ylabel="S")
plot!(p, sol_u1[:,4], sol_u1[:,5], color=:2, label="q > 0")
plot!(p, sol_u2[:,4], sol_u2[:,5], color=:1, label="q < 0")
savefig(p, "StommelL84_attractors_test_unscaled.png")