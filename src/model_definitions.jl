using StaticArrays

# Smooth absolute value function
function smoothabs(x, xi=10000)
    x*tanh(x*xi)
end

# L84 model equations
function L84(u,p,t)
    x, y, z = u
    a, b, F, G = p[1]

    Δ = y^2 + z^2
    dx = -Δ - a*(x - F)
    dy = x*y - b*x*z - (y - G)
    dz = b*x*y + x*z - z

    SA[dx, dy, dz]
end

# Stommel model equations following the Gottwald (2021) formulation
function stommel_gottwald(u,p,t)
    T, S = u
    μ, ϵa, θ0, σ0 = p[1]

    T_surf = θ0
    S_surf = σ0
    dT = -1/ϵa*(T - T_surf) - T - μ*smoothabs(S-T)*T
    dS = S_surf - S - μ*smoothabs(S-T)*S

    SA[dT, dS]
end

# Gottwald model functions without sea ice
function gottwald_noice(u,p,t; stochsys=true)
    x, y, z, T, S = u
    if stochsys
        pvec = p[1]
    else
        pvec = p
    end
    a, b, μ, ϵa, ϵf, F0, F1, G0, G1, θ0, θ1, σ0, σ1, x_mean, Δ_mean = pvec

    Δ = y^2 + z^2
    T_surf = θ0 + θ1 * (x - x_mean)/(sqrt(ϵf))
    S_surf = σ0 + σ1 * (Δ - Δ_mean)/(sqrt(ϵf))
    dx = 1/ϵf*(-Δ - a*(x - F0 - F1*T))
    dy = 1/ϵf*(x*y - b*x*z - (y - G0 + G1*T))
    dz = 1/ϵf*(b*x*y + x*z - z)
    dT = -1/ϵa*(T - T_surf) - T - μ*smoothabs(S-T)*T
    dS = S_surf - S - μ*smoothabs(S-T)*S

    SA[dx, dy, dz, dT, dS]
end

# Rescaled (x/y/z) for edge tracking
function gottwald_noice_scaled(u,p,t; stochsys=true)
    x, y, z, T, S = u
    if stochsys
        pvec = p[1]
    else
        pvec = p
    end
    a, b, μ, ϵa, ϵf, F0, F1, G0, G1, θ0, θ1, σ0, σ1, x_mean, Δ_mean, s = pvec

    Δ = y^2 + z^2
    T_surf = θ0 + θ1 * (s*x - x_mean)/(sqrt(ϵf))
    S_surf = σ0 + σ1 * (s^2*Δ - Δ_mean)/(sqrt(ϵf))
    dx = 1/ϵf*(-s*Δ - a/s*(s*x - F0 - F1*T))
    dy = 1/ϵf*(s*x*y - s*b*x*z - (s*y - G0 + G1*T)/s)
    dz = 1/ϵf*(s*b*x*y + s*x*z - z)
    dT = -1/ϵa*(T - T_surf) - T - μ*smoothabs(S-T)*T
    dS = S_surf - S - μ*smoothabs(S-T)*S

    SA[dx, dy, dz, dT, dS]
end

# Gottwald model functions without sea ice (inplace)
function gottwald_noice!(du,u,p,t)
    x, y, z, T, S = u
    a, b, μ, ϵa, ϵf, F0, F1, G0, G1, θ0, θ1, σ0, σ1, x_mean, Δ_mean = p[1]

    Δ = y^2 + z^2
    T_surf = θ0 + θ1 * (x - x_mean)/(sqrt(ϵf))
    S_surf = σ0 + σ1 * (Δ - Δ_mean)/(sqrt(ϵf))
    du[1] = 1/ϵf*(-Δ - a*(x - F0 - F1*T)) # x
    du[2] = 1/ϵf*(x*y - b*x*z - (y - G0 + G1*T)) # y
    du[3] = 1/ϵf*(b*x*y + x*z - z) # z
    du[4] = -1/ϵa*(T - T_surf) - T - μ*smoothabs(S-T)*T # T
    du[5] = S_surf - S - μ*smoothabs(S-T)*S # S
end

# Gottwald model functions without sea ice, but with seasonal cycle
function gottwald_noice_seasons(u,p,t)
    x, y, z, T, S = u
    a, b, μ, ϵa, ϵf, F0, F1, G0, G1, θ0, θ1, σ0, σ1, x_mean, Δ_mean, td, F2, G2 = p[1]

    ω = 2.0*pi*td # annual frequency
    Δ = y^2 + z^2
    T_surf = θ0 + θ1 * (x - x_mean)/(sqrt(ϵf))
    S_surf = σ0 + σ1 * (Δ - Δ_mean)/(sqrt(ϵf))
    dx = 1/ϵf*(-Δ - a*(x - F0 - F1*T - F2*cos(ω * t)))
    dy = 1/ϵf*(x*y - b*x*z - (y - G0 + G1*T - G2*cos(ω * t)))
    dz = 1/ϵf*(b*x*y + x*z - z)
    dT = -1/ϵa*(T - T_surf) - T - μ*smoothabs(S-T)*T
    dS = S_surf - S - μ*smoothabs(S-T)*S

    SA[dx, dy, dz, dT, dS]
end

# Default parameter values
function param_gwn_default(param_l84=nothing)
    if !isnothing(param_l84)
        # Long run of L84 model to determine mean and std of coupling quantities
        model_l84 = StochSystem(L84, param_l84, 3, 0.016)
        sim_l84 = relax(model_l84, [-1.5, 1.5, 1.5], tspan=(0.,5000.), solver=Tsit5(), reltol=1e-10, abstol=1e-10, saveat=range(0,5000,500001))
        x_mean = mean(sim_l84[1,50000:end])
        Δ_mean = mean(sim_l84[2,50000:end].^2 + sim_l84[3,50000:end].^2)
        x_std = std(sim_l84[1,50000:end])
        Δ_std = std(sim_l84[2,50000:end].^2 + sim_l84[3,50000:end].^2)
    else
        param_l84 = [0.25, 4., 8., 1.] # default L84 params [a, b, F0, G0]
        x_std = 0.513; Δ_std = 1.071; x_mean = 1.0147; Δ_mean = 1.7463 # long-run results from default L84
    end
    # Unpack L84 params
    a, b, F0, G0 = param_l84
    # Other default model params
    μ = 7.5; θ0 = 1.; σ0 = 0.9 # Stommel params
    ϵa = 0.34; ϵf = .0083 # timescales
    F1 = 0.; G1 = 0.; # coupling params
    perturb_scaling = 0.01 # coupling strength of L84 to Stommel

    # Coupling params from L84 run/default
    θ1 = min(θ0, perturb_scaling/x_std)
    σ1 = min(σ0, perturb_scaling/Δ_std)

    [a, b, μ, ϵa, ϵf, F0, F1, G0, G1, θ0, θ1, σ0, σ1, x_mean, Δ_mean]
end