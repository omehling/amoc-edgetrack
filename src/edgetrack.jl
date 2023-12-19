using DynamicalSystemsBase
using DelayEmbeddings: Dataset
using ChaosTools
using OrdinaryDiffEq
using LinearAlgebra: norm

inittest_default(D) = (state1, d0) -> state1 .+ d0/sqrt(D)

"""
    edgetracking(ds, u1, u2, attractors;
        eps1=1e-9, eps2=1e-8, converge=1e-5,
        dt=0.01, ϵ_mapper=0.1, diffeq = (;alg = Vern9()), maxit=100,
        verbose=true, output_level=0, kwargs...
    )

Apply the edge tracking algorithm as described in Skufca et al. (2006) and Lucarini & Bodai
(2017) to the dynamical system `ds` starting from two system states `u1` and `u2` which must
belong to different attractors.

`attractors` is a dictionary that maps a label onto an `N`-dimensional vector, where `N` is
the dimension of `ds`. Chaotic attractors must be approximated by, e.g., their mean.

IMPORTANT: this function is superseded by the `edgetracking` function
           in https://github.com/reykboerner/Attractors.jl/blob/edgetracking/src/boundaries/edgetracking.jl
           Please use this improved implementation for new code.

# Keyword arguments
* `eps1 = 1e-9`: The maximum (Euclidean) distance between two states after bisection,
                 before the start of an edge tracking step
* `eps2 = 1e-8`: The maximum (Euclidean) distance between two states at the end of an edge tracking step
* `converge = 1e-5`: Maximum distance between two successive edge tracking steps that terminates the algorithm.
                     Should not be used (i.e., set to very small values) for chaotic or periodic attractors. Use `maxit` instead.
* `dt = 0.01`: Time interval between two successive steps of the integrator. Not necessarily identical
               to the numerical timestep of the ODE system, see `DynamicalSystems.step!`.
* `ϵ_mapper = 0.01`: passed to `AttractorsViaProximity`
* `diffeq`: tuple of arguments passed to the ODE solver, such as the algorithm and timestepping arguments
* `maxit = 100`: Maximum number of edge tracking steps
* `verbose = true`: Print diagnostic output at every iteration
* `output_level = 0`: Switch between different outputs.
                      0 = only return final state vector
                      1 = return edge and bracketing trajectory after every iteration
                      2 = return edge and bracketing trajectory every Δt time units
* `kwargs...`: Keyword arguments passed on to `AttractorsViaProximity`
"""
function edgetracking(ds::DynamicalSystem, u1, u2, attractors::Dict;
    eps1=1e-9,
    eps2=1e-8,
    converge=1e-5,
    dt=0.01,
    ϵ_mapper=0.1,
    diffeq = (;alg = Vern9()),
    maxit=100,
    verbose=true,
    output_level=0,
    kwargs...)
    
    pinteg = parallel_integrator(ds, [u1, u2]; diffeq = diffeq)
    mapper = AttractorsViaProximity(ds, attractors, ϵ_mapper; diffeq = diffeq, Δt = dt, kwargs...)

    track_edge(pinteg, mapper;
        eps1=eps1, eps2=eps2, converge=converge, dt=dt, maxit=maxit,
        verbose=verbose, calc_lyapunovs=false, output_level=output_level)
end;


"""
    track_edge(pinteg, mapper;
        eps1=1e-9, eps2=1e-8, converge=1e-5, dt=0.01, maxit=100,
        verbose=false, output_level=0,
        calc_lyapunovs=false, lyapunov_spinup=200, d0=1e-9, lower_threshold=1e-11, upper_threshold=1e-7
    )

Main function for the edge tracking algorithm (Battelino et al. 1988, Skufca et al. 2006).

This function is called from `edgetracking`, see its documentation for an explanation of parameters.

Note: `calc_lyapunovs` and following keyword arguments are experimental,
      are NOT scientifically validated
      and do NOT reproduce the results from Mehling et al. 2023, arXiv:2308.16251
"""
function track_edge(pinteg, mapper::AttractorMapper;
    eps1=1e-9, eps2=1e-8, converge=1e-5, dt=0.01, maxit=100,
    verbose=false, output_level=0,
    calc_lyapunovs=false, lyapunov_spinup=200, d0=1e-9, lower_threshold=1e-11, upper_threshold=1e-7)
    
    lyapunov_spinup += 1
    verbose && println("=== Starting edge tracking algorithm ===")
    
    u1, u2 = bisect_to_edge(pinteg, mapper; abstol=eps1)
    edgestate = (u1 + u2)/2
    
    if output_level>=1
        track1 = [u1]
        track2 = [u2]
        edge = [edgestate]
    end
    # Initialize Lyapunov "test" trajectories
    pinteg_test = deepcopy(pinteg)
    inittest = inittest_default(size(get_state(pinteg, 1))[1])
    λ1 = zero(d0)
    λ2 = zero(d0)
    titer = zero(d0)
    # End Lyapunov init

    verbose && println("... Iteration 1: Edge at $(edgestate)")

    correction = converge + 1
    counter = 1
    while (correction > converge) & (maxit > counter)
        # Initialize/rescale Lyapunov calculation
        if calc_lyapunovs
            if counter == lyapunov_spinup
                println("Initializing Lyapunov")
                reinit!(pinteg_test, [inittest(u1, d0), inittest(u2, d0)])
                d_test1 = deepcopy(d0)
                d_test2 = deepcopy(d0)
            elseif counter > lyapunov_spinup
                # get distance between test and shadowing trajectories before rescaling the latter
                d_test1 = norm(get_state(pinteg_test, 1) - get_state(pinteg, 1))
                d_test2 = norm(get_state(pinteg_test, 2) - get_state(pinteg, 2))
                println("Re-initializing Lyapunov, distances:", d_test1, " ", d_test2)
                # re-adjust (but not re-scale) trajectory
                #reinit!(pinteg_test, [inittest(u1, d_test1), inittest(u2, d_test2)]) #OLD
                reinit!(pinteg_test, [
                    u1 + (pinteg_test.u[1] - pinteg.u[1])/(d_test1/d0),
                    u2 + (pinteg_test.u[2] - pinteg.u[2])/(d_test2/d0)
                ])
            end
        end
        # Re-initalize shawding trajectories
        reinit!(pinteg, [u1,u2])
        state = edgestate
        distance = norm(get_state(pinteg, 1)-get_state(pinteg, 2))
        while distance < eps2
            step!(pinteg, dt)
            if calc_lyapunovs && counter>=lyapunov_spinup
                step!(pinteg_test, dt)
                d_test1 = norm(get_state(pinteg_test, 1) - get_state(pinteg, 1))
                d_test2 = norm(get_state(pinteg_test, 2) - get_state(pinteg, 2))
                #println(titer, " ", d_test1, " ", d_test2)
                if (d_test1 < lower_threshold) || (d_test1 > upper_threshold)
                    # "threshold" rescale
                    a1 = d_test1/d0
                    #println("Rescaling state 1: a1 = ", a1)
                    λ1 += log(a1)
                    #pinteg_test.u[1] = inittest(get_state(pinteg, 1), d0)
                    pinteg_test.u[1] = pinteg.u[1] + (pinteg_test.u[1] - pinteg.u[1])/a1
                    u_modified!(pinteg_test, true)
                    d_test1 = norm(get_state(pinteg_test, 1) - get_state(pinteg, 1))
                    #println(d_test1)
                end
                if (d_test2 < lower_threshold) || (d_test2 > upper_threshold)
                    a2 = d_test2/d0
                    #println("Rescaling state 2: a2 = ", a2)
                    λ2 += log(a2)
                    #pinteg_test.u[2] = inittest(get_state(pinteg, 2), d0)
                    pinteg_test.u[2] = pinteg.u[2] + (pinteg_test.u[2] - pinteg.u[2])/a2
                    u_modified!(pinteg_test, true)
                    d_test2 = norm(get_state(pinteg_test, 2) - get_state(pinteg, 2))
                    #println(d_test2)
                end
                titer += 1
            end
            state1 = get_state(pinteg, 1)
            state2 = get_state(pinteg, 2)
            distance = norm(state1-state2)
            if output_level>=2
                push!(track1, state1)
                push!(track2, state2)
                push!(edge, (state1+state2)/2)
            end
        end
        u1, u2 = bisect_to_edge(pinteg, mapper; abstol=eps1)
        edgestate = (u1 + u2)/2
        correction = norm(edgestate - state)
        counter += 1

        if output_level==1
            push!(track1, u1)
            push!(track2, u2)
            push!(edge, edgestate)
        end
        
        verbose && println("... Iteration $(counter): Edge at $(edgestate)")
        (counter == maxit) && @warn("Reached maximum number of $(maxit) iterations; did not converge.")
    end

    (counter < maxit) && println("Edge-tracking converged after $(counter) iterations.")

    if calc_lyapunovs
        # Do finale rescale
        d_test1 = norm(get_state(pinteg_test, 1) - get_state(pinteg, 1))
        d_test2 = norm(get_state(pinteg_test, 2) - get_state(pinteg, 2))
        # store growth rates
        a1 = d_test1/d0
        a2 = d_test2/d0
        λ1 += log(a1)
        λ2 += log(a2)
        println("Lyapunov exponents:\n", λ1/(titer*dt), " ", λ2/(titer*dt), "(iters = ", titer, ")")
    end

    if output_level>=1
        return Dataset(edge), Dataset(track1), Dataset(track2)
    else
        return edgestate
    end
end;

"""
    bisect_to_edge(pinteg, mapper; abstol=1e-9)

Find the basin boundary between two states `u1, u2 = get_states(pinteg)` by bisecting along a
straight line in phase space. The states `u1` and `u2` must belong to different basins.
Returns two new states located on either side of the basin boundary at a maximum distance of
`abstol` between each other.

`pinteg` is a `parallel_integrator` with two states. The `mapper` must be an `AttractorMapper`
of subtype `AttractorsViaProximity` or `AttractorsViaRecurrences`.

Note: If the straight line between `u1` and `u2` intersects the basin boundary multiple
times, the method will find one of these intersection points. If more than two attractors
exist, one of the two returned states may belong to a different basin than the initial
conditions `u1` and `u2`. A warning is raised if the bisection involves a third basin.

# Keyword arguments
* `abstol = 1e-9`: The maximum (Euclidean) distance between the two returned states.
"""
function bisect_to_edge(pinteg, mapper::AttractorMapper; abstol=1e-9)
    u1, u2 = get_states(pinteg)
    idx1, idx2 = mapper(u1), mapper(u2)
    
    if idx1 == idx2
        error("Both initial conditions belong to the same basin of attraction: ", idx1)
    end
    
    distance = norm(u1-u2)
    while distance > abstol
        u_new = (u1 + u2)/2
        idx_new = mapper(u_new)
        # todo: what if u_new lies on boundary
        
        if idx_new == idx1
            u1 = u_new
        else 
            u2 = u_new
            if idx_new != idx2
                println("Warning: New bisection point belongs to a third basin of attraction or diverges: ", idx_new) # @warn
            end
        end    
        distance = norm(u1-u2)
    end
    [u1, u2]
end;
