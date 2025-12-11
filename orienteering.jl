using JuMP
using Gurobi

# Orienteering instance
struct Instance
    n::Int                     # Number of locations
    scores::Vector{Float64}    # Score s[i] for each location i
    dist::Matrix{Float64}      # Travel time/cost matrix d[i,j]
    Tmax::Float64              # Maximum available time
end

# Orienteering solution
struct Solution
    route::Vector{Int}      # Sequence of locations visited (including depot)
    score::Float64          # Total score reached
    distance::Float64       # Total distance traveled
end

# Create a small example instance
function example_instance()
    n = 5
    scores = [0.0, 10.0, 15.0, 20.0, 25.0]   # Location 1 = depot with score 0
    dist = [
        0.0   10.0  20.0  15.0  30.0;
        10.0  0.0   25.0  35.0  20.0;
        20.0  25.0  0.0   30.0  15.0;
        15.0  35.0  30.0  0.0   10.0;
        30.0  20.0  15.0  10.0  0.0
    ]
    Tmax = 60.0
    return Instance(n, scores, dist, Tmax)
end

# Helper function to extract a route from a starting point
function extract_tour(x_val::AbstractArray{<:Real,2}, L::AbstractVector{<:Integer}, depot::Int)
    current = depot
    route = [depot]
    while true
        succ = findfirst(j -> x_val[current, j] > 0.5, L)
        if succ == depot
            break
        end
        push!(route, succ)
        current = succ
    end
    return route
end

# Build and solve the problem as MIP
function solve_orienteering(inst::Instance)
    # Define set of all locations
    L = 1:inst.n

    # Create a Gurobi model
    model = Model(Gurobi.Optimizer)

    # Binary arc variables: x[i,j] = 1 if arc (i,j) is used
    @variable(model, x[i in L, j in L], Bin)

    # Binary location variables: y[i] = 1 if location i is visited
    @variable(model, y[i in L], Bin)

    # Objective: maximize collected scores
    @objective(model, Max, sum(inst.scores[i] * y[i] for i in L))

    # Linking x- and y-variables and degree balancing
    @constraint(model, [i in L], sum(x[i,j] for j in L if j != i) == y[i])
    @constraint(model, [i in L], sum(x[j,i] for j in L if j != i) == y[i])

    # Ensuring that the depot node is always visited
    @constraint(model, y[1] == 1)

    # Constraint ensuring the time limit is not exceeded
    @constraint(model, sum(inst.dist[i,j] * x[i,j] for i in L, j in L) <= inst.Tmax)
    
    # Lazy callbacks to eliminate subtours
    function lazy_cb(cb_data)
        # Extract current solution
        x_val = callback_value.(cb_data, x)
        y_val = callback_value.(cb_data, y)
        # Determine the set S of all locations in the same route as the depot
        S = extract_tour(x_val, L, 1)
        # If only a single route exists (i.e., S covers all visited locations)
        if (length(S) >= count(i -> y_val[i] > 0.5, L))
            return
        end
        # Add cuts if multiple routes exist
        S_comp = setdiff(L, S) # All remaining locations not in S
        for k in S_comp
            # Add connectivity constraints between S and S_comp (that are only enabled depending on y_k)
            con_out = @build_constraint(y[k] <= sum(x[i,j] for i in S, j in S_comp))
            con_in  = @build_constraint(y[k] <= sum(x[j,i] for j in S_comp, i in S))
            MOI.submit(model, MOI.LazyConstraint(cb_data), con_out)
            MOI.submit(model, MOI.LazyConstraint(cb_data), con_in)
        end
    end

    # Activate lazy callbacks
    MOI.set(model, MOI.LazyConstraintCallback(), lazy_cb)

    # Branching strategy focusing on the inclusion or exclusion of locations by setting branching priorities
    for i in L
        MOI.set(model, Gurobi.VariableAttribute("BranchPriority"), y[i], 10)
    end

    # Optimize
    optimize!(model)

    # Collect results
    status = termination_status(model)
    if status == MOI.OPTIMAL
        println("Optimal score: ", objective_value(model))
        println("Visited nodes: ", [i for i in L if value(y[i]) > 0.5])
        xsol = value.(x)
        tour = extract_tour(xsol, L, 1)
        print(tour)
        println(xsol)
        return model
    else
        println("Solver terminated with status: ", status)
        return nothing
    end
end

# Solve the example instance
inst = example_instance()
solve_orienteering(inst)

