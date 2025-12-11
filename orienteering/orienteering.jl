using JuMP
using Gurobi

# Orienteering instance
struct Instance
    n::Int                     # Number of locations
    Tmax::Float64              # Maximum available time
    dist::Matrix{Float64}      # Travel time/cost matrix d[i,j]
    scores::Vector{Float64}    # Score s[i] for each location i

end

# Orienteering solution
struct Solution
    route::Vector{Int}      # Sequence of locations visited (including depot)
    score::Float64          # Total score reached
    distance::Float64       # Total distance traveled
end

# Reads in an instance from a given file path
function read_instance(path::String)
    # Read the file header
    lines = readlines(path)
    n = parse(Int, lines[1])
    Tmax = parse(Float64, lines[2])
    # Read the distance matrix
    dist = Array{Float64}(undef, n, n)
    for i in 1:n
        row = split(lines[2 + i])
        dist[i, :] = parse.(Float64, row)
    end
    # Read the scores
    scores = Vector{Float64}(undef, n)
    for i in 1:n
        scores[i] = parse(Float64, lines[2 + n + i])
    end
    # Create the instance
    return Instance(n, Tmax, dist, scores)
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
        # Check if the solution is integer feasible
        if !all(val -> isapprox(val, round(val); atol=1e-6), x_val)
            return
        end
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

# Main function
function main()
    # Check the command line arguments
    if length(ARGS) < 1
        println("Usage: julia orienteering.jl <instance-file>")
        return
    end

    # Read the instance
    filepath = ARGS[1]
    inst = read_instance(filepath)
    
    # Print information about the instance
    println("Number of locations: ", inst.n)
    println("Time limit: ", inst.Tmax)
    
    # Solve the instance
    solve_orienteering(inst)

end

# Start the main function
main()

