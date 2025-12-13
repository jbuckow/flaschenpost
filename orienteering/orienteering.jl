using JuMP
using Gurobi

"""
    Instance

Data structure representing an Orienteering Problem instance.

# Fields
- `n::Int`: Number of locations in the instance.
- `Tmax::Float64`: Maximum available travel time.
- `dist::Matrix{Float64}`: Distance or travel time matrix `d[i,j]` between locations.
- `scores::Vector{Float64}`: Score `s[i]` associated with each location `i`.

# Notes
- The depot is typically included as one of the locations (index 1).
"""
struct Instance
    n::Int
    Tmax::Float64
    dist::Matrix{Float64}
    scores::Vector{Float64}
end

"""
    Solution

Represents a solution for the Orienteering Problem.

# Fields
- `route::Vector{Int}`: The ordered sequence of visited locations, including the depot (location 1) at the start and end of the route. Example: `[1, 3, 5, 1]` means the route starts at depot 1, visits locations 3 and 5, and returns to depot 1.
- `total_score::Float64`: The sum of the scores of all visited locations.
"""
struct Solution
    route::Vector{Int}
    total_score::Float64
end

"""
    validate_solution(inst::Instance, sol::Solution) -> Bool

Validates whether a given solution is feasible for a given Orienteering Problem instance.

# Checks
- The tour starts and ends at the depot (location 1).
- Each location is visited at most once.
- The total travel time of the route does not exceed the time limit.
- The total score of visited locations matches the score in the solution.

# Returns
- `true` if the solution is feasible, otherwise `false`.
"""
function validate_solution(inst::Instance, sol::Solution)::Bool
    # Check depot at start and end
    if sol.route[1] != 1 || sol.route[end] != 1
        println("Error: Route must start and end at depot location (1).")
        return false
    end
    # Check uniqueness of non-depot visits
    if length(sol.route) != length(Set(sol.route)) + 1
        println("Error: A location is visited more than once.")
        return false
    end
    # Check the tour length
    tour_length = 0.0
    for i in 1:(length(sol.route)-1)
        tour_length += inst.dist[sol.route[i], sol.route[i+1]]
    end
    if tour_length > inst.Tmax
        println("Error: Tour length = $tour_length exceeds Tmax = $(inst.Tmax)")
        return false
    end
    # Check the score
    visited_locations = Set(sol.route)
    computed_score = sum(inst.scores[j] for j in visited_locations)
    if computed_score != sol.total_score
        println("Error: Score mismatch. Calculated = $computed_score; Solution = $(sol.total_score)")
        return false
    end
    # Otherwise the solution is feasible
    return true
end

"""
    read_instance(path::String) -> Instance

Reads a problem instance from a text file.

# Arguments
- `path::String`: Path to the instance file.  
  Expected format:
  1. First line: number of locations `n` (Int)
  2. Second line: time limit `Tmax` (Float64)
  3. Next `n` lines: distance matrix (Float64 values separated by spaces)
  4. Next `n` lines: scores for each node (Float64)

# Returns
- `Instance`: An orienteering instance object.
"""
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

"""
    extract_tour(x_val::AbstractArray{<:Real,2}, L::AbstractVector{<:Integer}, depot::Int) -> Vector{Int}

Constructs a route starting from a given depot location based on the values of two-index variables.

# Arguments
- `x_val::AbstractArray{<:Real,2}`: A 2D array representing the values of the two-index variables
- `L::AbstractVector{<:Integer}`: Set of locations.
- `depot::Int`: Starting location (depot).

# Returns
- `Vector{Int}`: The sequence of locations visited in the tour.
"""
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

"""
    solve_orienteering(inst::Instance) -> Solution

Builds and solves the Orienteering Problem as a Mixed Integer Program (MIP).

# Arguments
- `inst::Instance`: The problem instance

# Returns
- `Bool`: `true` if the instance was solved successfully and `false` otherwise.
"""
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
    # Extract the tour
        xsol = value.(x)
        tour = extract_tour(xsol, L, 1)
        push!(tour, 1)
        # Create and return solution
        sol = Solution(tour, objective_value(model))
        return sol
    # If no solution was found
    else
        println("Solver terminated with status: ", status)
        return nothing
    end
end

"""
    main()

Entry point of the program. Handles command-line arguments, reads the problem instance, prints basic information, and calls the solver.

# Usage
Run the program from the command line: julia orienteering.jl <instance-file>

# Returns
- Nothing. Prints results or usage instructions to the console.
"""
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
    sol = solve_orienteering(inst)
    
    # Validate and print the solution
    if sol != nothing
        if !validate_solution(inst, sol)
            println("Error: solution is infeasible.")
            return
        end
        println("Total score: ", sol.total_score)
        println("Tour: ", sol.route)
    end
end

# Start the main function
main()

