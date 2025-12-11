using DataStructures
using Random, Dates

# Instance of P||Cmax
struct Instance
    m::Int                  # number of machines
    n::Int                  # number of jobs
    p::Vector{Int}          # processing times p_i
end

# Solution of P||Cmax
struct Solution
    assignment::Vector{Int} # machine assignment for each job (length n)
    makespan::Int           # resulting makespan
end

# Reads an instance from a file
function read_instance(path::String)
    lines = readlines(path)
    # Parse header
    m = parse(Int, lines[1])
    n = parse(Int, lines[2])
    # Parse processing times
    p = [parse(Int, lines[i]) for i in 3:(2+n)]
    # Create the instance
    return Instance(m, n, p)
end

# Decodes an indirect solution representation consisting of a list of job indices by applying list scheduling, i.e., scheduling the next job on the machine with the least load.
function decode(inst::Instance, perm::Vector{Int})
    # Create a priority queue to efficiently obtain the minimum load machine
    heap = BinaryMinHeap{Tuple{Int,Int}}()
    for machine in 1:inst.m
        push!(heap, (0, machine)) # At the beginning the load is 0
    end
    # Create an empty assignment of jobs to machines
    assignment = Vector{Int}(undef, inst.n)
    # Makespan counting variable
    makespan = 0
    # Process all jobs according the given permutation
    for j in perm
        # Select the minimum load machine
        (minload, machine) = pop!(heap)
        # Update load and job assignment
        newload = minload + inst.p[j]
        makespan = max(makespan, newload)
        assignment[j] = machine
        # Update priority queue
        push!(heap, (newload, machine))
    end
    # Return the created solution
    return Solution(assignment, makespan)
end

function simulated_annealing(inst::Instance; 
                             time_limit::Float64=10.0)

    # Initialisierung
    T0 = 100.0
    alpha = 0.95

    # Startzeit
    start_time = now()

    # Initial solution: Sort the jobs in decreasing order by their processing times and apply list scheduling, yielding a 4/3-approximation algorithm
    perm = Vector{Int}(1:inst.n)
    perm = sortperm(inst.p; rev=true)
    current_sol = decode(inst, perm)
    
    # Preparation for the main loop
    best_sol = current_sol
    T = T0

    # Main loop
    while (now() - start_time).value / 1e3 < time_limit * 1000
        # Swap two random jobs
        i, j = rand(1:inst.n, 2)
        perm[i], perm[j] = perm[j], perm[i]

        new_sol = decode(inst, perm)

        delta = new_sol.makespan - current_sol.makespan

        # Accept deteriorating solutions only with a given probability
        if delta < 0 || rand() < exp(-delta / T)
            current_sol = new_sol
            if new_sol.makespan < best_sol.makespan
                best_sol = new_sol
            end
        # Undo swap
        else
            perm[i], perm[j] = perm[j], perm[i]
        end

        # Cooling
        T *= alpha
    end

    # Return the best solution found
    return best_sol
end



function main()
    # Check the command line arguments
    if length(ARGS) < 1
        println("Usage: julia scheduling.jl <instance-file>")
        return
    end

    # Read the instance
    filepath = ARGS[1]
    inst = read_instance(filepath)
    
    # Print information about the instance
    println("Number of machines: ", inst.m)
    println("Number of jobs: ", inst.n)
    println("Processing times: ", inst.p)
    
    # Apply decode
    perm = Vector{Int}(1:inst.n)
    sol = decode(inst, perm)
    println(sol)
    
    # Apply simulated annealing
    sol = simulated_annealing(inst; time_limit=10.0)
    
end

# Start the main function
main()
