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
    makespan::Int           # resulting makespan
    assignment::Vector{Int} # assignment of jobs to machines (length n)
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
    return Solution(makespan, assignment)
end

# Timer struct for simulated annealing
mutable struct Timer
    time_limit::Float64   # Time limit in seconds
    start_time::Union{DateTime, Nothing} # The start time
    iterations::Int # The number of iterations performed
end

# Constructor for the timer
function Timer(time_limit::Float64)
    return Timer(time_limit, now(), 0)
end

# Checks if the time limit is reached
function expired(t::Timer)::Bool
    if t.start_time === nothing
        error("Timer not started!")
    end
    return elapsed(t) >= t.time_limit
end

# Calculated the remaining time
function elapsed(t::Timer)::Float64
    if t.start_time === nothing
        error("Timer not started!")
    end
    elapsed = (now() - t.start_time).value / 1000.0
    return elapsed
end

# Increment the iteration counter
function tick!(t::Timer)
    t.iterations += 1
end

# Simulated annealing main function
function simulated_annealing(inst::Instance; 
                             time_limit::Float64=10.0)

    # Initialisierung
    T0 = 100.0
    alpha = 0.999

    # Startzeit
    timer = Timer(time_limit)

    # Initial solution: Sort the jobs in decreasing order by their processing times and apply list scheduling, yielding a 4/3-approximation algorithm
    perm = Vector{Int}(1:inst.n)
    perm = sortperm(inst.p; rev=true)
    current_sol = decode(inst, perm)
    
    # Preparation for the main loop
    best_sol = current_sol
    T = T0

    # Main loop
    while !expired(timer)
        # Swap two random jobs
        i, j = rand(1:inst.n, 2)
        perm[i], perm[j] = perm[j], perm[i]

        # Create new solution
        new_sol = decode(inst, perm)
        
        # Accept deteriorating solutions only with a given probability
        delta = new_sol.makespan - current_sol.makespan
        rnd = rand()
        ex = exp(-delta / T)
        if delta < 0 || rnd < ex
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
        
        if (timer.iterations % 1000 == 0)
            println("iteration=", timer.iterations, " best=", best_sol.makespan, " curr=", current_sol.makespan, " new=", new_sol.makespan)
        end
        
        # Increment the iteration counter
        tick!(timer)
    end

    # Return the best solution found
    return best_sol
end

# Main function
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
    
    # Apply simulated annealing and print the solution
    sol = simulated_annealing(inst; time_limit=10.0)
    println("Assignment of jobs to machines:")
    for (job, machine) in enumerate(sol.assignment)
        println("Job $job -> Machine $machine")
    end
    println("Calculated makespan = ", sol.makespan)
end

# Start the main function
main()
