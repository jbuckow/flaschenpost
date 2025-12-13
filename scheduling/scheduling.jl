using Dates
using DataStructures
using Random

"""
    Instance

Represents an instance of the parallel machine scheduling problem (P||Cmax).

# Fields
- `m::Int`        : Number of machines
- `n::Int`        : Number of jobs
- `p::Vector{Int}`: Processing times `p_j` for each job
"""
struct Instance
    m::Int
    n::Int
    p::Vector{Int}
end

"""
    Solution

Represents a solution of the parallel machine scheduling problem (P||Cmax).

# Fields
- `makespan::Int`: The resulting makespan
- `assignment::Vector{Int}`: Assignment of jobs to machines
"""
struct Solution
    makespan::Int
    assignment::Vector{Int}
end

"""
    validate_solution(inst::Instance, sol::Solution) -> Bool

Validates if a given solution is feasible for the problem instance `inst`.

# Checks
- Each job is assigned to exactly one machine.
- Machine indices are within the valid range (1..m).
- The makespan equals the maximum load across all machines.

# Arguments
- `inst::Instance`: Problem instance.
- `sol::Solution`: Solution to be verified.

# Returns
- `Bool`: `true` if the solution is feasible, otherwise `false`.
"""
function validate_solution(inst::Instance, sol::Solution)::Bool
    # Check assignment length
    if length(sol.assignment) != inst.n
    	println("Error: assignment length = $(length(sol.assignment)), expected = $(inst.n)")
        return false
    end
    # Check machine indices
    if any(machine -> machine < 1 || machine > inst.m, sol.assignment)
        println("Error: invalid machine index found in assignment (valid range: 1..$(inst.m))")
        return false
    end
    # Compute load per machine
    loads = zeros(Int, inst.m)
    for (job, machine) in enumerate(sol.assignment)
        loads[machine] += inst.p[job]
    end
    # Check makespan
    computed_makespan = maximum(loads)
    if computed_makespan != sol.makespan
        println("Error: Makespan mismatch. Calculated = $computed_makespan; Solution = $(sol.makespan)")
    	return false
    end
    # Solution is feasible
    return true
end

"""
    read_instance(path::String) -> Instance

Reads an instance of the parallel machine scheduling problem (P||Cmax) from a text file.

# Arguments
- `path::String`: Path to the input file.  
  The file format is expected as:
  1. First line: number of machines (m)
  2. Second line: number of jobs (n)
  3. Next n lines: processing times `p_i` for each job

# Returns
- `Instance`: The instances red in as struct
"""

function read_instance(path::String)
    lines = readlines(path)
    # Parse header
    m = parse(Int, lines[1])
    n = parse(Int, lines[2])
    # Parse processing times
    p = [parse(Int, lines[j]) for j in 3:(2+n)]
    # Create the instance
    return Instance(m, n, p)
end

"""
    decode(inst::Instance, perm::Vector{Int}) -> Solution

Decodes an indirect solution representation for the parallel machine scheduling problem (P||Cmax).

The function applies **list scheduling**: each job in the given permutation is scheduled on a machine with the current minimum load. A priority queue is used to efficiently select the least loaded machine at each step.

# Arguments
- `inst::Instance`: Input instance
- `perm::Vector{Int}`: A permutation of job indices representing the order in which jobs are scheduled.

# Returns
- `Solution`: Decoded direct solution representation as a struct.
"""
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

"""
    Timer

A mutable struct representing a timer for simulated annealing or other iterative algorithms.

# Fields
- `time_limit::Float64` Time limit in seconds for the algorithm run.

- `start_time::Union{DateTime, Nothing}`  
  The recorded start time of the timer.

- `iterations::Int`  
  Counter for the number of iterations performed during the run.
"""
mutable struct Timer
    time_limit::Float64
    start_time::Union{DateTime, Nothing}
    iterations::Int
end

"""
    Timer(time_limit::Float64) -> Timer

Creates and initializes a `Timer` object with a given time limit.

# Arguments
- `time_limit::Float64`: The maximum allowed runtime in seconds.

# Returns
- `Timer`: A timer with the specified time limit.
"""
function Timer(time_limit::Float64)
    return Timer(time_limit, now(), 0)
end

"""
    expired(t::Timer) -> Bool

Checks whether the time limit has been reached.

# Arguments
- `t::Timer`: The timer.

# Returns
- `Bool`: `true` if the time limit is exceeded, otherwise `false`.
"""
function expired(t::Timer)::Bool
    return elapsed(t) >= t.time_limit
end

"""
    elapsed(t::Timer) -> Float64

Calculates the elapsed time since the start of the timer.

# Arguments
- `t::Timer`: The timer

# Returns
- `Float64`: The elapsed time in seconds since the timer was created.
"""
function elapsed(t::Timer)::Float64
    elapsed = (now() - t.start_time).value / 1000.0
    return elapsed
end

"""
    tick!(t::Timer)

Increments the iteration counter of the timer.

# Arguments
- `t::Timer`: The timer

# Returns
- `Nothing`: Just updates the timer.
"""
function tick!(t::Timer)
    t.iterations += 1
end

"""
    simulated_annealing(inst::Instance; 
                        time_limit::Float64=10.0,
                        start_temperature_factor::Float64=0.0001,
                        end_temperature_factor::Float64=0.000001) -> Solution

Performs the **simulated annealing metaheuristic** to solve the parallel machine scheduling problem (P||Cmax).

The algorithm works with an indirect solution representation in the shape of a permutation of all jobs. These indirect solutions are decoded into actual schedules by applying the list scheduling procedure, where the jobs are put one after the other to a machine with the currently least load. 

Since the procedure starts with an initial solution obtained by sorting jobs in decreasing order by their processing times, it yields a 4/3-approximation guarantee. Afterwards, in the main loop, neighbor solutions are created by swapping jobs in the permutation and applying a randomized simulated annealing acceptance criterion to escape local minima. A dynamic version of geometric cooling is applied to reach a given target temperature (scaled down by the objective value). After a given number of non-improving iterations, a random partially shuffling is applied to diversify. 

# Arguments
- `inst::Instance`: Problem instance.
- `time_limit::Float64=10.0`: Solving time limit in seconds.
- `seed::Int=425234`: The pseudo random number generator seed.
- `start_temperature_factor::Float64=0.0001`: Factor to initialize the starting temperature relative to the makespan.
- `end_temperature_factor::Float64=0.000001`: Factor to determine the target final temperature relative to the makespan.

# Returns
- `Solution`: The best solution found within the solving time limit.
"""
function simulated_annealing(inst::Instance; 
                             time_limit::Float64=10.0,
                             seed::Int32=425234,
                             start_temperature_factor::Float64=0.0001,
                             end_temperature_factor::Float64=0.000001)

    # Start the timer
    timer = Timer(time_limit)
    
    # Set the seed
    Random.seed!(seed)

    # Initial solution: Sort the jobs in decreasing order by their processing times and apply list scheduling, yielding a 4/3-approximation algorithm
    perm = Vector{Int}(1:inst.n)
    perm = sortperm(inst.p; rev=true)
    current_sol = decode(inst, perm)
    
    # Preparation for the main loop
    best_sol = current_sol
    temperature = start_temperature_factor * best_sol.makespan
    alpha = 0.9999

    # Main loop
    non_improving = 0
    while !expired(timer)
        # Shuffle part of the solution after many non-improving solutions
        if non_improving >= 1000000
            len = div(inst.n, 2)
            start = rand(1:(inst.n-len+1))
            stop = start + len - 1
            shuffle!(view(perm, start:stop))
            non_improving = 0
        end
        # Swap two random jobs
        i, j = rand(1:inst.n, 2)
        perm[i], perm[j] = perm[j], perm[i]

        # Create new solution
        new_sol = decode(inst, perm)
        
        # Accept deteriorating solutions only with a given probability
        delta = new_sol.makespan - current_sol.makespan
        if delta < 0 || rand() < exp(-delta / temperature)
            current_sol = new_sol
            if new_sol.makespan < best_sol.makespan
                best_sol = new_sol
                non_improving = 0
            else
                non_improving += 1
            end
        # Undo swap
        else
            perm[i], perm[j] = perm[j], perm[i]
        end

        # Adapt the cooling procedure dynamically such that the target temperature is reached in the end (independent of the time limit)
        if (timer.iterations > 5 && timer.iterations % 1000 == 0)
            elapsed_seconds = elapsed(timer)
            iters_per_second = timer.iterations / elapsed_seconds
            time_left_seconds = timer.time_limit - elapsed_seconds
            iters_remaining = time_left_seconds * iters_per_second
            target_temperature = best_sol.makespan * end_temperature_factor
            alpha = min(1.0, (target_temperature / temperature) ^ (1 / iters_remaining))
        end

        # Cooling
        temperature *= alpha
        
        # Print information about the solving process
        if (timer.iterations % 10000 == 0)
            println("iteration=", timer.iterations, " best=", best_sol.makespan, " curr=", current_sol.makespan, " new=", new_sol.makespan)
        end
        
        # Increment the iteration counter
        tick!(timer)
    end

    # Return the best solution found
    return best_sol
end

"""
    main()

Entry point of the scheduling program. Reads an instance from a file, applies the
simulated annealing algorithm, and prints the resulting solution.

# Usage
Run from the command line: julia scheduling.jl <instance-file> <time-limit-seconds>

# Arguments
- `ARGS[1]`: Path to the instance file.
- `ARGS[2]`: Solving time limit in seconds.
- `ARGS[3]`: Seed for the random number generator.

# Returns
- `Nothing`: Prints the solution and its makespan to the console.
"""
function main()
    # Check the command line arguments
    if length(ARGS) < 3
        println("Usage: julia scheduling.jl <instance-file> <time-limit-seconds> <seed>")
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
    sol = simulated_annealing(inst; time_limit=parse(Float64, ARGS[2]), seed=parse(Int32, ARGS[3]))
    
    # Check if the given solution is feasible
    if !validate_solution(inst, sol)
        println("Error: solution is infeasible.")
        return
    end
    
    # Print the solution, including its makespan
    println("Assignment of jobs to machines:")
    for (job, machine) in enumerate(sol.assignment)
        println("Job $job -> Machine $machine")
    end
    println("Calculated makespan: ", sol.makespan)
end

# Start the main function
main()
