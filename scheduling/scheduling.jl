using DataStructures

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

# Decodes an indirect solution representation consisting only of a list of processing times by applying list scheduling, i.e., scheduling the next job to the machine with the least load.
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
end

# Start the main function
main()
