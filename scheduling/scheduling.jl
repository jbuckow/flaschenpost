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
    # Parse first two lines
    m = parse(Int, lines[1])
    n = parse(Int, lines[2])
    # Parse processing times
    p = [parse(Int, lines[i]) for i in 3:(2+n)]
    # Create the instance
    return Instance(m, n, p)
end

function calculate_ub(m::Int, pjs::Vector{Int})
    # Sortiere Jobs absteigend (Largest Processing Time zuerst)
    #sort!(pjs, rev=true)

    # Min-Heap mit m Nullen (Startlasten der Maschinen)
    heap = BinaryMinHeap{Int}()
    for _ in 1:m
        push!(heap, 0)
    end

    # Laufendes Maximum der Maschinenlasten (Makespan)
    max_load = 0

    # Verteile Jobs: jeweils auf die Maschine mit kleinster Last
    for pj in pjs
        minload = pop!(heap)          # kleinste aktuelle Last
        newload = minload + pj        # aktualisierte Last dieser Maschine
        max_load = max(max_load, newload)
        push!(heap, newload)          # zurück in den Min-Heap
    end

    return max_load
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
end

# Start the main function
main()
