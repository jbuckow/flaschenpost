using DataStructures

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



jobs = [12, 4, 18, 5, 29]
m = 3

ub = calculate_ub(m, jobs)
println("Upper Bound: ", ub)
