using Distributed
addprocs(8)

@everywhere begin
    using Combinatorics
    using LinearAlgebra
    # using SparseArrays
    using Printf
    using BenchmarkTools
    using .Threads
    # import Pkg; Pkg.add("ThreadsX")
    using ThreadsX
    using MKL
end

# import Distributed # for addprocs() and nprocs()
# import CpuId
# Distributed.addprocs(4-Distributed.nprocs())


# function symm_matrix_quad_form(wedge_basis_idx::Vector{Vector{Int}}, q::Vector{Int})
#     elem_idx = collect(combinations(1:4, 2))
#     rows = Vector{Int}(undef, 3)
#     cols = Vector{Int}(undef, 3)
#     for i in 1:3
#         rows[i] = findfirst(x -> x == q[elem_idx[i]], wedge_basis_idx)
#         cols[i] = findfirst(x -> x == q[elem_idx[7-i]], wedge_basis_idx)
#     end
#     vals = [1,-1,1]
#     vals = [vals; vals]
#     rows_tmp = rows
#     rows = [rows; cols]
#     cols = [cols; rows_tmp]
#     return sparse(rows, cols, vals, length(wedge_basis_idx), length(wedge_basis_idx))
# end

function symm_matrix_quad_form(wedge_basis_idx::Vector{Vector{Int}}, q::Vector{Int})
    elem_idx = collect(combinations(1:4, 2))
    k = length(wedge_basis_idx)
    mat = zeros(Int, k, k)
    vals = [1, -1, 1]
    for i in 1:3
        row = findfirst(x -> x == q[elem_idx[i]], wedge_basis_idx)
        col = findfirst(x -> x == q[elem_idx[7 - i]], wedge_basis_idx)
        mat[row, col] = vals[i]
        mat[col, row] = vals[i]
    end
    return mat
end


@everywhere function signature_matrix(mat::Matrix{Int})
    evals = eigvals(Symmetric(mat))
    p = count(x -> x > 1e-8, evals)
    q = count(x -> x < -1e-8, evals)
    r = length(evals) - (p + q)
    return (p, q, r)
end

function read_signatures_file(filename::String)
    signatures = Vector{Tuple{Int, Int, Int}}()
    quad_forms = Vector{Vector{Vector{Int}}}()

    for line in eachline(filename)
        parts = split(line, ',')
        sig = Tuple(parse.(Int, parts[1:3]))

        # The rest are 4-forms encoded as stringified arrays: e.g., "[1 2 3 4]"
        forms = map(parts[4:end]) do s
            parse.(Int, split(strip(s, ['[', ']']), ' '))
        end

        push!(signatures, sig)
        push!(quad_forms, forms)
    end

    return signatures, quad_forms
end

function print_signatures(filename::String, signatures::Vector{Tuple{Int, Int, Int}}, quad_forms_basis_idx::Vector{Vector{Vector{Int}}}, mode::String = "w")
     # Output file
    open(filename, mode) do io
        for (sgn,combo) in zip(signatures, quad_forms_basis_idx)
            @printf(io, "%d,%d,%d", sgn...)
            for q in combo
                @printf(io, ",[%s]", join(q, " "))
            end
            write(io, "\n")
        end
    end
end

@everywhere function hash_sgn(sgn::NTuple{3, Int})
    # Convert the signature tuple to a 32-bit integer hash
    # This is a simple hash function that combines the three integers
    # using bitwise operations
    # Ensure the signature is a tuple of three integers
    # assumes sgn[i] < 256 i.e. n < 23 as 23 choose 2 = 253
    return sgn[1] + (sgn[2] << 8) + (sgn[3] << 16)
    # or just use Set{Tuple{Int,Int,Int}} or hash(sgn)
end

# Function to compute unique signatures of quadrilinear forms
function unique_signatures_quad_forms(n::Int, min_comb_size::Int=1, max_comb_size::Int=4, unique_signatures = Vector{Tuple{Int, Int, Int}}(), unique_quad_forms = Vector{Vector{Vector{Int}}}())
    wedge_basis_idx = collect(combinations(1:n, 2))
    quad_forms_basis_idx = collect(combinations(1:n, 4))
    # k = length(wedge_basis_idx)
    m = length(quad_forms_basis_idx)

    # Precompute symmetric matrices
    mats = [symm_matrix_quad_form(wedge_basis_idx, q) for q in quad_forms_basis_idx]

    # Signature cache
    seen_signatures = Set{UInt64}()
    for sgn in unique_signatures
        push!(seen_signatures, hash_sgn(sgn))
    end

    skipped_combos = Threads.Atomic{Int}(0)

    lk = ReentrantLock()
    # Find unique signatures
    for i in min_comb_size:max_comb_size
        @info "Processing combinations of size $i"
        combos = collect(combinations(1:m, i))

        # Multithreaded processing of combinations
        ThreadsX.foreach(combos) do combo
            # Early filter: check support spread
            all_indices = vcat(quad_forms_basis_idx[combo]...)
            if length(unique(all_indices)) < i + 3  # Heuristic: must span enough indices
                Threads.atomic_add!(skipped_combos, 1)
                return
            end


            sum_mat = zeros(Int, size(mats[1])...)
            for j in combo
                sum_mat .+= mats[j]
            end

            sgn = signature_matrix(sum_mat)
            key = hash_sgn(sgn)
            lock(lk) do
                if !(key in seen_signatures)
                    push!(seen_signatures, key)
                    push!(unique_signatures, sgn)
                    push!(unique_quad_forms, quad_forms_basis_idx[combo])
                end
            end
        end

        @info "Skipped $(skipped_combos[]) combinations of size $i due to index filter"

        # Save intermediate results
        dir = "data Julia/unique_sgns_$(n)"
        mkpath(dir)
        filename = "size_$(i).csv"
        path = joinpath(dir, filename)
        print_signatures(path, unique_signatures, unique_quad_forms)
    end
    return (unique_signatures, unique_quad_forms)
end


@everywhere function compute_signatures_batch(combos, mats, quad_forms_basis_idx)
    local_signatures = Vector{Tuple{Int, Int, Int}}()
    local_forms = Vector{Vector{Vector{Int}}}()
    seen_local = Set{UInt64}()

    for combo in combos
        # all_indices = vcat(quad_forms_basis_idx[combo]...)
        # if length(unique(all_indices)) < length(combo) + 3
        #     continue
        # end

        sum_mat = zeros(Int, size(mats[1])...)
        for j in combo
            sum_mat .+= mats[j]
        end

        sgn = signature_matrix(sum_mat)
        key = hash_sgn(sgn)
        if key in seen_local
            continue
        end

        push!(seen_local, key)
        push!(local_signatures, sgn)
        push!(local_forms, quad_forms_basis_idx[combo])
    end

    # @info "Processed $(length(combos)) combinations, found $(length(local_signatures)) signatures in this batch"
    return local_signatures, local_forms
end


# Function to compute unique signatures of quadrilinear forms
function unique_signatures_quad_forms_dist(
    n::Int,
    min_comb_size::Int=1,
    max_comb_size::Int=4,
    unique_signatures = Vector{Tuple{Int, Int, Int}}(),
    unique_quad_forms = Vector{Vector{Vector{Int}}}();
    BATCH_SIZE::Int = 5_000
)
    wedge_basis_idx = collect(combinations(1:n, 2))
    quad_forms_basis_idx = collect(combinations(1:n, 4))
    # k = length(wedge_basis_idx)
    m = length(quad_forms_basis_idx)

    # Precompute symmetric matrices
    mats = [symm_matrix_quad_form(wedge_basis_idx, q) for q in quad_forms_basis_idx]

    # Signature cache
    seen_signatures = Set{UInt64}()
    for sgn in unique_signatures
        push!(seen_signatures, hash_sgn(sgn))
    end

    for i in min_comb_size:max_comb_size
        @info "Processing combinations of size $i (distributed)"
        combos = collect(combinations(1:m, i))
        # BATCH_SIZE = div(length(combos), nprocs()) # Dynamic batch size based on number of processes
        batches = Iterators.partition(combos, BATCH_SIZE)

        results = pmap(batch -> compute_signatures_batch(batch, mats, quad_forms_basis_idx), batches)

        for (sgn_list, form_list) in results
            for (sgn, combo) in zip(sgn_list, form_list)
                key = hash_sgn(sgn)
                if !(key in seen_signatures)
                    push!(seen_signatures, key)
                    push!(unique_signatures, sgn)
                    push!(unique_quad_forms, combo)
                end
            end
        end
    end

    # Save intermediate results
    dir = "data Julia/unique_sgns_$(n)"
    mkpath(dir)
    filename = "size_$(max_comb_size)__distributed.csv"
    path = joinpath(dir, filename)
    print_signatures(path, unique_signatures, unique_quad_forms)
    return (unique_signatures, unique_quad_forms)
end

n = 8
known_size = 5
max_size = 6


# @profview unique_signatures_quad_forms(n)

file = "data Julia/unique_sgns_$(n)/size_$(known_size)__distributed.csv"
signatures, combos = read_signatures_file(file)

# wedge_basis_idx = collect(combinations(1:n, 2))
# for combo in combos
#     sum_mat = zeros(Int, length(wedge_basis_idx), length(wedge_basis_idx))
#     for quad_form in combo
#         sum_mat .+= symm_matrix_quad_form(wedge_basis_idx, quad_form)
#     end
#     eig_vals = eigvals(Symmetric(sum_mat))
#     println("Eigenvalues for combo $(combo): $eig_vals")
#     sgn = signature_matrix(sum_mat)
#     println("Signature for combo $(combo): $sgn")
# end

# @time unique_signatures_quad_forms(n, known_size + 1, max_size, signatures, combos)
# @time unique_signatures_quad_forms_dist(n, known_size + 1, max_size, signatures, combos)
# @time unique_signatures_quad_forms_dist2(n, known_size + 1, max_size, signatures, combos)

# Perfect benchmarking (around 30 seconds for n=8)
# @time unique_signatures_quad_forms(8,1,4)
@time unique_signatures_quad_forms_dist(8,1,4)
# @time unique_signatures_quad_forms_dist2(8,1,4)

# @time unique_signatures_quad_forms_dist(7,1,7)

# for B in (2_000, 5_000, 10_000, 20_000)
#     println("Trying BATCH_SIZE = $B")
#     GC.gc()
#     @time unique_signatures_quad_forms_dist(n, 6, 6; BATCH_SIZE=B)
# end


# unique_signature_forms = @benchmark unique_signatures_quad_forms(n,n)

# filename = "data Julia/unique_sgns_$(n)_parallel.csv"
# @btime print_signatures(filename, unique_signature_forms[1], unique_signature_forms[2])
