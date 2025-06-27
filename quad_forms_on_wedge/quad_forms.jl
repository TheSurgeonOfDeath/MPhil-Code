using Combinatorics
using LinearAlgebra
using SparseArrays
using Printf
using BenchmarkTools


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
#     return Matrix(sparse(rows, cols, vals, length(wedge_basis_idx), length(wedge_basis_idx)))
# end

function symm_matrix_quad_form(wedge_basis_idx::Vector{Vector{UInt8}}, q::Vector{UInt8})
    elem_idx = collect(combinations(1:4, 2))
    k = length(wedge_basis_idx)
    mat = zeros(Int8, k, k)
    vals = [1, -1, 1]
    for i in 1:3
        row = findfirst(x -> x == q[elem_idx[i]], wedge_basis_idx)
        col = findfirst(x -> x == q[elem_idx[7 - i]], wedge_basis_idx)
        mat[row, col] = vals[i]
        mat[col, row] = vals[i]
    end
    return mat
end


function signature_matrix(mat::Matrix{Int8, Int8})
    evals = eigvals(Symmetric(mat))
    p = count(x -> x > 1e-8, evals)
    q = count(x -> x < -1e-8, evals)
    r = length(evals) - (p + q)
    return (p, q, r)
end

function print_signatures(filename::String, signatures::Vector{Tuple{UInt8, UInt8, UInt8}}, quad_forms_basis_idx::Vector{Vector{Vector{UInt8}}})
     # Output file
    open(filename, "w") do io
        for (sgn,combo) in zip(signatures, quad_forms_basis_idx)
            @printf(io, "%d,%d,%d", sgn...)
            for q in combo
                @printf(io, ",[%s]", join(q, " "))
            end
            write(io, "\n")
        end
    end
end

function hash_sgn(sgn::NTuple{3, UInt32})
    # Convert the signature tuple to a 32-bit integer hash
    # This is a simple hash function that combines the three integers
    # using bitwise operations
    # Ensure the signature is a tuple of three integers
    # assumes sgn[i] < 256 i.e. n < 23 as 23 choose 2 = 253
    return sgn[1] + (sgn[2] << 8) + (sgn[3] << 16)
    # or just use Set{Tuple{Int,Int,Int}} or hash(sgn)
end

# # Function to compute unique signatures of quadrilinear forms
# function unique_signatures_quad_forms(n::Int)
#     wedge_basis_idx = collect(combinations(1:n, 2))
#     quad_forms_basis_idx = collect(combinations(1:n, 4))
#     # k = length(wedge_basis_idx)
#     m = length(quad_forms_basis_idx)

#     # Precompute symmetric matrices
#     mats = [symm_matrix_quad_form(wedge_basis_idx, q) for q in quad_forms_basis_idx]

#     # Signature cache
#     seen_signatures = Dict{Tuple{Int, Int, Int}, Bool}()

#     # Find unique signatures
#     unique_signatures = Vector{Tuple{Int, Int, Int}}()
#     unique_signatures_quad_forms = Vector{Vector{Vector{Int}}}()
#     for i in 1:n
#         for combo in combinations(1:m, i)
#             current_quad_forms = quad_forms_basis_idx[combo]
#             # Skip non-unique combinations
#             # num_common_elements = length(intersect(current_quad_forms...))
#             # if num_common_elements == 3
#             #     continue
#             # end

#             # Sum the matrices for the current combination
#             sum_mat = spzeros(Int, size(mats[1])...)
#             for idx in combo
#                 sum_mat .+= mats[idx]
#             end

#             # if unique signature, store it
#             sgn = signature_matrix(sum_mat)
#             if !haskey(seen_signatures, sgn)
#                 seen_signatures[sgn] = true
#                 push!(unique_signatures, sgn)
#                 push!(unique_signatures_quad_forms, current_quad_forms)
#             end
#         end
#     end
#     return (unique_signatures, unique_signatures_quad_forms)
# end


# Function to compute unique signatures of quadrilinear forms
function unique_signatures_quad_forms(n::UInt8, max_comb_size::UInt8=4)
    wedge_basis_idx = collect(combinations(1:n, 2))
    quad_forms_basis_idx = collect(combinations(1:n, 4))
    # k = length(wedge_basis_idx)
    m = length(quad_forms_basis_idx)

    # Precompute symmetric matrices
    mats = [symm_matrix_quad_form(wedge_basis_idx, q) for q in quad_forms_basis_idx]

    # Signature cache
    seen_signatures = Set{UInt64}()


    # Find unique signatures
    unique_signatures = Vector{Tuple{UInt8, UInt8, UInt8}}()
    unique_signatures_quad_forms = Vector{Vector{Vector{Int}}}()
    sum_mat = zeros(Int, size(mats[1])...) # preallocate sum matrix
    Threads.@threads for idx in 1:max_comb_size
        @info "Processing combinations of size $i"
        combos = collect(combinations(1:m, i))
        lock = ReentrantLock()

        Threads.@threads for idx in eachindex(combos)
            combo = combos[idx]
            local sum_mat = zeros(Int, size(mats[1])...)  # Initialize sum_mat for each thread

            # Sum matrices for the current combination
            for j in combo
                sum_mat .+= mats[j]  # Accumulate matrices for the current combination
            end

            # Optional early invariant filter
            tr = trace(sum_mat)
            if abs(tr) < 1  # Skip matrices with very low trace (likely degenerate)
                continue
            end

            sgn = signature_matrix(sum_mat)
            key = hash_sgn(sgn)

            lock(lock) do
                if !(key in seen_signatures)
                    push!(seen_signatures, key)
                    push!(unique_signatures, sgn)
                    push!(unique_signatures_quad_forms, quad_forms_basis_idx[combo])
                end
            end
        end

        # Save intermediate results
        filename = "data Julia/unique_sgns_$(n)_size_$(i)_julia.csv"
        print_signatures(filename, unique_signatures, unique_quad_forms)
    end
    return (unique_signatures, unique_signatures_quad_forms)
end

n = 8
# @profview unique_signatures_quad_forms(n)
unique_signature_forms = @benchmark unique_signatures_quad_forms(n)

filename = "data Julia/unique_sgns_$(n)_julia.csv"
@benchmark print_signatures(filename, unique_signature_forms[1], unique_signature_forms[2])
