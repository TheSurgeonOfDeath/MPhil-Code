using Combinatorics
using LinearAlgebra
# using SparseArrays
using Printf
using BenchmarkTools
using .Threads
# import Pkg; Pkg.add("ThreadsX")
using ThreadsX

import Distributed # for addprocs() and nprocs()
import CpuId
Distributed.addprocs(4-Distributed.nprocs())


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


function signature_matrix(mat::Matrix{Int})
    evals = eigvals(Symmetric(mat))
    p = count(x -> x > 1e-8, evals)
    q = count(x -> x < -1e-8, evals)
    r = length(evals) - (p + q)
    return (p, q, r)
end

function print_signatures(filename::String, signatures::Vector{Tuple{Int, Int, Int}}, quad_forms_basis_idx::Vector{Vector{Vector{Int}}})
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

function hash_sgn(sgn::NTuple{3, Int})
    # Convert the signature tuple to a 32-bit integer hash
    # This is a simple hash function that combines the three integers
    # using bitwise operations
    # Ensure the signature is a tuple of three integers
    # assumes sgn[i] < 256 i.e. n < 23 as 23 choose 2 = 253
    return sgn[1] + (sgn[2] << 8) + (sgn[3] << 16)
    # or just use Set{Tuple{Int,Int,Int}} or hash(sgn)
end

# Function to compute unique signatures of quadrilinear forms
function unique_signatures_quad_forms(n::Int, max_comb_size::Int=4)
    wedge_basis_idx = collect(combinations(1:n, 2))
    quad_forms_basis_idx = collect(combinations(1:n, 4))
    # k = length(wedge_basis_idx)
    m = length(quad_forms_basis_idx)

    # Precompute symmetric matrices
    mats = [symm_matrix_quad_form(wedge_basis_idx, q) for q in quad_forms_basis_idx]

    # Signature cache
    seen_signatures = Set{Int}()

    # Find unique signatures
    unique_signatures = Vector{Tuple{Int, Int, Int}}()
    unique_signatures_quad_forms = Vector{Vector{Vector{Int}}}()
    for i in 1:max_comb_size
        @info "Processing combinations of size $i"
        combos = collect(combinations(1:m, i))

        # Parallel map
        results = ThreadsX.mapreduce(vcat, combos; init=Vector{Tuple{UInt64, Tuple{Int, Int, Int}, Vector{Vector{Int}}}}()) do combo
            sum_mat = zeros(Int, size(mats[1])...)
            for j in combo
                sum_mat .+= mats[j]
            end

            tr = sum(diag(sum_mat))
            if abs(tr) < 1
                return []
            end

            sgn = signature_matrix(sum_mat)
            key = hash_sgn(sgn)
            return [(key, sgn, quad_forms_basis_idx[combo])]
        end

        # Deduplicate globally based on hash
        for (key, sgn, combo) in results
            if !(key in seen_signatures)
                push!(seen_signatures, key)
                push!(unique_signatures, sgn)
                push!(unique_quad_forms, combo)
            end
        end

        # Save intermediate results
        filename = "data Julia/unique_sgns_$(n)/size_$(i).csv"
        print_signatures(filename, unique_signatures, unique_signatures_quad_forms)
    end
    return (unique_signatures, unique_signatures_quad_forms)
end

n = 7
# @profview unique_signatures_quad_forms(n)
@time unique_signatures_quad_forms(n,n)
@time unique_signatures_quad_forms2(n,n)
# unique_signature_forms = @benchmark unique_signatures_quad_forms(n,n)


# filename = "data Julia/unique_sgns_$(n)_parallel.csv"
# @btime print_signatures(filename, unique_signature_forms[1], unique_signature_forms[2])
