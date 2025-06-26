using Combinatorics
using LinearAlgebra
using SparseArrays
using Printf

# Define parameters
const n = 6
const wedge_basis_idx = collect(combinations(1:n, 2))
const quad_forms_basis_idx = collect(combinations(1:n, 4))
const k = length(wedge_basis_idx)
const m = length(quad_forms_basis_idx)

# Construct symmetric matrices

function symm_matrix_quad_form(wedge_basis_idx::Vector{Vector{Int}}, q::Vector{Int})
    elem_idx = collect(combinations(1:4, 2))
    rows = Vector{Int}(undef, 3)
    cols = Vector{Int}(undef, 3)
    for i in 1:3
        rows[i] = findfirst(x -> x == q[elem_idx[i]], wedge_basis_idx)
        cols[i] = findfirst(x -> x == q[elem_idx[7-i]], wedge_basis_idx)
    end

    vals = [1,-1,1]
    vals = [vals; vals]
    rows_tmp = rows
    rows = [rows; cols]
    cols = [cols; rows_tmp]
    return sparse(rows, cols, vals, length(wedge_basis_idx), length(wedge_basis_idx))
end


function signature_matrix(mat::SparseMatrixCSC{Int, Int})
    evals = eigvals(Matrix(mat))
    p = count(x -> x > 1e-8, evals)
    q = count(x -> x < -1e-8, evals)
    r = length(evals) - (p + q)
    return (p, q, r)
end


# Precompute symmetric matrices
mats = [symm_matrix_quad_form(wedge_basis_idx, q) for q in quad_forms_basis_idx]

# Signature cache
seen_signatures = Dict{Tuple{Int, Int, Int}, Bool}()

# Output file
filename = "data Julia/unique_sgns_$(n)_julia.csv"
open(filename, "w") do io
    for i in 1:n
        for combo in combinations(1:m, i)

            # Sum the matrices for the current combination
            sum_mat = spzeros(Int, size(mats[1])...)
            for idx in combo
                sum_mat .+= mats[idx]
            end

            sgn = signature_matrix(sum_mat)
            if !haskey(seen_signatures, sgn)
                seen_signatures[sgn] = true
                @printf(io, "%d,%d,%d", sgn...)
                for idx in combo
                    @printf(io, ",[%s]", join(quad_forms_basis_idx[idx], " "))
                end
                write(io, "\n")
            end
        end
    end
end
