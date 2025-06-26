using Combinatorics
using LinearAlgebra
using SparseArrays
using Printf


# Function to compute unique signatures of quadrilinear forms
function unique_signatures_quad_forms(n::Int)
    wedge_basis_idx = collect(combinations(1:n, 2))
    quad_forms_basis_idx = collect(combinations(1:n, 4))
    # k = length(wedge_basis_idx)
    m = length(quad_forms_basis_idx)

    # Precompute symmetric matrices
    mats = [symm_matrix_quad_form(wedge_basis_idx, q) for q in quad_forms_basis_idx]

    # Signature cache
    seen_signatures = Dict{Tuple{Int, Int, Int}, Bool}()

    # Find unique signatures
    unique_signatures = Vector{Tuple{Int, Int, Int}}()
    unique_signatures_quad_forms = Vector{Vector{Vector{Int}}}()
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
                push!(unique_signatures, sgn)
                push!(unique_signatures_quad_forms, quad_forms_basis_idx[combo])
            end
        end
    end

    return (unique_signatures, unique_signatures_quad_forms)
end

function print_signatures(signatures::Vector{Tuple{Int, Int, Int}}, quad_forms_basis_idx::Vector{Vector{Vector{Int}}})
     # Output file
    filename = "data Julia/unique_sgns_$(n)_julia.csv"
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

n = 6
unique_signatures_forms = unique_signatures_quad_forms(n)
print_signatures(unique_signatures_forms[1], unique_signatures_forms[2])
