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

# Construct quad_forms and matrices
function q_form(q::Vector{Int}, b1::Vector{Int}, b2::Vector{Int})
    union_b = union(b1, b2)
    if length(union_b) != 4 || sort(union_b) != sort(q)
        return 0
    end
    x = alt_sum(map(i -> i in b1, q)) != 0
    return (-1)^x
end

function alt_sum(x::Vector{Bool})
    s = 0
    for (i, xi) in enumerate(x)
        s += xi ? (iseven(i) ? -1 : 1) : 0
    end
    return s
end


function symm_matrix(wedge_basis_idx::Vector{Vector{Int}}, q::Function)
    k = length(wedge_basis_idx)
    mat = zeros(Int, k, k)
    for i in 1:k, j in 1:i
        val = q(wedge_basis_idx[i], wedge_basis_idx[j])
        mat[i, j] = val
        mat[j, i] = val
    end
    return sparse(mat)
end


# function symm_matrix_quad_form(wedge_basis_idx::Vector{Vector{Int}}, q::Vector{Int})
#     b1 = q[1,2]
#     b2 = q[3,4]
#     b3 = q[1,3]
#     b4 = q[2,4]
#     row1 = findfirst(x -> x == b1, wedge_basis_idx)

# end

function signature_matrix(mat::SparseMatrixCSC{Int, Int})
    evals = eigvals(Matrix(mat))
    p = count(x -> x > 1e-8, evals)
    q = count(x -> x < -1e-8, evals)
    r = length(evals) - (p + q)
    return (p, q, r)
end


# Create all quad_forms as closures and precompute their matrices
quad_forms = Vector{Function}(undef, m)
mats = Vector{SparseMatrixCSC{Int, Int}}(undef, m)
for i in 1:m
    q_idx = quad_forms_basis_idx[i]
    quad_forms[i] = (b1, b2) -> q_form(q_idx, b1, b2)
    mats[i] = symm_matrix(wedge_basis_idx, quad_forms[i])
end

# Signature cache
seen_signatures = Dict{Tuple{Int, Int, Int}, Bool}()

# Output file
filename = "unique_sgns_$(n)_julia.csv"
open(filename, "w") do io
    for i in 1:n
        for combo in combinations(1:m, i)

            # Sum the matrices for the current combination
            sum_mat = sum(mats[combo])
            # sum_mat = copy(mats[combo[1]])
            # for j in 2:length(combo)
            #     sum_mat .+= mats[combo[j]]
            # end

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
