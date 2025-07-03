using LinearAlgebra
using Combinatorics
using StatsBase
using BenchmarkTools
using Printf


# Compute signature (n+, n−, n0) of symmetric matrix A
function signature_matrix(A; tol=1e-8)
    evals = eigvals(Symmetric(A))
    npos = count(x -> x > tol, evals)
    nneg = count(x -> x < -tol, evals)
    nzero = length(evals) - npos - nneg
    return (npos, nneg, nzero)
end

function symm_matrix_quad_form(Λ2_basis::Vector{Vector{Int}}, q::Vector{Int})
    elem_idx = collect(combinations(1:4, 2))
    k = length(Λ2_basis)
    mat = zeros(Int, k, k)
    vals = [1, -1, 1]
    for i in 1:3
        row = findfirst(x -> x == q[elem_idx[i]], Λ2_basis)
        col = findfirst(x -> x == q[elem_idx[7 - i]], Λ2_basis)
        mat[row, col] = vals[i]
        mat[col, row] = vals[i]
    end
    return mat
end

function random_qform_matrix(Λ4_basis_matrices)
    coeffs = randn(length(Λ4_basis_matrices))
    q_mat = zeros(size(Λ4_basis_matrices[1]))
    for i in eachindex(Λ4_basis_matrices)
        q_mat += coeffs[i] * Λ4_basis_matrices[i]
    end
    return q_mat
end

function collect_signatures(n::Int, num_trials::Int = 1000)
    # Define basis of Λ²(ℝⁿ) and Λ⁴(ℝⁿ)
    Λ2_basis = collect(combinations(1:n, 2))  # index pairs [i,j], i < j
    Λ4_basis = collect(combinations(1:n, 4))  # index vectors [i,j,k,l], i<j<k<l
    Λ4_basis_matrices = [symm_matrix_quad_form(Λ2_basis, q) for q in Λ4_basis]

    signatures = Vector{Tuple{Int, Int, Int}}()
    for _ in 1:num_trials
        Bq = random_qform_matrix(Λ4_basis_matrices)
        sig = signature_matrix(Bq)
        push!(signatures, sig)
    end
    return countmap(signatures)
end

function print_signature_counts(filename::String, sig_counts::Dict{Tuple{Int, Int, Int}}, mode::String = "w")
     # Output file
    open(filename, mode) do io
        for (sgn, count) in sig_counts
            @printf(io, "%d,%d,%d,%d", sgn..., count)
            write(io, "\n")
        end
    end
end


# q = random_4form()
# Bq = build_Bq(q)
# sig = signature(Bq)
# push!(signatures, sig)
# sig_counts = countmap(signatures)


# Example: run 1000 trials
n = 9
num_trials = 100_000
filename = "data Julia/unique_sgns_$(n)/random_samples_$(num_trials).csv"

@time sig_counts = collect_signatures(n, num_trials)
print_signature_counts(filename, sig_counts)
# for (sig, count) in sig_counts
#     println("Signature $sig occurred $count times.")
# end
