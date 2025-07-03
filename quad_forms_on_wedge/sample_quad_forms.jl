using LinearAlgebra
using Combinatorics
using StatsBase
using BenchmarkTools


# Compute sign of permutation needed to bring (i,j,k,l) to sorted order
function wedge_sign(i, j, k, l)
    v = [i, j, k, l]
    inv = 0
    for a in 1:3, b in a+1:4
        if v[a] > v[b]
            inv += 1
        end
    end
    return isodd(inv) ? -1 : 1
end

# Evaluate 4-form q on basis 4-tuple (i,j,k,l) using q_dict
function eval_q(q_dict, i,j,k,l)
    inds = sort([i,j,k,l])
    sgn = wedge_sign(i,j,k,l)
    return sgn * get(q_dict, inds, 0.0)
end

# Build the matrix of B_q for a given 4-form q_dict
function build_Bq(q_dict)
    dim = length(Λ2_basis)
    Bq = zeros(Float64, dim, dim)

    for (r, (i,j)) in enumerate(Λ2_basis)
        for (s, (k,l)) in enumerate(Λ2_basis)
            Bq[r,s] = eval_q(q_dict, i,j,k,l)
        end
    end

    # Enforce symmetry explicitly
    return 0.5 * (Bq + Bq')
end

# Compute signature (n+, n−, n0) of symmetric matrix A
function signature(A; tol=1e-8)
    evals = eigvals(Symmetric(A))
    npos = count(x -> x > tol, evals)
    nneg = count(x -> x < -tol, evals)
    nzero = length(evals) - npos - nneg
    return (npos, nneg, nzero)
end

# Generate a random 4-form as a Dict{NTuple{4,Int},Float64}
function random_4form()
    q = Dict{Vector{Int64},Float64}()
    for idx in Λ4_basis
        q[idx] = randn()
    end
    return q
end


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

# Sample many 4-forms and collect observed signatures
function collect_signatures(num_trials::Int)
    signatures = []
    for _ in 1:num_trials
        q = random_4form()
        Bq = build_Bq(q)
        sig = signature(Bq)
        push!(signatures, sig)
    end
    return countmap(signatures)
end

function random_qform_matrix()
    coeffs = randn(length(Λ4_basis_matrices))
    q = zeros(size(Λ4_basis_matrices[1]))
    for i in eachindex(Λ4_basis_matrices)
        q += coeffs[i] * Λ4_basis_matrices[i]
    end
    return q
end

function collect_signatures2(num_trials::Int)
    signatures = []
    for _ in 1:num_trials
        Bq = random_qform_matrix()
        sig = signature(Bq)
        push!(signatures, sig)
    end
    return countmap(signatures)
end


# Define basis of Λ²(ℝ⁷)
n = 7
Λ2_basis = collect(combinations(1:n, 2))  # index pairs (i,j), i < j
Λ4_basis = collect(combinations(1:n, 4))  # index tuples (i,j,k,l), i<j<k<l

# Map: (i,j) ↦ basis index in Λ²(ℝ⁷)
Λ2_index = Dict((i,j) => idx for (idx, (i,j)) in enumerate(Λ2_basis))

Λ4_basis_matrices = [symm_matrix_quad_form(Λ2_basis, q) for q in Λ4_basis]



# q = random_4form()
# Bq = build_Bq(q)
# sig = signature(Bq)
# push!(signatures, sig)
# sig_counts = countmap(signatures)

@benchmark collect_signatures(1_000)
@benchmark collect_signatures2(1_000)


# Example: run 1000 trials
sig_counts = collect_signatures(1_000_000)
for (sig, count) in sig_counts
    println("Signature $sig occurred $count times.")
end
