# # import Pkg; Pkg.add("Symbolics")
# # import Symbolics
# using Symbolics
# @variables x,y,z,g(..);

# import Pkg; Pkg.add("SymbolicUtils")
# using SymbolicUtils
# using SymbolicUtils: Add, Mul, Term, Symbolic

using SymbolicUtils, SymbolicUtils.Rewriters
# using BenchmarkTools

const f = SymbolicUtils.Sym{SymbolicUtils.FnType{Tuple, Number}}(:f)
r1 = @rule(f(~~a, +(~~b), ~~c) => sum([f(~~a..., x, ~~c...) for x in ~(~b)]))
r2 = @rule(f(~~a, ~c::(c -> c isa Integer || typeof(c) <: SymbolicUtils.BasicSymbolic{Int64}) * ~x, ~~b) => ~c * f(~~a..., ~x, ~~b...))
expand_multilinear = Fixpoint(Postwalk(Chain([r1, r2])))

r3 = @rule(f(~~a, ~x, ~y, ~~c) => f(~~a, ~y, ~x, ~~c))
# @syms u v w x y z
# expand_multilinear(f(u+2v, 3w+4x, 5y+6z))
@syms x b::Int y d::Int x1 b1::Int y1 d1::Int

f(2*x,3*y)
r3(f(2*x,3*y))
expand_multilinear(f(2*x,3*y))
f(b*x,3*y)
expand_multilinear(f(b*x,3*y))

typeof(x)
typeof(b)
x isa Integer || typeof(x) <: SymbolicUtils.BasicSymbolic{Number}
typeof(b) <: SymbolicUtils.BasicSymbolic{Int64}

exp1 = expand_multilinear(f(x + y, y1 + x1, b1*y + d*x1 + d1*x + b*y1))
exp2 = expand_multilinear(f(x + x1,y+y1, d1*x1 + b1*y1 + d*x + b*y))
exp1 + exp2

println(exp1)
println(exp2)
# using Symbolics
# Symbolics.variable(:x, T=Integer)

# f(3*a)
# expand_multilinear(f(3*a))
# expand_multilinear(f(x*a))
# typeof(x) <: SymbolicUtils.BasicSymbolic{Number}
# typeof(x)

@variables x y
typeof(x)
expand_multilinear(f(2*a,x*a))

@syms a b a1::Int b1::Int
expand_multilinear(f(a,b,b1*a + a1*b))
