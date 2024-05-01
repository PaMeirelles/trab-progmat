include("simplex.jl")

c = [ 4, 1 ]
A = [
    4 3
    1 2
]
s = [ '≥', '≤' ]
b = [ 6, 4 ]

x, z = simplex_method_two_stage(c, A, s, b)

@printf("z = %.3f\n", z)
for (i, value) in enumerate(x)
    @printf("x[%d] = %.3f\n", i, value)
end
