include("simplex.jl")

c = [ 5, 3 ]
A = [
    1 1
    0 1
    2 3
]
s = [ '≤', '=', '≥' ]
b = [ 4, 1, 6 ]

x, z = simplex_method(:Max, c, A, s, b)

@printf("z = %.3f\n", z)
for (i, value) in enumerate(x)
    @printf("x[%d] = %.3f\n", i, value)
end
