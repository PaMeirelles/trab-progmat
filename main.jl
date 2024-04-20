include("simplex.jl")

using Main.SimplexMethod

c = [ 4, 3 ]
A = [
    2 1
    1 2
]
b = [ 4, 4 ]

simplex_method(c, A, b)
