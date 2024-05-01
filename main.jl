include("simplex.jl")

c = [ -4, -3 ]
A = [
    2 1
    1 2
]
b = [ 4, 4 ]

simplex_method(:Min, c, A, b)
