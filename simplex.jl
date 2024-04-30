using LinearAlgebra, Combinatorics, Printf

@enum Direction Max Min 


mutable struct SimplexTableau
  z_c     ::Array{Float64} # z_j - c_j
  Y       ::Array{Float64} # inv(B) * A
  x_B     ::Array{Float64} # inv(B) * b
  obj     ::Float64        # c_B * x_B
  b_idx   ::Array{Int64}   # indices for basic variables x_B
end

function is_nonnegative(x::Vector)
  return length( x[ x .< 0] ) == 0
end

function is_nonpositive(z::Array)
  return length( z[ z .> 0] ) == 0
end

function initial_BFS(A, b)
  m, n = size(A)

  b_idx = (n-m+1):n
  B = A[:, b_idx]
  x_B = inv(B) * b
  if is_nonnegative(x_B)
    return b_idx, x_B, B
  end

  error("Infeasible")
end

function format_variable(m::Int, variable::Int)
  if variable <= m
    @sprintf("x[%2d]", variable)
  else
    @sprintf("s[%2d]", variable-m)
  end
end

function print_tableau(t::SimplexTableau, direction)
  m, n = size(t.Y)

  hline0 = repeat("-", 6)
  hline1 = repeat("-", 7*n)
  hline2 = repeat("-", 7)
  hline = join([hline0, "+", hline1, "+", hline2])

  println(hline)

  if direction == Max
    coef_var_obj = -1
    coef_obj = 1
  else
    coef_var_obj = 1
    coef_obj = -1
  end

  @printf("%6s|", "")
  for j in 1:length(t.z_c)
    @printf("%6.2f ", coef_var_obj * t.z_c[j])
  end
  @printf("| %6.2f\n", coef_obj * t.obj)

  println(hline)

  for i in 1:m
    @printf("%s |", format_variable(m, t.b_idx[i]))
    for j in 1:n
      @printf("%6.2f ", t.Y[i,j])
    end
    @printf("| %6.2f\n", t.x_B[i])
  end

  println(hline)
end

function pivoting!(t::SimplexTableau)
  m, n = size(t.Y)

  entering, exiting = pivot_point(t)
  @printf("Pivoting: entering = %s, exiting = %s\n", format_variable(m, entering), format_variable(m, t.b_idx[exiting]))

  # Pivoting: exiting-row, entering-column
  # updating exiting-row
  coef = t.Y[exiting, entering]
  t.Y[exiting, :] /= coef
  t.x_B[exiting] /= coef

  # updating other rows of Y
  for i in setdiff(1:m, exiting)
    coef = t.Y[i, entering]
    t.Y[i, :] -= coef * t.Y[exiting, :]
    t.x_B[i] -= coef * t.x_B[exiting]
  end

  # updating the row for the reduced costs
  coef = t.z_c[entering]
  t.z_c -= coef * t.Y[exiting, :]'
  t.obj += coef * t.x_B[exiting]

  # Updating b_idx
  t.b_idx[ findfirst(t.b_idx .== t.b_idx[exiting]) ] = entering
end

function pivot_point(t::SimplexTableau)
  # Finding the entering variable index
  entering = findfirst( t.z_c .> 0)[2]
  if entering == 0
    error("Optimal")
  end

  # min ratio test / finding the exiting variable index
  pos_idx = findall( t.Y[:, entering] .> 0 )
  if length(pos_idx) == 0
    error("Unbounded")
  end
  exiting = pos_idx[ argmin( t.x_B[pos_idx] ./ t.Y[pos_idx, entering] ) ]

  return entering, exiting
end

function initialize(c, A, b, direction)
  c = Array{Float64}(c)
  A = Array{Float64}(A)
  b = Array{Float64}(b)

  m, n = size(A)

  if direction == Max 
    for i in 1:size(c)[1]
      c[i] = -c[i]
    end
  end

  for i in 1:m 
    slack_column = zeros(1, m)
    slack_column[i] = 1
    A = hcat(A, slack_column')
  end
  
  # Finding an initial BFS
  b_idx, x_B, B = initial_BFS(A,b)

  Y = inv(B) * A
  c_B = c[b_idx]
  obj = dot(c_B, x_B)

  # z_c is a row vector
  z_c = zeros(1,n+m)
  n_idx = setdiff(1:n+m, b_idx)
  z_c[n_idx] = c_B' * inv(B) * A[:,n_idx] - c[n_idx]'

  return SimplexTableau(z_c, Y, x_B, obj, b_idx)
end

function is_optimal(t::SimplexTableau)
  return is_nonpositive(t.z_c)
end

function simplex_method(c, A, b, direction)
  m, _ = size(A)
  c = vcat(c, zeros(1,m)')

  tableau = initialize(c, A, b, direction)
  print_tableau(tableau, direction)

  while !is_optimal(tableau)
    pivoting!(tableau)
    print_tableau(tableau, direction)
  end

  opt_x = zeros(length(c))
  opt_x[tableau.b_idx] = tableau.x_B

  return opt_x, tableau.obj
end

