using LinearAlgebra, Combinatorics, Printf


mutable struct SimplexTableau
  z_c     ::Array{Float64} # z_j - c_j
  Y       ::Array{Float64} # inv(B) * A
  x_B     ::Array{Float64} # inv(B) * b
  obj     ::Float64        # c_B * x_B
  b_idx   ::Array{Int64}   # indices for basic variables x_B
  direction :: Symbol
end

function is_nonnegative(x)
  return length( x[ x .< 0] ) == 0
end

function is_nonpositive(z)
  return length( z[ z .> 0] ) == 0
end

function initial_BFS(A, b, b_idx=nothing)
  m, n = size(A)

  if isnothing(b_idx)
    b_idx = (n-m+1):n
  end
  B = A[:, b_idx]
  inv_b = inv(B)
  x_B = inv_b * b
  if is_nonnegative(x_B)
    return b_idx, x_B, inv_b
  end

  error("Infeasible")
end


function print_tableau(t::SimplexTableau)
  m, n = size(t.Y)

  hline0 = repeat("-", 6)
  hline1 = repeat("-", 7*n)
  hline2 = repeat("-", 7)
  hline = join([hline0, "+", hline1, "+", hline2])

  println(hline)

  @printf("%6s|", "")
  for j in 1:length(t.z_c)
    @printf("%6.2f ", t.z_c[j])
  end
  @printf("| %6.2f\n", t.obj)

  println(hline)

  for i in 1:m
    @printf("x[%d]  |", t.b_idx[i])
    for j in 1:n
      @printf("%6.2f ", t.Y[i,j])
    end
    @printf("| %6.2f\n", t.x_B[i])
  end

  println(hline)
end

function pivoting!(t::SimplexTableau)
  m, _ = size(t.Y)

  entering, exiting = pivot_point(t)
  @printf("Pivoting: entering = %s, exiting = %s\n", entering, t.b_idx[exiting])

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
  t.obj -= coef * t.x_B[exiting]

  # Updating b_idx
  t.b_idx[ findfirst(t.b_idx .== t.b_idx[exiting]) ] = entering
end

function pivot_point(t::SimplexTableau)
  # Finding the entering variable index
  if t.direction == :Min
    entering = findfirst( t.z_c .> 0)[2]
  else 
    entering = findfirst( t.z_c .< 0)[2]
  end
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

function initialize(direction, c, A, b, b_idx, s)
  c = Array{Float64}(c)
  A = Array{Float64}(A)
  b = Array{Float64}(b)

  m, n = size(A)

  for i in 1:m 
    slack_column = zeros(1, m)
    if s[i] == '≥'
      coef = -1
    elseif s[i] == '≤' 
      coef = 1
    else
      continue
    end
    slack_column[i] = coef
    A = hcat(A, slack_column')
  end

  b_idx, x_B, inv_b = initial_BFS(A, b, b_idx)

  B = A[:, b_idx]
  inv_b = inv(B)
  x_B = inv_b * b
  Y = inv_b * A


  c_B = c[b_idx]
  obj = dot(c_B, x_B)

  m, n = size(A)

  # z_c is a row vector
  z_c = zeros(1,n)
  n_idx = setdiff(1:n, b_idx)
  z_c[n_idx] = c_B' * inv_b * A[:,n_idx] - c[n_idx]'

  return SimplexTableau(z_c, Y, x_B, obj, b_idx, direction)
end

function initialize_two_stage(c, A, signals, b)
  c = Array{Float64}(c)
  A = Array{Float64}(A)
  b = Array{Float64}(b)

  m, n = size(A)
  b_idx = [] 
  c = vec(zeros(1,n))
  extra_variables = 0
  indexes_artificials = []

  for i in eachindex(signals)
    s = signals[i] 
    extra_variables += 1
    if s == '≤'
      extra_column = zeros(1, m)
      extra_column[i] = 1
      push!(b_idx, n + extra_variables)
      push!(c, 0)
    elseif s == '='
      extra_column = zeros(1, m)
      extra_column[extra_variables] = 1
      push!(b_idx, n + extra_variables)
      push!(c, -1)
      push!(indexes_artificials, n + extra_variables)
    else
      extra_column = zeros(1, m)
      extra_column[extra_variables] = 1
      push!(b_idx, n + extra_variables)
      push!(c, -1)
      push!(indexes_artificials, n + extra_variables)

      A = hcat(A, extra_column')

      extra_column = zeros(1, m)
      extra_column[extra_variables] = -1
      extra_variables += 1
      push!(c, 0)

    end
    A = hcat(A, extra_column')
  end
  B = A[:, b_idx]
  inv_b = inv(B)
  x_B = inv_b * b
  Y = inv_b * A
  c_B = c[b_idx]
  obj = -dot(c_B, x_B)

  # z_c is a row vector
  z_c = zeros(1,n+extra_variables)
  n_idx = setdiff(1:n+extra_variables, b_idx)
  z_c[n_idx] = c[n_idx]' - c_B' * inv_b * A[:,n_idx]

  return SimplexTableau(z_c, Y, x_B, obj, b_idx, :Min), indexes_artificials
end

function is_optimal(t::SimplexTableau)
  if t.direction == :Min
    return is_nonpositive(t.z_c)
  else
    return is_nonnegative(t.z_c)  
  end
end

function two_stage_is_optimal(t::SimplexTableau)
  return t.obj == 0
end

function subtract_elements(A, B)
  result = similar(A)
  for i in eachindex(A)
      n = sum(B .< A[i])
      result[i] = A[i] - n
  end
  return result
end

function simplex_method_two_stage(c, A, signals, b)
  m, _= size(A)
  c = vcat(c, zeros(1,m)')

  tableau, indexes_artificials = initialize_two_stage(c, A, signals, b)
  print_tableau(tableau)

  while !is_optimal(tableau)
    pivoting!(tableau)
    print_tableau(tableau)
  end

  new_base = subtract_elements(tableau.b_idx, indexes_artificials)

  return new_base
end

function simplex_method(direction, c, A, s, b)
  _, m = size(A)
  if any(x -> '=' in x || '≥' in x, s)
    b_idx = simplex_method_two_stage(c, A, s, b)
  else
    b_idx = nothing
  end
  c = vcat(c, zeros(1,m)')
  tableau = initialize(direction, c, A, b, b_idx, s)
  print_tableau(tableau)

  while !is_optimal(tableau)
    pivoting!(tableau)
    print_tableau(tableau)
  end

  opt_x = zeros(length(c))
  opt_x[tableau.b_idx] = tableau.x_B

  return opt_x, tableau.obj
end

