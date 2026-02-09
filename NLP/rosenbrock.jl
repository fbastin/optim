"""
    rosenbrock(x) 

Calculate the Rosenbrock function value at point x.

The Rosenbrock function is a non-convex optimization test function defined as:
f(x, y) = (1 - x)² + 100(y - x²)²

It has a global minimum at (1, 1) with f(1, 1) = 0. The function is also known
as Rosenbrock's valley or Rosenbrock's banana function due to its characteristic
curved, narrow valley shape.

# Arguments
- `x::Vector`: Input vector of length 2, where x[1] and x[2] represent the 
  coordinates in the optimization space

# Returns
- `Float64`: The function value at point x
"""
function rosenbrock(x::Vector)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end

"""
    rosenbrock_gradient!(storage, x)

Compute the gradient of the Rosenbrock function at point x.

The gradient is calculated as:
∂f/∂x = -2(1 - x) - 400x(y - x²)
∂f/∂y = 200(y - x²)

# Arguments
- `storage::Vector`: Pre-allocated vector of length 2 where the gradient will be stored
- `x::Vector`: Input vector of length 2 representing the current point

# Modifies
- `storage[1]`: Partial derivative with respect to x[1]
- `storage[2]`: Partial derivative with respect to x[2]
"""
function rosenbrock_gradient!(storage::Vector, x::Vector)
    # ∂f/∂x: derivative with respect to first variable
    storage[1] = -2.0 * (1.0 - x[1]) - 400.0 * (x[2] - x[1]^2) * x[1]
    # ∂f/∂y: derivative with respect to second variable
    storage[2] = 200.0 * (x[2] - x[1]^2)
end

"""
    rosenbrock_hessian!(storage, x)

Compute the Hessian matrix (second derivatives) of the Rosenbrock function at point x.

The Hessian is a 2×2 symmetric matrix:
H = [∂²f/∂x²    ∂²f/∂x∂y]
    [∂²f/∂y∂x   ∂²f/∂y²  ]

# Arguments
- `storage::Matrix`: Pre-allocated 2×2 matrix where the Hessian will be stored
- `x::Vector`: Input vector of length 2 representing the current point

# Modifies
- `storage[1,1]`: Second partial derivative ∂²f/∂x² = 2 - 400y + 1200x²
- `storage[1,2]`: Mixed partial derivative ∂²f/∂x∂y = -400x
- `storage[2,1]`: Mixed partial derivative ∂²f/∂y∂x = -400x (symmetric)
- `storage[2,2]`: Second partial derivative ∂²f/∂y² = 200
"""
function rosenbrock_hessian!(storage::Matrix, x::Vector)
    # ∂²f/∂x²: second derivative with respect to x[1]
    storage[1, 1] = 2.0 - 400.0 * x[2] + 1200.0 * x[1]^2
    # ∂²f/∂x∂y: mixed partial derivative
    storage[1, 2] = -400.0 * x[1]
    # ∂²f/∂y∂x: mixed partial derivative (Hessian is symmetric)
    storage[2, 1] = -400.0 * x[1]
    # ∂²f/∂y²: second derivative with respect to x[2]
    storage[2, 2] = 200.0
end