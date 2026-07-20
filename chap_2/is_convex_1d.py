import numpy as np

def is_convex_1d(f_str, a, b, n_samples=1000, tol=1e-8):
    """
    Numerically verify the convexity of a one-dimensional function on the interval [a, b].
    Input:
        f_str is the function
        [a,b] is the interval
        n_samples is the number of samples
        tol is the tolerance for floating-point comparison
    Output:
        (is_convex, violations, max_violation)
        is_convex: bool
        violations: the number of violations
        max_violation: the maximum violation value 
    """
    # Samples
    np.random.seed(66)
    x = np.random.uniform(a, b, n_samples)
    y = np.random.uniform(a, b, n_samples)
    theta = np.random.uniform(0.0, 1.0, n_samples)

    mid_point = theta * x + (1 - theta) * y

    # namespace for eval
    namespace = {'x': 0, 'sin': np.sin, 'exp': np.exp}

    # Define the function
    def f(x_arr):
        namespace['x'] = x_arr
        return eval(f_str, {"__builtins__": None}, namespace)
    
    # Compute both sides of the inequality
    lhs = f(mid_point)
    rhs = theta * f(x) + (1 - theta) * f(y)

    # Count the violations and find maximum violation
    diff = lhs - rhs
    violation_mask = diff > tol
    violations = np.sum(violation_mask)
    max_violation = np.max(diff) if violations > 0 else 0.0

    is_convex = (violations == 0)
    return bool(is_convex), int(violations), float(max_violation)

# Test cases
test_case = [
    ("x ** 2", -5, 5),
    ("x ** 3", -5, 5),
    ("exp(x)", -10, 2),
    ("sin(x)", 0, np.pi),
    ("-exp(-x ** 2)", -1.5, 1.5)
]

print(f"{'Expr':<15} {'Interval':<20} {'Convex':<10} {'Violations':<12} {'Max Viol':<20}")
for expr, a, b in test_case:
    is_cvx, viol, max_v = is_convex_1d(expr, a, b)
    interval = f"[{a}, {b}]"
    print(f"{expr:<15} {interval:<20} {is_cvx:<10} {viol:<12} {max_v:<20}")
    
