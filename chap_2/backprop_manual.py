import math

def f_and_grad(x_1,x_2):
    # Forward Pass
    omega_1 = x_1
    omega_2 = x_2
    omega_3 = omega_1 * omega_2
    omega_4 = math.sin(omega_1)
    omega_5 = omega_3 + omega_4
    f = omega_5

    # Backward Pass
    df_domega_5 = 1
    df_domega_4 = df_domega_5 * 1
    df_domega_3 = df_domega_5 * 1
    df_domega_2 = df_domega_3 * omega_1
    df_domega_1 = df_domega_3 * omega_2 + df_domega_4 * math.cos(omega_1)

    return f, df_domega_1, df_domega_2

if __name__ == "__main__":
    x_1 = math.pi / 2
    x_2 = 1
    f, df_dx1, df_dx2 = f_and_grad(x_1, x_2)
    print(f"f({x_1}, {x_2}) = {f}")
    print(f"df/dx1 = {df_dx1}")
    print(f"df/dx2 = {df_dx2}")