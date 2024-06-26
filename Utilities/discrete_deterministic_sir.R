## Core equations for transitions between compartments:
update(S) <- S - beta * S * I / N
update(I) <- I + beta * S * I / N - gamma * I
update(R) <- R + gamma * I

## Total population size (odin will recompute this at each timestep:
## automatically)
N <- S + I + R

## Initial states:
initial(S) <- 1 - I_ini
initial(I) <- I_ini # will be user-defined
initial(R) <- 0

## User defined parameters - default in parentheses:
I_ini <- user(1e-3)
beta <- user(0.2)
gamma <- user(0.1)
