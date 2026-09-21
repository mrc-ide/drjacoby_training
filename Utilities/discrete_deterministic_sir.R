## Core equations for transitions between compartments:
update(S) <- S - beta * S * I / N
update(I) <- I + beta * S * I / N - gamma * I
update(R) <- R + gamma * I

## Total population size (odin will recompute this at each timestep:
## automatically)
N <- S + I + R

## Initial states:
initial(S) <- 1 - I_init
initial(I) <- I_init # will be user-defined
initial(R) <- 0

## User defined parameters - default in parentheses:
I_init <- user(1e-3)
beta <- user(0.2)
gamma <- user(0.1)
