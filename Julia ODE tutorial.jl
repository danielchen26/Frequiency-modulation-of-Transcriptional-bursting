#=
The following is the simulation test for the paper " Frequency Modulation of Transcriptional Bursting Enables Sensitive and Rapid Gene Regulation "
=#

# ---- Importing modules

using DynamicalSystems
using DifferentialEquations

# -------------------------------------
# ---- Defining scalar equation problem
# -------------------------------------

f(u, p, t) = 1.01* u
u0 = 0.5
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)

# ---- Solving problem

# 1.
sol = solve(prob)
# 2. control lower relative tolerance using "reltol"
sol = solve(prob, reltol = 1e-6)
# 3. control the solver save time points
sol = solve(prob, reltol = 1e-6, saveat = 0.1)
# 4. control the saved points and speed up solution. Only save the endpoint.
sol = solve(prob, reltol = 1e-6, save_everystep = false)


# ---- Choosing a solver algorithm

# 1. provide a hint for solver to choose the right algorithm
# if we have a stiff problem where we need high accuracy, but don't know the best stiff algorithm for this problem, we can use:
sol = solve(prob, alg_hints = [:stiff], reltol = 1e-8, abstol = 1e-8)

# 2. to choose a specific algorithm to use, EX. choose 5th order Tsitouras method:
sol = solve(prob, Tsit5())

# 3. solver controls can be combined with algorithm choice. We can solve problem with Tsit5() with lower tolerance via:
sol = solve(prob, Tsit5(), reltol =1e-8, abstol =1e-8)





#=
some good "go-to" choices for ODEs are:
* AutoTsit5(Rosenbrock23()) : handles both stiff and non-stiff equations. This is a good algorithm to use if you know nothing about the equation

* BS3() for fast low accuracy non-stiff.[5]

* Tsit5() for standard non-stiff. This is the first algorithm to try in most cases

* Vern7() for high accuracy non-stiff

* Rodas4() for stiff equations with Julia-defined types, events, etc.

* radau() for really high accuracy stiff equations (requires installing ODEInterfaceDiffEq.jl)
=#




# ---- Analyzing the solution

# 1. The 5th value of the solution:
sol[5]

# 2. 8th timestep:
sol.t[8]

# 3. Build an array :
[t + u for (u,t) in tuple(sol)]# or
[t + 2u for (u,t) in zip(sol.u, sol.t)]

# 4. We can access the interpolated values by treating sol as a function:
sol(0.45)




# ---- Plotting solutions

using Plots # to open an interactive GUi, type gui()
plot(sol)

plot(sol, linewidth = 5, title  = " Solution to linear ODE with thick line", xaxis = "Time(t)", yaxis = "u(t) (in μm)", label = "My Thick Line!")
# Add another plot on top of the above Plot
plot!(sol.t, t->0.5*exp(1.01t), lw = 3, ls =:dash, label = "True solution")





# -------------------------------------
# Solving systems of equations
# -------------------------------------
#=
Defining your ODE function to be in-place updating can have performance benefits:

Instead of writing a function which outputs its solution, you write a function which updates a vector that is designated to hold the solution.

DifferentialEquations.jl's solver packages are able to reduce the amount of array allocations and achieve better performance.

=#

# Defining Lorenz function
function lorenz(du, u, p, t)
    du[1] = 10.0* (u[2] - u[1])
    du[2] = u[1]* (28.0 - u[3]) -u[2]
    du[3] = u[1]* u[2] - (8/3)* u[3]
end

# Using lorenz function in a probelm
u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz, u0, tspan)
sol = solve(prob)

# We can choose to do a 3D phase space plot between the different variables:
plot(sol, vars = (1,2,3))
# We can plot the time series just for the second component(note that "variable 0" is time):
plot(sol, vars=(0,2))



# -------------------------------------
# Defining Parameterized Functions
# -------------------------------------

function parameterized_lorenz(du,u,p,t)
 du[1] = p[1]*(u[2]-u[1])
 du[2] = u[1]*(p[2]-u[3]) - u[2]
 du[3] = u[1]*u[2] - p[3]*u[3]
end

# we add the parameters to the ODEProblem:
u0 = [1.0;0.0;0.0]
tspan = (0.0,1.0)
p = [10.0,28.0,8/3]
prob = ODEProblem(parameterized_lorenz,u0,tspan,p)


#=
There exists a @ode_def macro allows for "defining your ODE in pseudocode" and getting a function which is efficient and runnable:

To use the macro, you write out your system of equations with the left-hand side being d_ and those variables will be parsed as the dependent variables. The independent variable is t, and the other variables are parameters which you pass at the end
=#

# Rewrite lorenz system using @ode_def macro

g = @ode_def LorenzExample begin
    dx = σ*(y-x)
    dy = x*(ρ-z) - y
    dz = x*y -β*z
end σ ρ β

u0 = [1.0;0.0;0.0]
tspan = (0.0,1.0)
p = [10.0,28.0,8/3]
prob = ODEProblem(g,u0,tspan,p)
