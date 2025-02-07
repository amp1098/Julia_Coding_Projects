using DifferentialEquations, Plots

alpha = 50
beta = 2
gamma = 50
delta = 1

# u1 is the prey population density
# u2 is the predator population density

function LotkaVolterra(du, u, p, t)
    du[1] = alpha * u[1] - beta * u[1] * u[2]
    du[2] = -gamma * u[2] + delta * u[1] * u[2]
end

u0 = [30, 30]
tspan = (0.0, 50.0)
prob = ODEProblem(LotkaVolterra, u0, tspan)

sol = solve(prob, saveat=1e-3)


# plot(sol, title = "Prey/Predators Population Density", label = ["Prey" "Predator"])

s1 = [sol.u[i][1] for i in range(1, length(sol.u))]
s2 = [sol.u[i][2] for i in range(1, length(sol.u))]

plot(s1, s2, title = "Phase Space", xlim=(0,alpha*2), ylim=(0,gamma*2))