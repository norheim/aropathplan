using Random, Plots
include("dynamics.jl")
include("solver.jl")
include("draw.jl")
include("utils.jl")

polygons = [[(0.9523419261326787, 1.2051055699042035),(1.1513383271320627, 3.0614593163069603),(2.9865055038717943, 2.8621257285224893),(2.84258294443677, 0.9874795990287133)],[(1.4262579525737413, 3.2640242774485237),(-1.3456141131402157, 4.068355254376987),(-0.3956746428512883, 6.807941685418777),(2.284510716217514, 5.9667168544219455)],[(-2.114905763318293, 3.1850944708391964),(-4.407663506321082, 2.64865839466308),(-4.924879935747185, 4.9082032393347115),(-2.5770483511800366,5.426346932419806)],[(-2.9636891980808455, 6.289460805548596),(-5.753741251290793, 6.983503604854553),(-5.0240409494617895, 9.759506208685838),(-2.3074113187047347, 9.083727693572143)],[(0.08071408668041968, 7.24971257168146),(1.2142517614497157, 9.916073272561658),(3.8581895411215674, 8.78076536717065),(2.7429442494432914, 6.151186124763792)]]
P = [-1 0; 0 1]
q = [-0.015; 0.03]
#P,q = makePq(polygons);
obj = [[1,2]]
N = 40
T = 1

x0 = [0;0;0;0]
xN = [0;0;10;0]

ϵ2 = 0.01 # Probability of failure
max_wind = 0.1
ρ = max_wind*sqrt(2*log(1/ϵ2))

x_k_det, simulate_stoc, Ai, B, Bw = get_dynamic_matrices(T, N)
x_k_det_x0 = (k, α) -> x_k_det(x0, k, α)
problem = (Ai, Bw, B, N, x_k_det_x0, P, q, obj, xN)
ϵmin = pathplanner((problem..., 0, ρ), true)
ϵ = 1.1*ϵmin
K_out, α_out = pathplanner((problem..., ϵ, ρ), false)

iterations = 1000
Γ = 1
posx = zeros(iterations, N)
posy = zeros(iterations, N)
for i = 1:iterations
        w = Γ * max_wind * randn(Float64, (2, N - 1))
        x = simulate_stoc(x0, α_out, K_out, w)
        posx[i,:] = x[1,:]
        posy[i,:] = x[3,:]
end

faults = check_collisions(posx, posy, P, q, obj)
p_faults = sum(faults)/(N*iterations)

pl = plot()
plot_polygons(polygons)
plot_goal(xN, ϵ)
plot!(posx[:,:]', posy[:,:]', alpha=0.05, linestyle=:dash, markershape=:circle, ms=2)
display(pl)
