using Random
include("dynamics.jl")
include("nksolver.jl")
include("draw.jl")
include("utils.jl")

polygons = [[(0.015, 0), (0.015, 0.03), (0.035, 0.03),(0.035,0)]]
P = [-1 0; 0 1]
q = [-0.015; 0.03]
#P,q = makePq(polygons);
obj = [[1,2]]
N = 10
npastK = 1
T = 1

x0 = [0;0;0;0]
xN = [0.03;0;0.035;0]

ϵ2 = 0.01 # Probability of failure
max_wind = 0.01
ρ = max_wind*sqrt(2*log(1/ϵ2))

x_k_det, simulate_stoc, Ai, B, Bw = get_dynamic_matrices(T, N)
x_k_det_x0 = (k, α) -> x_k_det(x0, k, α)
problem = (Ai, Bw, B, N, x_k_det_x0, P, q, obj, xN)
ϵmin = pathplanner((problem..., 0, ρ), npastK, true)
ϵ = 1.1*ϵmin
K_out, α_out = pathplanner((problem..., ϵ, ρ), npastK, false)

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
