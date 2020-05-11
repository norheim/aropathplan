using Random
include("dynamics.jl")
include("solver.jl")
include("detsolver.jl")
include("draw.jl")
include("utils.jl")
include("polygonlib.jl")

#polygons = [[(0.015, 0), (0.015, 0.03), (0.035, 0.03),(0.035,0)]]
# P = [-1 0; 0 1]
# q = [-0.015; 0.03]
polygons = blackmore1
P, q = makePq(polygons)
obj = [[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16],[17,18,19,20]]
N = 10
npastK = 2
T = 20

x0 = [0;0;0;0]
xN = [0;0;10;0]

ϵ2 = 0.01 # Probability of failure
max_wind = 0.01
ρ = max_wind*sqrt(2*log(1/ϵ2))

x_k_det, simulate_stoc, Ai, B, Bw = get_dynamic_matrices(T, N)
x_k_det_x0 = (k, α) -> x_k_det(x0, k, α)
problem = (Ai, Bw, B, N, x_k_det_x0, P, q, obj, xN)
ϵmin = pathplanner((problem..., 0, ρ), true)
ϵ = 2*ϵmin
ϵ[[1,3]] = 7.5/2*ϵmin[[1,3]]
xdet, ydet = detpathplanner((problem...,  ϵ, ρ))

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

pl = plot(aspect_ratio=:equal)
plot_polygons(polygons)
plot_goal(xN, ϵ)
plot!(xdet, ydet, linestyle=:dash, markershape=:circle, ms=2)
plot!(posx[:,:]', posy[:,:]', alpha=0.05, linestyle=:dash, markershape=:circle, ms=2)
display(pl)
