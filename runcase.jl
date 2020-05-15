using Random
include("dynamics.jl")
include("nksolver.jl")
include("detsolver.jl")
include("draw.jl")
include("utils.jl")
include("polygonlib.jl")

# polygons = [[(0.015, 0), (0.015, 0.03), (0.035, 0.03),(0.035,0)]]
# polygons = [[(0.25, 0), (0.25, 0.75), (1.1, 0.75),(1.1,0)]]
# P = [-1 0; 0 1]
# q = [-0.015; 0.03]
# q = [-0.25; 0.75]
# obj = [[1, 2]]
polygons = blackmore1
obj = blackmore1obj
P, q = makePq(polygons)
N = 10
npastK = 0
T = 20

x0 = [0;0;0;0.5];
xN = [0;0;10;0];
# x0 = [0;0;0;0]
# xN = [1;0;1;0]

ϵ2 = 0.01 # Probability of failure
max_wind = 0.1/N
max_au = 2.5/N
ρ = max_wind*sqrt(2*log(1/ϵ2))

x_k_det, simulate_stoc, Ai, B, Bw = get_dynamic_matrices(T, N)
x_k_det_x0 = (k, α) -> x_k_det(x0, k, α)
problem = (Ai, Bw, B, N, x_k_det_x0, P, q, obj, xN, max_au)
ϵmin = pathplanner((problem..., 0, ρ), npastK, true)
ϵmindet = pathplanner((problem..., 0, 0), 0, true)
#ϵ = 1.1*ϵmin
#ϵ = [0.2, 0.1, 0.2, 0.1]
ϵdet = [1e-3, 1e-1, 1e-3, 7e-1]
ϵ = [1, 1, 1, 1]

xdet, ydet, aud, audi = detpathplanner((problem...,  ϵdet, ρ))

K_out, α_out, aur = pathplanner((problem..., ϵ, ρ), npastK, false)

iterations = 1000
Γ = 1
posx = zeros(iterations, N)
posy = zeros(iterations, N)
aus = zeros(iterations)
for i = 1:iterations
        w = Γ * max_wind * randn(Float64, (2, N - 1))
        x, au = simulate_stoc(x0, α_out, K_out, w)
        aus[i] = sum(au)
        posx[i,:] = x[1,:]
        posy[i,:] = x[3,:]
end

faults = check_collisions(posx, posy, P, q, obj)
p_faults = sum(faults)/(N*iterations)

pl = plot(aspect_ratio=:equal, dpi=160)
plot_polygons(polygons)
plot_goal(xN, ϵ)
plot!(xdet, ydet, linestyle=:dash, markershape=:circle, ms=2)
scatter!(posx[:,:]', posy[:,:]', linestyle=:dash, alpha=0.05,  ms=2)
plot!(mean(posx, dims=1)', mean(posy, dims=1)', linestyle=:dash, markershape=:circle, ms=2)
display(pl)
savefig(pl, "test_plot.png");
