# Robust plan
using JuMP, Gurobi, Plots, LinearAlgebra, Random
GRB_ENV = Gurobi.Env();

# Create model for dynamics
dt = 0.1
A = zeros(4,4)+I
A[1,2] = A[3,4] = dt
B = zeros(4,2)
B[2,1] = B[4,2] = dt
C = zeros(2,4)
C[1,2] = C[2,4] = 1
D = zeros(2,2)
Bw = copy(B);

# Create model
T = 0.9
N = Int(T/dt);
# Choose n-K parameter(how many past "winds" we can use in our control)
npastK = N-2
x0 = [0;0;0;0]
xN = [0.03;0;0.035;0]
ϵ = [1e-4;
     0.00255;
     1e-4;
     0.00255];
    # how close to the goal we want to get

# Obstacles
# P = [-1 0; 0 1]
# q = [-0.03; 0.03] # this should be feasible
P = [-1 0; 0 1]
q = [-0.015; 0.03]
m = length(q)
obj = [[1,2]] # Each entry corresponds to the lines defining an object
Nobj = length(obj)
M = 100; #Big M

# Helper variables
Q = zeros(2,4)
Q[1,1] = Q[2,3] = -1
N_K = Int((N-2)*(N-1)/2) # number of K's
#N_s = Int((N-3)*(N-2)/2)
Ak = zeros(N+1,4,4)
Ak[1,:,:] = zeros(Int, size(A))+I
for i=2:N+1
    Ak[i,:,:] = A*Ak[i-1,:,:]
end
Ai = i -> Ak[i+1,:,:]
accum_counter = n -> Int(n*(n+1)/2)
get_cone = (i,j) -> accum_counter(i-3)+j;
get_ncone = n -> (i,j) -> accum_counter(i-3)+j-accum_counter(max(3+n,i)-(3+n));
get_nkcone = get_ncone(npastK)
N_K = get_ncone(npastK)(N,npastK)

# Robustness parameters
ϵ2 = 0.05; # Probability of failure
max_wind = 0.01;
ρ = max_wind*sqrt(2*log(1/ϵ2));

# Model
x_k_det = (k,α) -> Ai(k-1)*x0+sum(Ai(k-i-2)*B*α[:,i+1] for i=0:k-2)
model = Model(optimizer_with_attributes(
        () -> Gurobi.Optimizer(GRB_ENV), #"OutputFlag" => 0
    ));
@variable(model, K[1:N_K,1:2,1:2])   # Control law linear
@variable(model, α[1:2, 1:N-1])      # Control law constant
@variable(model, au[1:2, 1:N-1] >=0) # Absolute value of u(control)
@variable(model, z[1:N-1, 1:m] >=0, Bin) # Binary vars for obstacles
#@variable(model, minϵ[1:4] >=0)

# Variables for auxiliary second order cones:
@variable(model, t[1:N_K,1:m] >= 0) # (we need one for each halfplane)
@variable(model, s[1:N_K,1:2] >= 0); # (one for each control vector)
@variable(model, r[1:npastK,1:4] >= 0); # (one for each state variable)

@constraint(model, [k=1:N-1, i=1:Nobj],
   sum(z[k,j] for j in obj[i]) == length(obj[i])-1); # Obstacles


function get_K_independent_vars(k, S=I, size=2)
   if npastK+1<=k-2
           return ρ*mapslices(norm, sum(
                   (S*Ai(k-i-2)*Bw) for i=npastK+1:k-2), dims=2)[:]
   end
   return zeros(size)
end

# Obstacle free path
# Time step k=2
@constraint(model, (P*Q*x_k_det(2,α)
        +ρ*mapslices(norm, P*Q*Bw, dims=2)[:]
        ) .<= -q+M*z[1,:]);

# Time steps k=3 to N:
@constraint(model, [k=3:N], (P*Q*x_k_det(k,α)
        + get_K_independent_vars(k,P*Q, 2)
        +sum(ρ*t[get_nkcone(k,i),:] for i=1:min(k-2,npastK))
        +ρ*mapslices(norm, P*Q*Bw, dims=2)[:]
        ) .<= -q+M*z[k-1,:])



# Maximum control (2 inequalities per k)
@constraint(model, α[:,1] .<= au[:,1])
@constraint(model, -au[:,1] .<= α[:,1])
@constraint(model, [k=2:N-1], (α[:,k]
    +ρ*sum(s[get_nkcone(k,i),:] for i=1:min(npastK,k-2))
    ) .<= au[:,k])
@constraint(model, [k=2:N-1],  (-α[:,k]
    +ρ*sum(s[get_nkcone(k,i),:] for i=1:min(npastK,k-2))
    ) .<= au[:,k]);

# Goals (4 inequalities)
@constraint(model, (x_k_det(N, α)
        + get_K_independent_vars(N, I, 4)
        + sum(ρ * r[i,:] for i=1:min(npastK, N-2))
        + ρ * mapslices(norm, Bw, dims=2)[:]
        ) .<= xN+ϵ);
@constraint(model, (-x_k_det(N, α)
        + get_K_independent_vars(N, I, 4)
        + sum(ρ * r[i,:] for i=1:min(npastK, N-2))
        + ρ * mapslices(norm, Bw, dims=2)[:]
        ) .<= -xN+ϵ);

# # Second order constraints

# Obstacles
@constraint(model, [k=3:N,i=1:min(k-2,npastK),l=1:m], vcat(t[get_nkcone(k,i),l],
      (P*Q*Ai(k-i-2)*Bw)[l,:]+sum((P*Q*Ai(k-j-1)*B*K[get_nkcone(j,i),:,:])[l,:]
                for j=(i+1):(k-1))) in SecondOrderCone())
# Control
@constraint(model, [k=3:N-1,i=1:min(k-2,npastK),l=1:2],
        vcat(s[get_nkcone(k,i),l],
        sum(K[get_nkcone(j+1,i),l,:] for j=(i+1):(k-1))) in SecondOrderCone())

# Goals
@constraint(model, [i=1:min(npastK,N-2),l=1:4], vcat(r[i,l],
        (Ai(N-i-2)*Bw)[l,:]+sum(
                (Ai(N-j-1)*B*K[get_nkcone(j+1,i),:,:])[l,:] for j=(i+1):(N-1))
        ) in SecondOrderCone());
# Objective
@objective(model, Min, sum(au));

# Optimize
optimize!(model)

K_out = value.(K)
α_out = value.(α)
s_out = value.(s)
au_out = value.(au)

rectangle(w, h, x, y) = Shape(ones(4)*x + [0,w,w,0], ones(4)*y + [0,0,h,h])

iterations = 100
Γ = 1
xs2 = zeros(Float64, iterations, N)
ys2 = zeros(Float64, iterations, N)
faults = zeros(Int, 1)
for i = 1:iterations
        w = Γ * ρ * randn(Float64, (2, N - 1))
        y = zeros(4, N)
        y[:, 1] = x0
        for j = 2:N
                u = α_out[:, j - 1]
                #         u = maximum.([α_out[:, j - 1], -α_out[:, j - 1]])
                for k = 1:min(j-2,npastK)
                        #u += ρ*s_out[get_cone(j,k),:]
                        u += K_out[get_nkcone(j, k), :, :]*w[:, k]
                end
                y[:,j] = A*y[:, j-1] + Bw*w[:, j-1] + B*u
                for o in obj
                        o_size_temp = size(o)
                        o_size = o_size_temp[1]
                        rhs = zeros(Float64, o_size, 1)
                        for o2 = 1:o_size
                                rhs[o2,1] = q[o[o2]]
                        end
                        y_vec = zeros(Float64, 2, 1)
                        y_vec[1,1] = y[1,j]
                        y_vec[2,1] = y[3,j]
                        if sum(P[o,:]*y_vec.<= rhs) == o_size
                                faults[1] += 1
                        end
                end
        end
        xs2[i, :] = y[1,:]
        ys2[i, :] = y[3,:];
end
println(faults)

pl = plot(rectangle(0.035,0.03,0.015,0), opacity=.5, legend=:none)
plot!(xs2[:,:]', ys2[:,:]', linestyle=:dash, markershape=:circle, ms=2)
display(pl)

# Polygon constructor
ap = (p1,p2) -> [p1[2]-p2[2],p2[1]-p1[1]]
bp = (p1,p2) -> p2[1]*p1[2]-p1[1]*p2[2]
function makePq(points)
        s = size(points)
        s2 = zeros(Int, s[1])
        for i = 1:s[1]
                s3 = size(points[i])
                s2[i] = s3[1]
        end
        P = zeros(sum(s2), 2)
        q = zeros(sum(s2))
        idx = 1
        for i = 1:s[1]
                for (p1,p2) in zip(points[i],
                        vcat(points[i][2:end],points[i][1]))
                        println(p1,p2)
                        P[idx,:] = ap(p1, p2)
                        q[idx] = bp(p1, p2)
                        idx += 1
                end
        end
        return (P, q)
end

# polygons = [[(0.015,0),(0.015,0.03),(0.05,0.03),(0.05,0)]]
#
# pl = plot()
# for p in polygons
#         display(plot!(Shape(p), fillcolor = plot_color(:yellow, 0.3)))
# end
#
# Pp,qp = makePq(polygons)
