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
npastK = N
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
N_K = get_ncone(npastK)(N,N-2)

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
@variable(model, r[1:min(npastK,N-2),1:4] >= 0); # (one for each state variable)

@constraint(model, [k=1:N-1, i=1:Nobj],
   sum(z[k,j] for j in obj[i]) == length(obj[i])-1); # Obstacles

# Obstacle free path
# Time step k=2
@constraint(model, (P*Q*x_k_det(2,α)
        +ρ*mapslices(norm, P*Q*Bw, dims=2)[:]
        ) .<= -q+M*z[1,:]);
# Time steps k=3 to N:
@constraint(model, [k=3:N], (P*Q*x_k_det(k,α)
        +sum(ρ*t[get_cone(k,i),:] for i=1:k-2)
        +ρ*mapslices(norm, P*Q*Bw, dims=2)[:]
        ) .<= -q+M*z[k-1,:]);

# Limited K's
# Obstacles
# @constraint(model, [k=3:N], (P*Q*x_k_det(k,α)
#        +ρ*mapslices(norm, sum((P*Q*Ai(k-i-2)*Bw) for i=npastK+1:k-2), dims=2)[:]
#         +sum(ρ*t[get_nkcone(k,i),:] for i=1:min(k-2,npastK))
#         +ρ*mapslices(norm, P*Q*Bw, dims=2)[:]
#         ) .<= -q+M*z[k-1,:])
#
# #Control
# @constraint(model, [k=2:N-1], (α[:,k]
#     +ρ*sum(s[get_nkcone(k,i),:] for i=1:max(3,k-2))
#     ) .<= au[:,k])
#
# #Goal
# @constraint(model, (x_k_det(N,α)
#        +ρ*mapslices(norm, sum(Ai(k-i-2)*Bw for i=npastK+1:k-2), dims=2)[:]
#         +sum(ρ*r[i,:] for i=1:min(npastK,N-2))
#         +ρ*mapslices(norm, Bw, dims=2)[:]
#         ) .<= xN+ϵ);
#
# @constraint(model, [k=3:N,i=1:min(k-2,npastK),l=1:m], vcat(t[get_nkcone(k,i),l],
#       (P*Q*Ai(k-i-2)*Bw)[l,:]+sum((P*Q*Ai(k-j-1)*B*K[get_nkcone(j,i),:,:])[l,:]
#                 for j=(i+1):(k-1))) in SecondOrderCone())
#
# @constraint(model, [k=3:N-1,i=1:min(k-2,npastK),l=1:2],
#         vcat(s[get_nkcone(k,i),l],
#         sum(K[get_nkcone(j+1,i),l,:] for j=(i+1):(k-1))) in SecondOrderCone())
#
# @constraint(model, [i=1:min(npastK,N-2),l=1:4], vcat(r[i,l],
#         (Ai(N-i-2)*Bw)[l,:]+sum(
#                 (Ai(N-j-1)*B*K[get_nkcone(j+1,i),:,:])[l,:] for j=(i+1):(N-1))
#         ) in SecondOrderCone());

# Maximum control (2 inequalities per k)
@constraint(model, α[:,1] .<= au[:,1])
@constraint(model, -au[:,1] .<= α[:,1])
@constraint(model, [k=2:N-1], (α[:,k]
    +ρ*sum(s[get_cone(k,i),:] for i=1:k-2)
    ) .<= au[:,k])
@constraint(model, [k=2:N-1],  (-α[:,k]
    +ρ*sum(s[get_cone(k,i),:] for i=1:k-2)
    ) .<= au[:,k]);

# Goals (4 inequalities)
@constraint(model, (x_k_det(N,α)
        +sum(ρ*r[i,:] for i=1:N-2)
        +ρ*mapslices(norm, Bw, dims=2)[:]
        ) .<= xN+ϵ);
@constraint(model, (-x_k_det(N,α)
        +sum(ρ*r[i,:] for i=1:N-2)
        +ρ*mapslices(norm, Bw, dims=2)[:]
           ) .<= -xN+ϵ);

# # Second order constraints

# k=N
# i=N-2
# [get_cone(j+1,i) for j=(i+1):(k-1)]
# get_cone(3,1)

# Obstacles
@constraint(model, [k=3:N,i=1:k-2,l=1:m], vcat(t[get_cone(k,i),l], (P*Q*Ai(k-i-2)*Bw)[l,:]+sum(
            (P*Q*Ai(k-j-1)*B*K[get_cone(j+1,i),:,:])[l,:] for j=(i+1):(k-1))) in SecondOrderCone());
# Control
@constraint(model, [k=3:N-1,i=1:k-2,l=1:2], vcat(s[get_cone(k,i),l], sum(
            K[get_cone(j+1,i),l,:] for j=(i+1):(k-1))) in SecondOrderCone());
# Goals
@constraint(model, [i=1:N-2,l=1:4], vcat(r[i,l], (Ai(N-i-2)*Bw)[l,:]+sum(
            (Ai(N-j-1)*B*K[get_cone(j+1,i),:,:])[l,:] for j=(i+1):(N-1))) in SecondOrderCone());
# Objective
@objective(model, Min, sum(au));

# Optimize
optimize!(model)

K_out = value.(K)
α_out = value.(α)
s_out = value.(s)
au_out = value.(au)

rectangle(w, h, x, y) = Shape(ones(4)*x + [0,w,w,0], ones(4)*y + [0,0,h,h])

iterations = 10
Γ = .125
xs2 = zeros(Float64, iterations, N)
ys2 = zeros(Float64, iterations, N)
faults = zeros(Int, 1)
for i = 1:iterations
        w = Γ * max_wind * randn(Float64, (2, N - 1))
        y = zeros(4, N)
        y[:, 1] = x0
        for j = 2:N
                u = α_out[:, j - 1]
                #         u = maximum.([α_out[:, j - 1], -α_out[:, j - 1]])
                for k = 1:j-2
                        #u += ρ*s_out[get_cone(j,k),:]
                        u += K_out[get_cone(j, k), :, :]*w[:, k]
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

# polygons = [[(-0.05,0),(-1,0),(-1,3.5),(-.05,3.5)],[(0,-1),(0,-0.01),(4,-0.01),(4,-1)],[(0,2),(0,4),(3,2)],[(4,3),(1,5),(4,6)],[(4,0),(4,7),(5,7),(5,0)]];
polygons = [[(-0.0125,0),(0.0375,4),(2,2)],[(0,5),(2,8),(2,4)]]

# Create model
dt = 10/30;
T = 10;
N = Int(T/dt);
x0 = [1;0;-1;0];
xN = [1;0;9;0];
ϵ = [0.5;
     10;
     2;
     10];
    # how close to the goal we want to get

pl = plot()
pl = plot(rectangle(2*ϵ[1], 2*ϵ[3], xN[1]-ϵ[1], xN[3]-ϵ[3]), opacity=.5, legend=:none)
for p in polygons
        display(plot!(Shape(p), fillcolor = plot_color(:yellow, 0.3)))
end

P,q = makePq(polygons);

# Starting here, we will be doing an example and vary different parts

obj = [[1,2,3],[4,5,6]];

# Create model for dynamics
A = zeros(4,4)+I;
A[1,2] = A[3,4] = dt;
B = zeros(4,2);
B[2,1] = B[4,2] = dt;
C = zeros(2,4);
C[1,2] = C[2,4] = 1;
D = zeros(2,2);
Bw = copy(B);

# display(pl)

# Obstacles
m = length(q);
Nobj = length(obj);
M = 100; #Big M

# Helper variables
Q = zeros(2,4);
Q[1,1] = Q[2,3] = -1;
N_K = Int((N-2)*(N-1)/2); # number of K's
#N_s = Int((N-3)*(N-2)/2)
Ak = zeros(N+1,4,4);
Ak[1,:,:] = zeros(Int, size(A))+I;
for i=2:N+1
    Ak[i,:,:] = A*Ak[i-1,:,:];
end
Ai = i -> Ak[i+1,:,:];
get_cone = (i,j) -> Int((i-2)*(i-3)/2)+j;

# Robustness parameters
ϵ2 = 0.05; # Probability of failure
# max_wind = 0.03;
max_wind = 1;
# max_au = 10;
ρ = max_wind*sqrt(2*log(1/ϵ2));
ρ
# Model
x_k_det = (k,α) -> Ai(k-1)*x0+sum(Ai(k-i-2)*B*α[:,i+1] for i=0:k-2);
model = Model(optimizer_with_attributes(
        () -> Gurobi.Optimizer(GRB_ENV), #"OutputFlag" => 0
    ));
@variable(model, K[1:N_K,1:2,1:2]);   # Control law linear
@variable(model, α[1:2, 1:N-1]);      # Control law constant
@variable(model, au[1:2, 1:N-1] >=0); # Absolute value of u(control)
@variable(model, z[1:N-1, 1:m] >=0, Bin); # Binary vars for obstacles
#@variable(model, minϵ[1:4] >=0)

# Variables for auxiliary second order cones:
@variable(model, t[1:N_K,1:m] >= 0); # (we need one for each halfplane)
@variable(model, s[1:N_K,1:2] >= 0); # (one for each control vector)
@variable(model, r[1:N-2,1:4] >= 0); # (one for each state variable)

@constraint(model, [k=1:N-1, i=1:Nobj],
   sum(z[k,j] for j in obj[i]) == length(obj[i])-1); # Obstacles

# Obstacle free path
# Time step k=2
@constraint(model, (P*Q*x_k_det(2,α)
        +ρ*mapslices(norm, P*Q*Bw, dims=2)[:]
        ) .<= -q+M*z[1,:]);
# Time steps k=3 to N:
@constraint(model, [k=3:N], (P*Q*x_k_det(k,α)
        +sum(ρ*t[get_cone(k,i),:] for i=1:k-2)
        +ρ*mapslices(norm, P*Q*Bw, dims=2)[:]
        ) .<= -q+M*z[k-1,:]);

# Maximum control (2 inequalities per k)
@constraint(model, α[:,1] .<= au[:,1]);
@constraint(model, -au[:,1] .<= α[:,1]);
@constraint(model, [k=2:N-1], (α[:,k]
    +ρ*sum(s[get_cone(k,i),:] for i=1:k-2)
    ) .<= au[:,k]);
@constraint(model, [k=2:N-1],  (-α[:,k]
    +ρ*sum(s[get_cone(k,i),:] for i=1:k-2)
    ) .<= au[:,k]);

# Goals (4 inequalities)
@constraint(model, (x_k_det(N,α)
        +sum(ρ*r[i,:] for i=1:N-2)
        +ρ*mapslices(norm, Bw, dims=2)[:]
        ) .<= xN+ϵ);
@constraint(model, (-x_k_det(N,α)
        +sum(ρ*r[i,:] for i=1:N-2)
        +ρ*mapslices(norm, Bw, dims=2)[:]
           ) .<= -xN+ϵ);

# # Second order constraints
# @constraint(model, [k=1:N-1], au[1,k] ^ 2 + au[2,k] ^ 2 <= max_au ^ 2);

# k=N
# i=N-2
# [get_cone(j+1,i) for j=(i+1):(k-1)]
# get_cone(3,1)

# Obstacles
@constraint(model, [k=3:N,i=1:k-2,l=1:m], vcat(t[get_cone(k,i),l], (P*Q*Ai(k-i-2)*Bw)[l,:]+sum(
            (P*Q*Ai(k-j-1)*B*K[get_cone(j+1,i),:,:])[l,:] for j=(i+1):(k-1))) in SecondOrderCone());
# Control
@constraint(model, [k=3:N-1,i=1:k-2,l=1:2], vcat(s[get_cone(k,i),l], sum(
            K[get_cone(j+1,i),l,:] for j=(i+1):(k-1))) in SecondOrderCone());
# Goals
@constraint(model, [i=1:N-2,l=1:4], vcat(r[i,l], (Ai(N-i-2)*Bw)[l,:]+sum(
            (Ai(N-j-1)*B*K[get_cone(j+1,i),:,:])[l,:] for j=(i+1):(N-1))) in SecondOrderCone());

# Objective
@objective(model, Min, sum(au));

# Optimize
set_time_limit_sec(model, 600)
optimize!(model)

α_out = value.(α)

x = zeros(Float64,4,N);

x[:,1] = x0;
for k=2:N
        x[:,k] = A * x[:,k-1] + B * α_out[:,k-1];
end

α_out

xs = zeros(Float64, N)
ys = zeros(Float64, N)
for k=1:N
        xs[k] = x[1,k];
        ys[k] = x[3,k];
end
plot!(xs, ys, linestyle=:dash, markershape=:circle, ms=2)

# print(x)
