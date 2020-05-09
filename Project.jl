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
x0 = [0;0;0;0]
xN = [0.03;0;0.035;0]
ϵ = [0.01;
        0.003;
        0.01;
        0.003]; # how close to the goal we want to get

# Obstacles
# P = [-1 0; 0 1]
# q = [-0.03; 0.03] # this should be feasible
P = [-1 0; 0 1]
q = [-0.02; 0.03]
m = length(q)
obj = [[1,2]] # Each entry corresponds to the lines defining an object
Nobj = length(obj)
M = 100; #Big M

# Helper variables
Q = zeros(2,4)
Q[1,1] = Q[2,3] = -1
N_K = Int((N-2)*(N-1)/2) # number of K's
Ak = zeros(N+1,4,4)
Ak[1,:,:] = zeros(Int, size(A))+I
for i=2:N+1
    Ak[i,:,:] = A*Ak[i-1,:,:]
end
Ai = i -> Ak[i+1,:,:]
get_cone = (i,j) -> Int((i-2)*(i-3)/2)+j

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
# Variables for auxiliary second order cones:
@variable(model, t[1:N_K,1:m] >= 0) # (we need one for each halfplane)
@variable(model, s[1:N_K,1:2] >= 0); # (one for each control vector)
@variable(model, r[1:N_K,1:4] >= 0); # (one for each state variable)

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

# Maximum control
@constraint(model, α[:,1] .<= au[:,1])
@constraint(model, -au[:,1] .<= α[:,1])
@constraint(model, [k=2:N-1], (α[:,k]
    +ρ*sum(s[get_cone(k,i),:] for i=1:k-2) 
    ) .<= au[:,k])
@constraint(model, [k=2:N-1],  (-α[:,k]
    +ρ*sum(s[get_cone(k,i),:] for i=1:k-2) 
    ) .<= au[:,k]);

# Goals
@constraint(model, (x_k_det(N,α)
        +sum(ρ*r[get_cone(N,i),:] for i=1:N-2)
        +ρ*mapslices(norm, Bw, dims=2)[:]
        ) .<= xN+ϵ);
@constraint(model, (-x_k_det(N,α)
        +sum(ρ*r[get_cone(N,i),:] for i=1:N-2)
        +ρ*mapslices(norm, Bw, dims=2)[:] 
           ) .<= -xN+ϵ);

# # Second order constraints

# Obstacles
@constraint(model, [k=3:N,i=1:(k-2),l=1:m], vcat(t[get_cone(k,i),l], (P*Q*Ai(k-i-2)*Bw)[l,:]+sum(
            (P*Q*Ai(k-j-1)*B*K[get_cone(j,i),:,:])[l,:] for j=(i+1):(k-1))) in SecondOrderCone());
# Control
@constraint(model, [k=3:N-1,i=1:(k-2),l=1:2], vcat(s[get_cone(k,i),l], sum(
            K[get_cone(j,i),l,:] for j=(i+1):(k-1))) in SecondOrderCone());
# Goals
@constraint(model, [k=3:N,i=1:(k-2),l=1:4], vcat(r[get_cone(k,i),l], (Ai(k-i-2)*Bw)[l,:]+sum(
            (Ai(k-j-1)*B*K[get_cone(j,i),:,:])[l,:] for j=(i+1):(k-1))) in SecondOrderCone());
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
for i = 1:iterations
        w = Γ * ρ * randn(Float64, (2, N - 1))
        y = zeros(4, N)
        y[:, 1] = x0
        for j = 2:N
                #u = α_out[:, j - 1]
                #         u = maximum.([α_out[:, j - 1], -α_out[:, j - 1]])
        #         for k = 1:j-2
        #                 u += ρ*s_out[get_cone(j,k),:] 
        #                 #u += K_out[get_cone(j, k), :, :]*w[:, k]
        #         end
                y[:,j] = A*y[:, j-1] + Bw*w[:, j-1] + B*au_out[:,j-1]
        end
        xs2[i, :] = y[1,:]
        ys2[i, :] = y[3,:];
end

plot(rectangle(0.035,0.03,0.015,0), opacity=.5, legend=:none)
plot!(xs2[:,:]', ys2[:,:]', linestyle=:dash, markershape=:circle, ms=2)
