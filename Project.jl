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
ϵ = 0.0001 # how close to the goal we want to get

# Obstacles
# P = [-1 0; 0 1]
# q = [-0.03; 0.03] # this should be feasible
P = [-1 0; 0 1]
q = [-0.02; 0.03]
m = length(q)
obj = [[1,2]] # Each entry corresponds to the lines defining an object
Nobj = length(obj)
M = 1000; #Big M

# Helper variables
Q = zeros(2,4)
Q[1,2] = Q[2,4] = 1
N_K = Int((N-1)*N/2) # number of K's
Ak = zeros(N,4,4)
Ak[1,:,:] = zeros(Int, size(A))+I
for i=2:N
    Ak[i,:,:] = A*Ak[i-1,:,:]
end
Ai = i -> Ak[i,:,:]
get_cone = (i,j) -> Int((i-2)*(i-3)/2)+j;

# Robustness parameters
ρ = 0.001;

# Model
model = Model(optimizer_with_attributes(
        () -> Gurobi.Optimizer(GRB_ENV), #"OutputFlag" => 0
    ));
@variable(model, K[1:N_K,1:2,1:2])   # Control law linear
@variable(model, α[1:2, 1:N])      # Control law constant
@variable(model, au[1:2, 1:N] >=0) # Absolute value of u(control)
@variable(model, z[1:N, 1:m] >=0, Bin) # Binary vars for obstacles
# Variables for auxiliary second order cones:
@variable(model, t[1:N_K,1:m] >= 0) # (we need one for each halfplane)
@variable(model, s[1:N_K,1:2] >= 0); # (one for each control vector)
@variable(model, r[1:N_K,1:4] >= 0); # (one for each state variable)

@constraint(model, cobj[k=2:N, i=1:Nobj],
   sum(z[k,j] for j in obj[i]) == length(obj[i])-1); # Obstacles

# Obstacle free path
# Time step k=2
@constraint(model, (P*Q*Ai(2)*x0+P*Q*B*α[:,2]
        +ρ*mapslices(norm, P*Q*Bw, dims=2)[:]) .<= q+M*z[2,:]);
# All time steps after:
@constraint(model, [k=3:N], (P*Q*Ai(k)*x0+sum(P*Q*Ai(k-i)*B*α[:,i] for i=1:k-1)
        +sum(ρ*t[get_cone(k,i),:] for i=1:k-2)
        +ρ*mapslices(norm, P*Q*Bw, dims=2)[:]) .<= q+M*z[k,:]);

# Maximum control
@constraint(model, α[:,2] .<= au[:,2])
@constraint(model, -au[:,2] .<= α[:,2])
@constraint(model, [k=3:N-1], α[:,k]+ρ*sum(s[get_cone(k,i),:] for i=1:k-2) .<= au[:,k])
@constraint(model, [k=3:N-1],  -α[:,k]+ρ*sum(s[get_cone(k,i),:] for i=1:k-2) .<= au[:,k]);

# Goals
@constraint(model, (Ai(N)*x0+sum(Ai(N-i)*B*α[:,i] for i=1:N-1)
        +sum(ρ*r[get_cone(N,i),:] for i=1:N-2)
        +ρ*mapslices(norm, Bw, dims=2)[:]) .<= xN+ϵ);
@constraint(model, (-Ai(N)*x0+sum(Ai(N-i)*B*-α[:,i] for i=1:N-1)
        +sum(ρ*r[get_cone(N,i),:] for i=1:N-2)
        +ρ*mapslices(norm, Bw, dims=2)[:]) .<= -xN+ϵ);

# # Second order constraints

# Obstacles
@constraint(model, obs[k=3:N,i=1:(k-2),l=1:m], vcat(t[get_cone(k,i),l], (P*Q*Ai(k-i)*Bw)[l,:]+sum(
            (P*Q*Ai(k-j)*B*K[get_cone(j,i),:,:])[l,:] for j=(i+1):(k-1))) in SecondOrderCone());

# Control
@constraint(model, [k=3:N-1,i=1:(k-2),l=1:2], vcat(s[get_cone(k,i),l], sum(
            K[get_cone(j,i),l,:] for j=(i+1):(k-1))) in SecondOrderCone());

# Goals
@constraint(model, [k=3:N,i=1:(k-2),l=1:4], vcat(r[get_cone(k,i),l], (Ai(k-i)*Bw)[l,:]+sum(
            (Ai(k-j)*B*K[get_cone(j,i),:,:])[l,:] for j=(i+1):(k-1))) in SecondOrderCone());

# Objective
@objective(model, Min, sum(au));

# Optimize
optimize!(model)

function getx(α, r)
    x = zeros(4,N)
    x[:,1] = x0
    x[:,2] = Ai(2)*x0+ρ*mapslices(norm, Bw, dims=2)[:]
    for k=3:N
         x[:,k] = (Ai(k)*x0+sum(Ai(k-i)*B*α[:,i] for i=1:k-1)
                    +sum(ρ*r[get_cone(k,i),:] for i=1:k-2)
                    +ρ*mapslices(norm, Bw, dims=2)[:])
    end
    return x
end

x = getx(value.(α), value.(r));
K_out = value.(K)
α_out = value.(α)
value.(z)

xs = x[1,:]
ys = x[3,:];

rectangle(w, h, x, y) = Shape(ones(4)*x + [0,w,w,0], ones(4)*y + [0,0,h,h])

plot(xs, ys, seriestype = :scatter)

iterations = 20
Γ = 2
xs2 = zeros(Float64, iterations, N)
ys2 = zeros(Float64, iterations, N)
for i = 1:iterations
        w = Γ * ρ * randn(Float64, (2, N - 1))
        y = zeros(4, N)
        y[:, 1] = x0
        for j = 2:N
                u = α_out[:, j - 1]
                for k = 1:j-1
                        u += K_out[get_cone(j + 1, k), :, :] * w[:, k]
                end
                y[:, j] = A * y[:, j - 1] + Bw * w[:, j - 1] + B * u
        end
        xs2[i, :] = y[1,:]
        ys2[i, :] = y[3,:];
end

xs2t = transpose(xs2) # Each row is an iteration. Column is time stamp
ys2t = transpose(ys2)
scatter(xs2t, ys2t) # This is not actually how it works...
