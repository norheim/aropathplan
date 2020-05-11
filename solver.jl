using JuMP, Gurobi, LinearAlgebra
GRB_ENV = Gurobi.Env();
include("utils.jl")

Q = zeros(2,4)
Q[1,1] = Q[2,3] = -1
bigM = 100; #Big M

function pathplanner(problem, feasability=false)
        #feasability=true
        #ϵ_in = ϵ
        Ai, Bw, B, N, x_k_det, P, q, obj, xN, ϵ_in, ρ = problem
        PQ = P*Q
        Nobj = length(obj)
        m = length(q)
        npastK = N
        get_nkcone = get_ncone(npastK)
        N_K = get_ncone(npastK)(N, N-2)

        model = Model(optimizer_with_attributes(
                () -> Gurobi.Optimizer(GRB_ENV), #"OutputFlag" => 0
            ));

        @variable(model, K[1:N_K, 1:2, 1:2])   # Control law linear
        @variable(model, α[1:2, 1:N-1])      # Control law constant
        @variable(model, au[1:2, 1:N-1] >=0) # Absolute value of u(control)
        @variable(model, z[1:N-1, 1:m] >=0, Bin) # Binary vars for obstacles
        if feasability
                @variable(model, minϵ[1:4] >=0)
                ϵ = minϵ
        else
                ϵ = ϵ_in
        end

        # Variables for auxiliary second order cones:
        @variable(model, t[1:N_K, 1:m] >= 0) # (we need one for each halfplane)
        @variable(model, s[1:N_K, 1:2] >= 0); # (one for each control vector)
        @variable(model, r[1:min(npastK,N-2), 1:4] >= 0) # (one for each state variable)

        @constraint(model, [k=1:N-1, i=1:Nobj],
           sum(z[k,j] for j in obj[i]) == length(obj[i])-1) # Obstacles

        # Obstacle free path
        # Time step k=2
        @constraint(model, (PQ*x_k_det(2,α)
                +ρ*mapslices(norm, PQ*Bw, dims=2)[:]
                ) .<= -q+bigM*z[1,:])
        # Time steps k=3 to N:
        @constraint(model, [k=3:N], (PQ*x_k_det(k,α)
                +sum(ρ*t[get_cone(k, i),:] for i=1:k-2)
                +ρ*mapslices(norm, PQ*Bw, dims=2)[:]
                ) .<= -q+bigM*z[k-1,:])

        # Maximum control (2 inequalities per k)
        @constraint(model, α[:, 1] .<= au[:, 1])
        @constraint(model, -au[:, 1] .<= α[:, 1])
        @constraint(model, [k=2:N-1], (α[:, k]
            +ρ*sum(s[get_cone(k, i),:] for i=1:k-2)
            ) .<= au[:,k])
        @constraint(model, [k=2:N-1],  (-α[:,k]
            +ρ*sum(s[get_cone(k, i),:] for i=1:k-2)
            ) .<= au[:,k])

        # Goals (4 inequalities)
        @constraint(model, (x_k_det(N, α)
                +sum(ρ*r[i,:] for i=1:N-2)
                +ρ*mapslices(norm, Bw, dims=2)[:]
                ) .<= xN+ϵ)
        @constraint(model, (-x_k_det(N, α)
                +sum(ρ*r[i,:] for i=1:N-2)
                +ρ*mapslices(norm, Bw, dims=2)[:]
                   ) .<= -xN+ϵ)

        # Obstacles SOC
        @constraint(model, [k=3:N, i=1:k-2, l=1:m], vcat(t[get_cone(k, i), l],
                (PQ * Ai(k-i-1) * Bw)[l,:]+sum(
                (PQ * Ai(k-j-1) * B * K[get_cone(j + 1, i),:,:])[l,:]
                for j=(i+1):(k-1))) in SecondOrderCone())

        # Control SOC
        @constraint(model, [k=3:N-1, i=1:k-2, l=1:2], vcat(
                s[get_cone(k, i), l], sum(
                K[get_cone(j + 1, i), l, :] for j=(i+1):(k-1)))
                in SecondOrderCone())

        # Goals SOC
        @constraint(model, [i=1:N-2, l=1:4], vcat(r[i, l],
                (Ai(N-i-1)*Bw)[l,:]+sum(
                (Ai(N-j-1)*B*K[get_cone(j + 1, i),:,:])[l,:] for j=(i+1):(N-1)))
                in SecondOrderCone())

        # Objective
        if feasability
                @objective(model, Min, sum(minϵ))
                optimize!(model)
                return value.(minϵ)
        else
                @objective(model, Min, sum(au))
                optimize!(model)
                K_val = value.(K)
                K_out = (i,j) -> K_val[get_cone(i, j), :, :]
                return K_out, value.(α)
        end
end
