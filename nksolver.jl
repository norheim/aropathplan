using JuMP, Gurobi, LinearAlgebra
GRB_ENV = Gurobi.Env();
include("utils.jl")

Q = zeros(2,4) # Matrix to recover -x and -y coordinates from state
Q[1,1] = Q[2,3] = -1
bigM = 100; # Big M (do not make too big, otherwise need to adjust
            # tolerance of sovler)

function pathplanner(problem, npastK=0, feasability=false)
        Ai, Bw, B, N, x_k_det, P, q, obj, xN, max_au, ϵ_in, ρ = problem
        PQ = P*Q
        Nobj = length(obj)
        m = length(q)
        if npastK == 0
                npastK = N-2
        end
        get_nkcone = get_ncone(npastK)
        N_K = get_ncone(npastK)(N, npastK)

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
        @variable(model, r[1:npastK, 1:4] >= 0) # (one for each state variable)

        @constraint(model, [k=1:N-1, i=1:Nobj],
           sum(z[k,j] for j in obj[i]) <= length(obj[i])-1) # Obstacles

        function get_K_independent_vars(k, S=I, size=2)
                if npastK+1<=k-2
                      return ρ*mapslices(norm, sum(
                              (S*Ai(k-i-2)*Bw) for i=npastK+1:k-2), dims=2)[:]
                end
                return zeros(size)
        end

        # Obstacle free path
        # Time step k=2
        @constraint(model, (PQ*x_k_det(2,α)
                + ρ*mapslices(norm, PQ*Bw, dims=2)[:]
                ) .<= -q+bigM*z[1,:])
        # Time steps k=3 to N:
        @constraint(model, [k=3:N], (PQ*x_k_det(k,α)
                + get_K_independent_vars(k, PQ, m)
                + sum(ρ*t[get_nkcone(k, i),:] for i=1:min(k-2,npastK))
                + ρ*mapslices(norm, PQ*Bw, dims=2)[:]
                ) .<= -q+bigM*z[k-1,:])

        # Maximum control (2 inequalities per k)
        @constraint(model, au .<= max_au)
        @constraint(model, α[:, 1] .<= au[:, 1])
        @constraint(model, -au[:, 1] .<= α[:, 1])
        @constraint(model, [k=2:N-1], (α[:, k]
            +ρ*sum(s[get_nkcone(k, i),:] for i=1:min(npastK,k-2))
            ) .<= au[:,k])
        @constraint(model, [k=2:N-1],  (-α[:,k]
            +ρ*sum(s[get_nkcone(k, i),:] for i=1:min(npastK,k-2))
            ) .<= au[:,k])

        # Goals (4 inequalities)
        @constraint(model, (x_k_det(N, α)
                + get_K_independent_vars(N, I, 4)
                + sum(ρ*r[i,:] for i=1:min(npastK, N-2))
                + ρ*mapslices(norm, Bw, dims=2)[:]
                ) .<= xN+ϵ)
        @constraint(model, (-x_k_det(N, α)
                + get_K_independent_vars(N, I, 4)
                + sum(ρ*r[i,:] for i=1:min(npastK, N-2))
                + ρ*mapslices(norm, Bw, dims=2)[:]
                   ) .<= -xN+ϵ)

        # Obstacles SOC
        @constraint(model, [k=3:N, i=1:min(k-2, npastK), l=1:m], vcat(
                t[get_nkcone(k, i), l],
                (PQ * Ai(k-i-1) * Bw)[l,:]+sum(
                (PQ * Ai(k-j-1) * B * K[get_nkcone(j + 1, i),:,:])[l,:]
                for j=(i+1):(k-1))) in SecondOrderCone())

        # Control SOC
        @constraint(model, [k=3:N-1, i=1:min(k-2, npastK), l=1:2], vcat(
                s[get_nkcone(k, i), l], sum(
                K[get_nkcone(j + 1, i), l, :] for j=(i+1):(k-1)))
                in SecondOrderCone())

        # Goals SOC
        @constraint(model, [i=1:min(N-2, npastK), l=1:4], vcat(r[i, l],
                (Ai(N-i-1)*Bw)[l,:]+sum(
                (Ai(N-j-1)*B*K[get_nkcone(j + 1, i),:,:])[l,:] for j=(i+1):(N-1)))
                in SecondOrderCone())

        # Time-out to return result
        set_time_limit_sec(model, 800)

        # Objective
        if feasability
                @objective(model, Min, sum(minϵ))
                optimize!(model)
                return value.(minϵ)
        else
                @objective(model, Min, sum(au))
                optimize!(model)
                K_val = value.(K)
                K_out = (i,j) -> (j <= min(i-2, npastK) ?
                        K_val[get_nkcone(i, j), :, :] :
                        zeros(size(K_val)[[2, 3]]))
                return K_out, value.(α), objective_value(model)
        end
end
