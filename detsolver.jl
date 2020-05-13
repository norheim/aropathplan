Q = zeros(2,4)
Q[1,1] = Q[2,3] = -1
bigM = 100; #Big M

function detpathplanner(problem)
    Ai, Bw, B, N, x_k_det, P, q, obj, xN, ϵ_in, ρ = problem
    PQ = P*Q
    Nobj = length(obj)
    m = length(q)

    model = Model(optimizer_with_attributes(
            () -> Gurobi.Optimizer(GRB_ENV), #"OutputFlag" => 0
        ));
    @variable(model, x[1:4, 1:N])
    @variable(model, u[1:2, 1:N-1])
    @variable(model, au[1:2, 1:N-1] >=0)
    @variable(model, z[1:N-1,1:m] >=0, Bin)
    @constraint(model, [k=2:N], PQ*x[:,k] .<= -q + bigM*z[k-1,:] )

    @constraint(model, [k=1:N-1, i=1:Nobj],
           sum(z[k,j] for j in obj[i]) <= length(obj[i])-1) # Obstacles

    @constraint(model, x[:,1] .== x0)
    @constraint(model, x[:,N] .== xN)
    @constraint(model, [k=1:N-1], x[:,k+1] .== Ai(1)*x[:,k]+B*u[:,k])
    @constraint(model, u .<= au)
    @constraint(model, -au .<= u)
    @objective(model, Min, sum(au));
    optimize!(model)

    return value.(x)[1,:],value.(x)[3,:], objective_value(model)
end
