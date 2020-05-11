using LinearAlgebra
# Create model for dynamics
function get_dynamic_matrices(T, N)
    dt = T/N
    A = zeros(4,4)+I
    A[1,2] = A[3,4] = dt
    B = zeros(4,2)
    B[2,1] = B[4,2] = dt
    #C = zeros(2,4)
    #C[1,2] = C[2,4] = 1
    #D = zeros(2,2)
    Bw = copy(B);

    Ak = zeros(N+1,size(A)...)
    Ak[1,:,:] = zeros(Int, size(A))+I
    for i=2:N+1
        Ak[i,:,:] = A*Ak[i-1,:,:]
    end
    Ai = i -> Ak[i+1,:,:]

    x_k_det = (x0, k, α) -> Ai(k-1)*x0+sum(Ai(k-i-2)*B*α[:,i+1] for i=0:k-2)
    x_k_stoc = (x0, k, α, K, w) -> (x_k_det(k,α)
        +sum((Ai(k-i-1) * Bw
            +sum(Ai(k-j-1) * B * K(j+1, i) for j=i+1:k-1))*w[:,i] for i=1:k-2)
        +Bw*w[:,k-1])

    function simulate_stoc(x0, α, K, w)
        x = zeros(4,N)
        x[:, 1] = x0
        for k = 2:N
            u = α[:, k - 1]
            for i = 1:k-2
                u += K(k, i)*w[:,i]
            end
            x[:, k] = A*x[:, k-1] + Bw*w[:,k-1] + B*u
        end
        return x
    end

    return x_k_det, simulate_stoc, Ai, B, Bw
end
