accum_counter = n -> Int(n*(n+1)/2)
get_cone = (i,j) -> accum_counter(i-3)+j;
get_ncone = n -> (i,j) -> accum_counter(i-3)+j-accum_counter(max(3+n,i)-(3+n));

function check_collisions(x, y, P, q, obj)
        iterations = size(x, 1)
        Np = size(x, 2)
        faults = zeros(iterations)
        for i = 1:iterations
                for j = 1:Np
                        y_vec = [x[i,j]; y[i,j]]
                        for o in obj
                                if sum(P[o,:]*y_vec .<= q[o]) == size(o,1)
                                        faults[i] += 1
                                end
                        end
                end
        end
        return faults
end

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
