using Plots
rectangle(w, h, x, y) = Shape(ones(4)*x + [0,w,w,0], ones(4)*y + [0,0,h,h])

function plot_goal(xN, ϵ)
    plot!(rectangle(2*ϵ[1], 2*ϵ[3], xN[1]-ϵ[1], xN[3]-ϵ[3]),
        opacity=.5, legend=:none)
end

function plot_polygons(polygons)
    for p in polygons
        plot!(Shape(p), fillcolor = plot_color(:yellow, 0.3))
    end
end
