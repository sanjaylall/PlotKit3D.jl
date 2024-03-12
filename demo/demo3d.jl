
module Demo3D

using PlotKit

plotpath(x) = joinpath(ENV["HOME"], "plots/", x)

function main1()
    x = range(-6, 6, length = 101)
    y = range(-6, 6, length = 101)
    qf(x,y) = sqrt(x*x + y*y) - 3
    zmin = -3
    zmax = maximum([qf(a,b) for a in x for b in y])
    function cfun(x,y)
        r=(qf(x,y)-zmin)/zmax
        return  Color(r,  0.5,  0.9-r/2)
    end
    d = mesh_height_color_fn(x, y, qf, cfun; xlineindices=3:4:100,
         ylineindices=3:4:100, azimuth=-70, elevation=60)
    save(d, plotpath("mesh2_x.pdf"))
end    


end

