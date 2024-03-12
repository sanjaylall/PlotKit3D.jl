
module Demo3D

using PlotKit

plotpath(x) = joinpath(ENV["HOME"], "plots/", x)


##############################################################################
# for raytracing

function getsurf1()
    function plh(x)
        alpha = 1
        if abs(x)<alpha
            return x*x
        end
        return alpha*alpha*(1-2*log(alpha) + 2*log(abs(x)) )
    end
    
    function dplh(x)
        alpha = 1
        if abs(x) < alpha
            return 2*x
        end
        return 2*alpha*alpha/x
    end

    z_plh(x,y) = 1/3*( plh(x+3) + plh(2*y+6) + plh(x+y-1))

    function dz_plh(x,y)
        dzdx = 1/3*(dplh(x+3) + dplh(x+y-1))
        dzdy = 1/3*(2*dplh(2*y+6) + dplh(x+y-1))
        return dzdx, dzdy
    end
    return z_plh, dz_plh
end

##############################################################################
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
    save(d, plotpath("mesh1.pdf"))
end    


function main2()
    # 3d peaks plot, ala Matlab
    x = range(-3, 3, length = 49)
    y = range(-3, 3, length = 49)
    zf(x,y) =  3*(1-x)^2*exp(-(x^2) - (y+1)^2) - 10*(x/5 - x^3 - y^5)*exp(-x^2-y^2) - 1/3*exp(-(x+1)^2 - y^2)
    d = mesh_height_fn(x, y, zf)
    save(d, plotpath("mesh2.pdf"))
end



function main3()
    f, df = getsurf1()
    d=surf(f, df; xmin = -5, xmax = 5, ymin = -5, ymax = 5,
           axisoptions3_azimuth = -25, axisoptions3_elevation = 40)
    save(d, plotpath("raytrace3.pdf"))
end



end

