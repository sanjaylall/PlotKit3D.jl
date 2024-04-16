
module Surface1





using PlotKit
using PlotKit3D
using LinearAlgebra
using Zygote


plotpath(x) = joinpath(ENV["HOME"], "plots/", x)




peaks(x,y) =  3*(1-x)^2*exp(-(x^2) - (y+1)^2) - 10*(x/5 - x^3 - y^5)*exp(-x^2-y^2) - 1/3*exp(-(x+1)^2 - y^2)

function main()
    zfn(x,y) =  peaks(x,y)
    dzfn(x,y) = Zygote.gradient(zfn, x,y)
    rt = surf(zfn, dzfn; xmin = -3, xmax = 3, ymin = -3, ymax = 3,
              axisoptions3_azimuth = -25, axisoptions3_elevation = 40)
    d = draw(rt)
    save(d, plotpath("surface1.pdf"))
end








end
