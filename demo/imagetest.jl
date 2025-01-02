
module ImageTest

using PlotKit
using PlotKit3D
using LinearAlgebra
using Cairo

const pk3 = PlotKit3D

plotpath(x) = joinpath(ENV["HOME"], "plots/", x)


function getsurf2()
    oldzfun(x,y) =  3*(1-x)^2*exp(-(x^2) - (y+1)^2) - 10*(x/5 - x^3 - y^5)*exp(-x^2-y^2) - 1/3*exp(-(x+1)^2 - y^2)
    function olddzfun(x,y)
        dzdx = -(2/3)*exp(-2x-x*x-(1+y)^2)*(-exp(2y)*(1+x) + 9*exp(2x)*(1 - 2*x*x + x*x*x) + exp(1 + 2x + 2y)*(3 - 51x*x + 30*x^4 + 30*x*y^5))
        dzdy =  (1/3)*exp(-2x-x*x-(1+y)^2)*(2*exp(2y)*y - 18*exp(2x)*(-1+x)^2*(1+y) -6*exp(1+2x+2y)*y*(-2x + 10*x^3 + 5*y^3*(-5+2*y*y)))
        return dzdx, dzdy
    end
    return oldzfun, olddzfun
end


function rt1()
    f, df = getsurf2()
    rt = pk3.surf(f, df; xmin = -3, xmax = 3, ymin = -3, ymax = 3,
              axisoptions3_azimuth = -25,
              axisoptions3_elevation = atan(1/sqrt(2))*180/pi)
    d = draw(rt)
    save(d, plotpath("imagetest1.pdf"))
end


function rt2()
    f, df = getsurf2()
    rt = pk3.surf(f, df; xmin = -3, xmax = 3, ymin = -3, ymax = 3,
              axisoptions3_azimuth = -25,
              axisoptions3_elevation = atan(1/sqrt(2))*180/pi)
    d = draw(rt)
    save(d, plotpath("imagetest1.png"))
end

#
# Cairo ARGB32 uses premultiplied alpha!!!
#
# https://www.cairographics.org/documentation/cairomm/reference/classCairo_1_1Surface.html
#

fftoi(a::Float64) = round(UInt32, a * 255)
pack(r::Float64, g::Float64, b::Float64, a::Float64) = fftoi(a) << 24 + fftoi(a*r) << 16 + fftoi(a*g) << 8 + fftoi(a*b)

function rt3(; fname = "imagetest2.pdf")
    width = 800
    height = 600
    dw = Drawable(width, height; fname = plotpath(fname))
    rect(dw, Point(0,0), Point(800, 600); fillcolor = Color(1,0.8,0.9))
    line(dw, Point(10, 20), Point(800, 450); linestyle = LineStyle(Color(:blue), 5))
    line(dw, Point(10, 500), Point(600, 50); linestyle = LineStyle(Color(0,0.7,0.3), 5))


    tri_mask = [Point(50,50), Point(50,300), Point(400,175)]
    #line(dw, tri_mask; closed=true, linestyle = LineStyle(Color(:black), 5), fillcolor = Color(:red))

    iw = 400
    ih = 300
    # x, y are pixel coords
    # returns a UInt32 for each coord
    function colfn(x, y)
        r = x / iw
        if y / ih < 0.5
            g = 0.2
        else
            g = 0.8
        end
        b = 0.0
        a = 0.5
        return pack(r,g,b,a)
    end
    Z = UInt32[colfn(x,y) for x = 1:iw, y=1:ih]
    pik = Pik(Z, iw, ih)
    drawimage_to_mask(dw.ctx, pik, tri_mask, 1, 1; operator = Cairo.OPERATOR_OVER)
    close(dw)
    return pik
end




function main1()
    rt1()
    rt2()
    rt3()
    rt3(; fname = "imagetest2.png")
end



end
