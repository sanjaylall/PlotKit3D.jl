
module Demo3D

using PlotKit
using PlotKit3D
using LinearAlgebra

const pk3 = PlotKit3D

plotpath(x) = joinpath(ENV["HOME"], "plots/", x)

eye(n) = Matrix{Float64}(I(n))
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

function getsurf2()
    oldzfun(x,y) =  3*(1-x)^2*exp(-(x^2) - (y+1)^2) - 10*(x/5 - x^3 - y^5)*exp(-x^2-y^2) - 1/3*exp(-(x+1)^2 - y^2)
    function olddzfun(x,y)
        dzdx = -(2/3)*exp(-2x-x*x-(1+y)^2)*(-exp(2y)*(1+x) + 9*exp(2x)*(1 - 2*x*x + x*x*x) + exp(1 + 2x + 2y)*(3 - 51x*x + 30*x^4 + 30*x*y^5))
        dzdy =  (1/3)*exp(-2x-x*x-(1+y)^2)*(2*exp(2y)*y - 18*exp(2x)*(-1+x)^2*(1+y) -6*exp(1+2x+2y)*y*(-2x + 10*x^3 + 5*y^3*(-5+2*y*y)))
        return dzdx, dzdy
    end
    return oldzfun, olddzfun
end

##############################################################################

function main()
    main1()
    main2()
    main3()
    main4()
    main5()
    main6()
    main7()
    main8()
    main9()
    main10()
    main11()
    main12()
    main13()
    main14()
end
    

function main1()
    x = range(-6, 6, length = 101)
    y = range(-6, 6, length = 101)
    f(x,y) = sqrt(x*x + y*y) - 3
    zmin = -3
    zmax = maximum([f(a,b) for a in x for b in y])
    
    function patchcolor(x,y)
        r=(f(x,y)-zmin)/zmax
        return  Color(r,  0.5,  0.9-r/2)
    end
    me = pk3.Mesh(x, y, f; patchcolor,
                  xlineindices=3:4:100,
                  ylineindices=3:4:100, azimuth=-70, elevation=60)
    d = draw(me)
    save(d, plotpath("mesh1.pdf"))
end    


function main2()
    # 3d peaks plot, ala Matlab
    x = range(-3, 3, length = 49)
    y = range(-3, 3, length = 49)
    f(x,y) =  3*(1-x)^2*exp(-(x^2) - (y+1)^2) - 10*(x/5 - x^3 - y^5)*exp(-x^2-y^2) - 1/3*exp(-(x+1)^2 - y^2)
    me = pk3.Mesh(x, y, f)
    d = draw(me)
    save(d, plotpath("mesh2.pdf"))
end



function main3()
    f, df = getsurf1()
    rt = pk3.surf(f, df; xmin = -5, xmax = 5, ymin = -5, ymax = 5,
              axisoptions3_azimuth = -25, axisoptions3_elevation = 40)
    d = draw(rt)
    save(d, plotpath("raytrace3.pdf"))
end


function main4()
    f, df = getsurf2()
    rt = pk3.surf(f, df; xmin = -3, xmax = 3, ymin = -3, ymax = 3,
              axisoptions3_azimuth = -25,
              axisoptions3_elevation = atan(1/sqrt(2))*180/pi)
    d = draw(rt)
    save(d, plotpath("raytrace4.pdf"))
end



function main5()
    ell = Ellipsoid(eye(3),[2,-2,1], pk3.Uniform(:red))
    box = Box3(-3,3,-3,3,-3,3)
    rt = Raytrace([ell], box)
    d = draw(rt)
    save(d, plotpath("raytrace5.pdf"))
end

function main6()
    A = [1.0  0.0  0.0;
        -1.0  0.0  0.0;
         0.0  1.0  0.0;
         0.0 -1.0  0.0;
         0.0  0.0  1.0;
         0.0  0.0 -1.0;]
    b = 0.75* [1.0,1.0,1.0,1.0,1.0,1.0]
    poly1 = Polytope(A, b, fill(Uniform(:blue), 6))
    
    A2 =[ 1  0  0;
         -1  0  0;
          0  1  0;
          0 -1  0;
          0  0  1;
          0  0 -1;]
    b2 = [1,1,1,1,4,4]
    poly2 = Polytope(A2, b2, fill(Uniform(:green), 6))

    A3 =[ 1.0  1.0  0.0;
         -1.0 -1.0  0.0;
          1.0 -1.0  0.0;
         -2.0  1.0  0.2;
          0.0  -0.2  1.0;
          0.0  0.0 -1.0;]
    b3 = [1,1,1,1,3,3]
    poly3 = Polytope(A3, b3, fill(Uniform(:red), 6))

    A4 = [1 1 1;
          1 1 -1;
          1 -1 1 ;
          1 -1 -1;
          -1 1 1;
          -1 1 -1;
          -1 -1 1 ;
          -1 -1 -1]
    b4 = [1,1,1,1,1,1,1,1]
    x4center = [1,0,1]
    b4 = b4 + A4*x4center
    poly4 = Polytope(A4, b4, fill(Uniform(:green), 8))
          
    ell = Ellipsoid(eye(3),[2,-2,1], Uniform(:white))
    ell2 = Ellipsoid(eye(3),[0,-2,-1], Uniform(:white))
    
    shapes = [ poly1, poly4, ell, ell2]
    box = Box3(-1,2,-3,3,-3,2)
    rt = Raytrace(shapes, box;  shadows = true)
    d = draw(rt)
    save(d, plotpath("raytrace6.pdf"))
end





function getsurf7()
    zfun(x,y) =  3*(1-x)^2*exp(-(x^2) - (y+1)^2) - 10*(x/5 - x^3 - y^5)*exp(-x^2-y^2) - 1/3*exp(-(x+1)^2 - y^2)
    function dzfun(x,y)
        dzdx = -(2/3)*exp(-2x-x*x-(1+y)^2)*(-exp(2y)*(1+x) + 9*exp(2x)*(1 - 2*x*x + x*x*x) + exp(1 + 2x + 2y)*(3 - 51x*x + 30*x^4 + 30*x*y^5))
        dzdy =  (1/3)*exp(-2x-x*x-(1+y)^2)*(2*exp(2y)*y - 18*exp(2x)*(-1+x)^2*(1+y) -6*exp(1+2x+2y)*y*(-2x + 10*x^3 + 5*y^3*(-5+2*y*y)))
        return dzdx, dzdy
    end
    return zfun, dzfun
end

# main
function main7()
    # 3d peaks plot, ala Matlab
    f, df = getsurf7()
    rt = surf(f, df; xmin = -3, xmax = 3, ymin = -3, ymax = 3,
              azimuth = -25,elevation =  atan(1/sqrt(2))*180/pi)
    d = draw(rt)
    save(d, plotpath("raytrace7.pdf"))
end


# draw a Polytope
function main8()
    A = [1.0 1 1;
         1 1 -1;
         1 -1 1 ;
         1 -1 -1;
         -1 1 1;
         -1 1 -1;
         -1 -1 1 ;
         -1 -1 -1]
    b = [1,1,1,1,1,1,1,1]
    xcenter = [1,0,1]
    b = b + A*xcenter
    poly = Polytope(A, b, fill(Uniform(:red), 8))
    rt = Raytrace([poly], Box3(-1,2,-3,3,-3,2); shadows = true)
    d = draw(rt)
    save(d, plotpath("raytrace8.pdf"))
end


# draw an ellipsoid
function main9()
    Z = [4.8 1.2 -3.6; 1.2 4.8 -1.4; -3.6 -1.4 4.4]
    ell = Ellipsoid(Z,[0,-2,1], Uniform(:red))
    rt = Raytrace([ell], Box3(-1,2,-3,3,-3,2) ;  shadows = true)
    d = draw(rt)
    save(d, plotpath("raytrace9.pdf"))
end




##############################################################################
# test8

function gettv()
    tvinside(a) = a.x^4 + a.y^4 + a.z^4 <= 1
    function tvnormal(a)
        q =  Vec3(4*a.x^3, 4*a.y^3, 4*a.z^3)
        return q, 1
    end
    return tvinside, tvnormal
end


# draw a tvscreen
function main10()
    f, df = gettv()
    arb = ArbitrarySolid(f, df, 1e-3, Uniform(:red), Uniform(:green))
    rt = Raytrace([arb], Box3(-1,2,-3,3,-3,2); shadows = true)
    d = draw(rt)
    save(d, plotpath("raytrace10.pdf"))
end


# high res
function main11()
    f, df = getsurf1()
    rt = surf(f, df; xmin = -5, xmax = 5, ymin = -5, ymax = 5,
              raytrace_renderwidth = 1600, raytrace_renderheight = 1200,
              axisoptions3_azimuth = -25, axisoptions3_elevation = 40)
    d = draw(rt)
    save(d, plotpath("raytrace11.pdf"))
end


# aspect ratio
function main12()
    f, df = getsurf1()
    rt = surf(f, df; xmin = -5, xmax = 5, ymin = -5, ymax = 5,
           axisoptions3_azimuth = -25, axisoptions3_elevation = 40,
           raytrace_renderwidth = 1600, raytrace_renderheight = 1200,
           axisoptions3_cube = Box3(-2,2, -2,2, -1, 1), axisoptions3_scale3=0.5
              )
    d = draw(rt)
    save(d, plotpath("raytrace12.pdf"))
end

# simple plane
function main13()
    # plot a'x = 0
    a = -2*[-0.3, 0.3, -1]
    c1 = -a[1]/a[3]
    c2 = -a[2]/a[3]
   
    f(x, y) =  c1*x + c2*y
    df(x, y) = c1, c2

    w = 6
    b3 = Box3(-w, w, -w, w, -w, w)
    t1 = Grid(Material(:red), Material(:white), 0.02, 1/3, 1/3)
    t2 = Grid(Material(:red), Material(:white), 0.02, 1/3, 1/3)

    surface = Surface(f, df; texture1=t1, texture2=t2)
    rt = Raytrace([surface], b3;
                  axisoptions3_azimuth = -25,
                  axisoptions3_elevation = 40,
                  renderwidth = 1600,
                  renderheight = 1600,
                  axisoptions3_cube = Box3(-1,1, -1,1, -1,1),
                  axisoptions3_scale3=0.72
                  )
    d = draw(rt)

    am3 = rt.axis3.axismap3
    f(q)  = ctxfromaxis(am3, q)
    v3(x,y,z) = Vec3(x, y, z)
    circle(d.ctx, f(v3(0,0,0)), 3; fillcolor = Color(:black))
    arrow = TriangularArrow(size = 5, fillcolor = Color(:white))
    arrows = ((1, arrow),)
    p = Path(; arrows, points = [f(v3(0,0,0)), f(v3(a...))])
    draw(d, p)

    save(d, plotpath("raytrace13.pdf"))
    return rt
end





# simple plane, non-raytrace
function main14()
    # plot a'x = 0
    a = -4*[-0.3, 0.3, -1]
    c1 = -a[1]/a[3]
    c2 = -a[2]/a[3]
    gs(v) = v - a*(a'*v)/(norm(a)^2)
    v1 = gs(2*[2.5,-2.5,0])
    v2 = gs(2*[1,2,0])
    
    f(x, y) =  c1*x + c2*y
    w = 6
    x = range(-w, w, step = 2)
    y = range(-w, w, step = 2)

    rt = Mesh(x, y, f; patchcolor = Color(0.8,0.8,1),
                  width = 600, height = 600,
                  azimuth = -25, elevation = 40,windowbackgroundcolor = Color(:white),
                  drawbackground = false,
                  axisstyle3_axisbackgroundcolor = Color(0.9176, 0.9176, 0.95),
                  axisstyle3_gridlinestyle = LineStyle(0.8 * Color(:white), 1),
                  axisbox_zmin = -w, axisbox_zmax = w,
                  tickbox_zmin = -w, tickbox_zmax = w)
    d = draw(rt)

    v3(x) = Vec3(x[1], x[2], x[3])

    drawvec(v, c) = draw(d, arrow([Vec3(0,0,0), v3(v)]; center = true,
                                  linestyle = LineStyle(c,1), fillcolor = c))
    drawvec(a, Color(:blue))
    drawvec(v1, Color(:red))
    drawvec(v2, Color(:red))
    circle(d, Vec3(0,0,0), 2; fillcolor = Color(:black))

    save(d, plotpath("raytrace14.pdf"))
    return rt
end


    

end

