
module Demo3D

using PlotKit
using LinearAlgebra

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


function main4()
    f, df = getsurf2()
    d = surf(f, df; xmin = -3, xmax = 3, ymin = -3, ymax = 3,
             axisoptions3_azimuth = -25, axisoptions3_elevation = atan(1/sqrt(2))*180/pi)

    save(d, plotpath("raytrace4.pdf"))
   
end



function main5()
    ell = ellipsoid(eye(3),[2,-2,1], uniform(:red))
    box = box3(-3,3,-3,3,-3,3)
    d = raytrace([ell], box)
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
    poly1 = polytope(A, b, fill(uniform(:blue), 6))
    
    A2 =[ 1  0  0;
         -1  0  0;
          0  1  0;
          0 -1  0;
          0  0  1;
          0  0 -1;]
    b2 = [1,1,1,1,4,4]
    poly2 = polytope(A2, b2, fill(uniform(:green), 6))

    A3 =[ 1.0  1.0  0.0;
         -1.0 -1.0  0.0;
          1.0 -1.0  0.0;
         -2.0  1.0  0.2;
          0.0  -0.2  1.0;
          0.0  0.0 -1.0;]
    b3 = [1,1,1,1,3,3]
    poly3 = polytope(A3, b3, fill(uniform(:red), 6))

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
    poly4 = polytope(A4, b4, fill(uniform(:green), 8))
          
    ell = ellipsoid(eye(3),[2,-2,1], uniform(:white))
    ell2 = ellipsoid(eye(3),[0,-2,-1], uniform(:white))
    
    shapes = [ poly1, poly4, ell, ell2]
    box = box3(-1,2,-3,3,-3,2)
    d = raytrace(shapes, box;  shadows = true)
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
    d=surf(f, df; xmin = -3, xmax = 3, ymin = -3, ymax = 3,
           azimuth = -25,elevation =  atan(1/sqrt(2))*180/pi)
    save(d, plotpath("raytrace7.pdf"))
end


# draw a polytope
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
    poly = polytope(A, b, fill(uniform(:red), 8))
    d = raytrace([poly], box3(-1,2,-3,3,-3,2); shadows = true)
    save(d, plotpath("raytrace8.pdf"))
end


# draw an ellipsoid
function main9()
    Z = [4.8 1.2 -3.6; 1.2 4.8 -1.4; -3.6 -1.4 4.4]
    ell = ellipsoid(Z,[0,-2,1], uniform(:red))
    d = raytrace([ell], box3(-1,2,-3,3,-3,2) ;  shadows = true)
    save(d, plotpath("raytrace9.pdf"))
end




##############################################################################
# test8

function gettv()
    tvinside(a) = a.x^4 + a.y^4 + a.z^4 <= 1
    function tvnormal(a)
        q =  vec3(4*a.x^3, 4*a.y^3, 4*a.z^3)
        return q, 1
    end
    return tvinside, tvnormal
end


# draw a tvscreen
function main10()
    f, df = gettv()
    arb = arbitrarysolid(f, df, 1e-3, uniform(:red), uniform(:green))
    d = raytrace([arb], box3(-1,2,-3,3,-3,2); shadows = true)
    save(d, plotpath("raytrace10.pdf"))
end


# high res
function main11()
    f, df = getsurf1()
    d=surf(f, df; xmin = -5, xmax = 5, ymin = -5, ymax = 5,
           raytrace_renderwidth = 1600, raytrace_renderheight = 1200,
           axisoptions3_azimuth = -25, axisoptions3_elevation = 40)
    save(d, plotpath("raytrace11.pdf"))
end


# aspect ratio
function main12()
    f, df = getsurf1()
    d=surf(f, df; xmin = -5, xmax = 5, ymin = -5, ymax = 5,
           axisoptions3_azimuth = -25, axisoptions3_elevation = 40,
           raytrace_renderwidth = 1600, raytrace_renderheight = 1200,
           axisoptions3_cube = box3(-2,2, -2,2, -1, 1), axisoptions3_scale3=0.5
           )
    save(d, plotpath("raytrace12.pdf"))
end


end

