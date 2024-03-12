
##############################################################################
# the abstract type which we can see

module RayTracer

# from here
export Camera, raytrace_main


#export Shape, Texture, Grid, transform

using PlotKit

using ..Axes3D: Box3, Box3f, AxisMap3
using ..CoreRT
using ..RTloop

##############################################################################

##############################################################################

##############################################################################
# Camera

struct Camera
    am3::AxisMap3
    xmin::Int64
    xmax::Int64
    ymin::Int64
    ymax::Int64
    width::Int64    # of the context ctx
    height::Int64
    renderwidth::Int64   # of the pixel array used for raytracing
    renderheight::Int64
end

function cubefromarray(camera::Camera, z::Point)
    w = Point(z.x * camera.width / camera.renderwidth,
              z.y * camera.height / camera.renderheight)
    return cubefromctx(camera.am3, w)
end

# this returns a box containing the 2d region
# we want to render, in array coords
function get_bounding_box(am3::AxisMap3, width, height, renderwidth, renderheight)
    hexagon = get_hexagon(am3)
    xvals = [a.x * renderwidth / width for a in hexagon]
    xmin = Int64(floor(minimum(xvals)))
    xmax = Int64(ceil(maximum(xvals)))
    ymin = 2
    ymax = renderheight-1
    return xmin, xmax, ymin, ymax
end

function Camera(am3::AxisMap3, width, height, renderwidth, renderheight)
    xmin, xmax, ymin, ymax = get_bounding_box(am3::AxisMap3, width, height,
                                              renderwidth, renderheight)
    return Camera(am3, xmin, xmax, ymin, ymax, width, height,
                  renderwidth, renderheight)
end


# create ray originating at x,y
function makeray(camera::Camera, x, y)
    # 1.74 is the distance from the cube corner to the center
    raydirection = camera.am3.depth
    rayorigin = cubefromarray(camera, Point(x,y)) - 1.74*raydirection
    return rayorigin, raydirection
end

        
##############################################################################



                            
# compute the color seen by a ray
# also returns needsrefinement
function tracefromxy(x, y, camera::Camera,
                     surfacelist, lighting, shadows)
    # in cube coordinates
    rayorigin, raydirection = makeray(camera, x, y)
    cube = camera.am3.cube
    return traceray(rayorigin, raydirection, cube, surfacelist, lighting, 1, false, shadows)
end

function subcheck(X,x,y,n)
    d = 0.1
    if abs(X[x,y,n] - X[x-1,y,n]) > d || abs(X[x,y,n] - X[x+1,y,n]) > d ||
        abs(X[x,y,n] - X[x,y-1,n]) > d || abs(X[x,y,n] - X[x,y+1,n]) > d
        return true
    end
    return false
end

checkvariation(X,x,y) = subcheck(X,x,y,1) || subcheck(X,x,y,2) || subcheck(X,x,y,3)


function raytrace_main(camera::Camera,
                       lighting,
                       surfacelist,
                       refine, shadows) 


    X = ones(camera.renderwidth, camera.renderheight, 3)
    Y = falses(camera.renderwidth, camera.renderheight)


    for x = camera.xmin:camera.xmax
        for y = camera.ymin:camera.ymax
            needsrefinement, c = tracefromxy(x-1, y-1, camera, surfacelist, lighting, shadows)
            Y[x, y] = needsrefinement
            X[x, y, 1] = c.r
            X[x, y, 2] = c.g
            X[x, y, 3] = c.b
        end
    end

    if refine
        for x = camera.xmin:camera.xmax
            for y = camera.ymin:camera.ymax
                if checkvariation(X,x,y)
                    # refine this one
                    Y[x,y] = true
                end
            end
        end

        m = 1
        m1 = 1/(2*m+1)
        m2 = m1*m1

        for x = camera.xmin:camera.xmax
            for y = camera.ymin:camera.ymax
                if Y[x,y] 
                    # center point known already
                    csum = Color(X[x,y,1], X[x,y,2], X[x,y,3])
                    for i = -m:m
                        for j = -m:m
                        # already done center point
                            if i != 0 || j != 0 
                                xx = x-1 + m1*i
                                yy = y-1 + m1*j

                                needsrefinement, c = tracefromxy(xx, yy, camera, surfacelist,
                                                                 lighting, shadows)
                                csum = csum + c
                            end
                        end
                    end
                    X[x, y, 1] = m2*csum.r
                    X[x, y, 2] = m2*csum.g
                    X[x, y, 3] = m2*csum.b
                end
            end
        end
    end
    return X
end

end
