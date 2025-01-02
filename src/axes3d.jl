
module Axes3D

using PlotKit
using LinearAlgebra

using ..Basic3D: Vec3, Array23, Array32, dot, vec3_hadamard, vec3_hadamarddiv

export Axis3, AxisMap3, Box3, Box3f, drawaxis3, ctxfromcube, axisfromcube
export cubefromctx, ctxfromaxis, AxisOptions3
export get_hexagon

const pk = PlotKit

ifnotnothing(x::Nothing, y) = y
ifnotnothing(x, y) = x

hadamard(p::Vec3, q::Vec3) = vec3_hadamard(p,q)
hadamarddiv(p::Vec3, q::Vec3) = vec3_hadamarddiv(p,q)

mutable struct Box3 
    xmin
    xmax
    ymin
    ymax
    zmin
    zmax
end

Box3(;
     xmin=missing, xmax=missing,
     ymin=missing, ymax=missing,
     zmin=missing, zmax=missing) = Box3(xmin, xmax, ymin, ymax, zmin, zmax)

function PlotKit.ifnotmissing(a::Box3, b::Box3)
    return Box3(ifnotmissing(a.xmin, b.xmin),
               ifnotmissing(a.xmax, b.xmax),
               ifnotmissing(a.ymin, b.ymin),
               ifnotmissing(a.ymax, b.ymax),
               ifnotmissing(a.zmin, b.zmin),
               ifnotmissing(a.zmax, b.zmax)
               )
end

#box3(a...; kw...) = Box3(a...; kw...)

struct Box3f
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    zmin::Float64
    zmax::Float64
end

Box3f(a::Box3) = Box3f(a.xmin, a.xmax, a.ymin, a.ymax, a.zmin, a.zmax)
Base.convert(::Type{Box3f}, x) = Box3f(x)
Base.convert(::Type{Box3f}, x::Box3f) = x




##############################################################################


#
# contains the functional info needed to draw on the axes
#
struct AxisMap3
    depth::Vec3
    up::Vec3
    origin::Vec3              # axis origin in cube coordinates
    scale::Vec3               # ratio of  cube coords to axis coords 
    projector::Array23        # projector matrix P maps cube to unit square 2x3
    projectorinv::Array32
    axis::Axis                # 2d axis from PlotKit
    cube::Box3
end




Base.@kwdef mutable struct Ticks3
    xticks = missing
    xtickstrings = missing
    yticks = missing
    ytickstrings = missing
    zticks = missing
    ztickstrings = missing
end

function PlotKit.ifnotmissing(a::Ticks3, b::Ticks3)
    return Ticks3(ifnotmissing(a.xticks, b.xticks),
                  ifnotmissing(a.xtickstrings, b.xtickstrings),
                  ifnotmissing(a.yticks, b.yticks),
                  ifnotmissing(a.ytickstrings, b.ytickstrings),
                  ifnotmissing(a.zticks, b.zticks),
                  ifnotmissing(a.ztickstrings, b.ztickstrings))                  
end


Base.@kwdef mutable struct AxisStyle3
    drawaxisbackground = true
    axisbackgroundcolor = Color(:white)
    axisedgelinestyle = LineStyle(Color(:black), 1)
    faraxisedgelinestyle = LineStyle(0.87 * Color(:white), 1)
    drawaxisedges = true
    drawfaraxisedges = false
    gridlinestyle = LineStyle(0.87 * Color(:white), 1)
    fontsize = 10
    fontname = "Sans"
    xtickoffset = nothing
    ytickoffset = nothing
    ztickoffset = nothing
    fontcolor = Color(:black)
    labelfontsize = nothing
    labelfontname = nothing
    labelfontcolor = nothing
    xlabel = "x"
    xlabeloffset = nothing
    ylabel = "y"
    ylabeloffset = nothing
    zlabel = "z"
    zlabeloffset = nothing
    
end
#
# we use this to draw the window, in addition to the axis.
#
mutable struct Axis3
    width
    height
    axismap3::AxisMap3
    ax2::Axis
    box3::Box3
    ticks3::Ticks3
    as3::AxisStyle3
    windowbackgroundcolor
    drawwindowbackground
end


Base.@kwdef mutable struct AxisOptions3
    xmin = -Inf
    xmax = Inf
    ymin = -Inf
    ymax = Inf
    zmin = -Inf
    zmax = Inf
    width = 800
    height = 600
    lmargin = 80
    rmargin = 80
    tmargin = 80
    bmargin = 80
    xidealnumlabels = 6
    yidealnumlabels = 6
    zidealnumlabels = 6
    windowbackgroundcolor = Color(0.9176, 0.9176, 0.95)
    drawbackground = true
    ticks = Ticks3()
    axisstyle3 = AxisStyle3()
    tickbox = Box3()
    axisbox = Box3()
    azimuth = -37.5
    elevation = 30
    scale3 = 0.7
    cube = Box3(-1,1,-1,1,-1,1)
end




##############################################################################

Ticks3(b::Box3, xidl, yidl, zidl) = Ticks3(b.xmin, b.xmax,
                                         b.ymin, b.ymax,
                                         b.zmin, b.zmax,
                                         xidl, yidl, zidl)

function Ticks3(xmin, xmax, ymin, ymax, zmin, zmax,
                xidealnumlabels, yidealnumlabels, zidealnumlabels)
    xt = best_ticks(xmin, xmax, xidealnumlabels)
    yt = best_ticks(ymin, ymax, yidealnumlabels)
    zt = best_ticks(zmin, zmax, zidealnumlabels)
    xl = best_labels(xt)
    yl = best_labels(yt)
    zl = best_labels(zt)
    ticks =  Ticks3(xt, xl, yt, yl, zt, zl)
    return  ticks
end


function get_tick_extents3(t::Ticks3)
    xmin, xmax, ymin, ymax, zmin, zmax = minimum(t.xticks), maximum(t.xticks), minimum(t.yticks), maximum(t.yticks), minimum(t.zticks), maximum(t.zticks)
    return Box3(xmin, xmax, ymin, ymax, zmin, zmax)
end

# compute c, t so that
#
#  cubecoords = c + t * axiscoords
#
function cube_axis(axismin, axismax, cubemin, cubemax)
    t = (cubemax - cubemin) / (axismax - axismin)
    c = cubemin - t * axismin
    return c, t
end

function compute_coeffs_for_map_between_cube_and_axis_coords(axisbox, cube::Box3)
    cx, tx = cube_axis(axisbox.xmin, axisbox.xmax, cube.xmin, cube.xmax)
    cy, ty = cube_axis(axisbox.ymin, axisbox.ymax, cube.ymin, cube.ymax)
    cz, tz = cube_axis(axisbox.zmin, axisbox.zmax, cube.zmin, cube.zmax)

    origin = Vec3(cx, cy, cz)
    scale = Vec3(tx, ty, tz)
    return origin, scale
end

function old_compute_coeffs_for_map_between_cube_and_axis_coords(axisbox)
    tx = 2/(axisbox.xmax-axisbox.xmin)
    cx = -tx*axisbox.xmin-1
    ty = 2/(axisbox.ymax-axisbox.ymin)
    cy = -ty*axisbox.ymin-1
    tz = 2/(axisbox.zmax-axisbox.zmin)
    cz = -tz*axisbox.zmin-1

    origin = Vec3(cx, cy, cz)
    scale = Vec3(tx, ty, tz)
    return origin, scale
end

# follow matlab:    y left,  z up, x in
function AxisMap3(axis2, axisbox::Box3, azimuth, elevation, scale3, cube)

    # map 3d unit cube to 2d unit plane
    P, Pinv, depth, up = projectors(azimuth, elevation, scale3)

    # maps between axis coords and unit cube
    origin, scale = compute_coeffs_for_map_between_cube_and_axis_coords(axisbox, cube)
    
    return AxisMap3(depth, up, origin, scale, P, Pinv, axis2, cube)

end

# databox defines the minimum area which the ticks
# are guaranteed to contain.
function Axis3(databox::Box3, ao3::AxisOptions3)
    #println("databox = ", databox)
    #println("ao3.tickbox = ", ao3.tickbox)

    tickbox = ifnotmissing(ao3.tickbox, databox)
    #println("tickbox = ", tickbox)
    #println("ao3.ticks = ", ao3.ticks)
    ticks3 = ifnotmissing(ao3.ticks, Ticks3(tickbox,
                                            ao3.xidealnumlabels,
                                            ao3.yidealnumlabels,
                                            ao3.zidealnumlabels))

    #println((;ticks3))
    axisbox3 = ifnotmissing(ao3.axisbox, get_tick_extents3(ticks3))
    #println("axisbox3 = ", axisbox3)
    
    axis2 = Axis(; axisequal = true,
                 width = ao3.width, height = ao3.height,
                 windowbackgroundcolor = 0.9*Color(:white),
                 lmargin = ao3.lmargin,
                 rmargin = ao3.rmargin,
                 tmargin = ao3.tmargin,
                 bmargin = ao3.bmargin,
                 axisbox_xmin = -1, axisbox_xmax = 1,
                 axisbox_ymin = -1, axisbox_ymax =1)

    axismap3 = AxisMap3(axis2, axisbox3, ao3.azimuth, ao3.elevation, ao3.scale3, ao3.cube)
    
    axis3 = Axis3(ao3.width, ao3.height, axismap3, axis2,
                  axisbox3, ticks3, ao3.axisstyle3,
                  ao3.windowbackgroundcolor,
                  ao3.drawbackground)
                  
    return axis3
end


function parse_axisoptions3(; kw...)
    ao3 = AxisOptions3()
    setoptions!(ao3, "", kw...)
    setoptions!(ao3.tickbox, "tickbox_", kw...)
    setoptions!(ao3.axisbox, "axisbox_", kw...)
    setoptions!(ao3.ticks, "ticks_", kw...)
    setoptions!(ao3.axisstyle3, "axisstyle3_", kw...)
    return ao3
end
Axis3(p; kw...) = Axis3(p, parse_axisoptions3(; kw...))
Axis3(; kw...) = Axis3(missing; kw...)


    


function drawaxis3(ctx, axis3::Axis3)
    # if axis3.drawwindowbackground
    #     rect(ctx, Point(0,0), Point(axis3.width, axis3.height);
    #          fillcolor=axis3.windowbackgroundcolor)
    # end
    drawaxis3(ctx, axis3.axismap3, axis3.ax2, axis3.box3, axis3.ticks3,
              axis3.as3)
end

##############################################################################



function azel(azimuth, elevation)
    deg = pi/180
    # azimuth rotate about z
    Raz = [ cos(azimuth*deg) sin(azimuth*deg) 0  ;
           -sin(azimuth*deg) cos(azimuth*deg) 0  ;
           0 0 1 ]

    # elevation rotate about new horizontal, would be x if az=0
    Rel = [1 0 0 ;
           0 cos(elevation*deg) -sin(elevation*deg);
           0 sin(elevation*deg) cos(elevation*deg)]

    return Raz, Rel
end

# we have three coordinate systems
#
#   1.  axis coordinates (x,y,z as marked on the 3d axis)
#   2.  cube coordinates (between -1 and 1)
#   3.  screenaligned 3d coordinates, with x horizonal, z vertical, and y into the screen
#   4.  square coordinates in the 2d plane of the figure, between -1,1, axis equal
#   5.  ctx coordinates. cairo coordinates
#
# The important thing about cube coordinates is that they are "axis_equal"
# That is, axis3.projector is a constant multiplied by
# a 2 by 3 matrix with orthonormal rows
# 
#
#   cubefromaxis                 maps 1 --> 2    scale then shift (origin, scale)
#   R                            maps 2 --> 3    rotation matrix
#   P2 = scale3*[1 0 0; 0 0 1]   maps 3 --> 4    just keep scaled x, z components
#   axis3.projector              maps 2 --> 4    P = scale3 * P2 * R, matrix multiplication
#   axis3.axis.ax.f              maps 4 --> 5    the usual 2d axis
#


# P maps cube coordinates to the 2d unit square 
#
function projectors(azimuth, elevation, scale3)
    Raz, Rel = azel(azimuth, elevation)
    
    # The rotation matrix R maps from (scaled) cube coordinates to 
    # 3d coordinates aligned with the square coordinates
    # with x horizontal, z vertical, and y into the screen
    R = Rel * Raz
    
    # The square coordinates are just the x-z components of the
    # above aligned system
    P2 = [1 0 0 ;
          0 0 1]

    # The vector [0,1,0] points into screen in screen-aligned coords
    # The vector 'depth' is in cube coordinates
    depth = R'*[0,1,0]


    # The vector [0,0,1] points up in screen-aligned coords
    # The vector 'up' is in cube coordinates
    up = R'*[0,0,1]

    
    P = scale3 * P2 * R
    
    # Pinv is really just P' scaled, since P is a 
    # (scaled) matrix  with orthonormal rows
    Pinv = pinv(P)
   
    return P, Pinv, depth, up

end



cubefromaxis(a::AxisMap3, q::Vec3) = hadamard(a.scale, q) + a.origin
axisfromcube(a::AxisMap3, q::Vec3) = hadamarddiv(q - a.origin, a.scale)

cubefromctx(a::AxisMap3, q::Point) = a.projectorinv * a.axis.ax.finv(q)
ctxfromcube(a::AxisMap3, q::Vec3) = a.axis.ax(a.projector * q)

ctxfromaxis(a::AxisMap3, q::Vec3) = ctxfromcube(a, cubefromaxis(a, q))

(a::AxisMap3)(p::Vec3) = ctxfromaxis(a, p)
(a::AxisMap3)(plist::Array{Vec3}) = a.(plist)


function get_hexagon(axismap3)
    f(xa, ya, za)  = ctxfromcube(axismap3, Vec3(xa, ya, za))
    cube = axismap3.cube
    q1 = f(cube.xmin, cube.ymin, cube.zmin)
    q2 = f(cube.xmax, cube.ymin, cube.zmin)
    q3 = f(cube.xmax, cube.ymin, cube.zmax)
    q4 = f(cube.xmax, cube.ymax, cube.zmax)
    q5 = f(cube.xmin, cube.ymax, cube.zmax)
    q6 = f(cube.xmin, cube.ymax, cube.zmin)
    q7 = f(cube.xmin, cube.ymin, cube.zmax)
    q8 = f(cube.xmax, cube.ymax, cube.zmin)
    if q5.x > q1.x
        return [q1, q2, q8, q4, q5, q7]
    end
    return [q1, q2, q3, q4, q5, q6]
end

function drawaxis3(ctx, axismap3, axis2, box, ticks, as3::AxisStyle3)

    f(a, b, c)  = ctxfromaxis(axismap3, Vec3(a,b,c))
    g(z) = ctxfromaxis(axismap3, z)
    xmin = box.xmin
    xmax = box.xmax
    ymin = box.ymin
    ymax = box.ymax
    zmin = box.zmin
    zmax = box.zmax

    xticks = ticks.xticks
    xtickstrings = ticks.xtickstrings
    yticks = ticks.yticks
    ytickstrings = ticks.ytickstrings
    zticks = ticks.zticks
    ztickstrings = ticks.ztickstrings

    # zfar is the z coordinate of the bottom face of the cube
    # i.e., the one behind the plot

    if axismap3.depth.z > 0
        zfar = zmax
        znear = zmin
    else
        zfar = zmin
        znear = zmax
    end

    if axismap3.depth.y > 0
        yfar = ymax
        ynear = ymin
    else
        yfar = ymin
        ynear = ymax
    end

    if axismap3.depth.x > 0
        xfar = xmax
        xnear = xmin
    else
        xfar = xmin
        xnear = xmax
    end

    # each axis must be at the intersection of a back face and a front face
    # pick leftmost vertical for z axis
    xnearleft = false
    if f(xnear,yfar,zmin).x < f(xfar,ynear,zmin).x
        # xnear is to the left
        zaxis_x = xnear
        zaxis_y = yfar
        xnearleft = true
    else
        zaxis_x = xfar
        zaxis_y = ynear
    end
   
    
    # draw the faces of the cube
    if as3.drawaxisbackground
        zface = Vec3[(xmin,ymin,zfar), (xmax,ymin,zfar), (xmax,ymax,zfar), (xmin,ymax,zfar)]
        line(ctx, g.(zface); closed=true, fillcolor = as3.axisbackgroundcolor)
        
        yface = Vec3[(xmin,yfar,zmin), (xmax,yfar,zmin), (xmax,yfar,zmax), (xmin,yfar,zmax)]
        line(ctx, g.(yface); closed=true, fillcolor = as3.axisbackgroundcolor)

        xface = Vec3[(xfar,ymin,zmin), (xfar,ymax,zmin), (xfar,ymax,zmax), (xfar,ymin,zmax)]
        line(ctx, g.(xface); closed=true, fillcolor = as3.axisbackgroundcolor)
    end


    function tickline(x1, y1, z1, x2, y2, z2)
        line(ctx, g.(Vec3[(x1,y1,z1), (x2,y2,z2)]), linestyle = as3.gridlinestyle)
    end

    
    # draw the lines on the faces of the cube that correspond to ticklines
    for xt in xticks
        if xmin < xt < xmax
            tickline(xt, ymax, zfar, xt, ymin, zfar)
            tickline(xt, yfar, zmin, xt, yfar, zmax)
        end
    end
    for yt in yticks
        if ymin < yt < ymax
            tickline(xmin, yt, zfar, xmax, yt, zfar)
            tickline(xfar, yt, zmin, xfar, yt, zmax)
        end
    end
    
    for zt in zticks
        if zmin < zt < zmax
            tickline(xmin, yfar, zt, xmax, yfar, zt)
            tickline(xfar, ymin, zt, xfar, ymax, zt)
        end
    end
  
            
    if as3.drawaxisedges
        linestyle = as3.axisedgelinestyle
        # x axis is drawn at nearest y
        line(ctx, g.(Vec3[(xmin, ynear, zfar), (xmax, ynear, zfar)]); linestyle)
        # y
        line(ctx, g.(Vec3[(xnear, ymin, zfar), (xnear, ymax, zfar)]); linestyle)
        # z
        line(ctx, g.(Vec3[(zaxis_x, zaxis_y, zmin), (zaxis_x, zaxis_y, zmax)]); linestyle)
    end

    if as3.drawfaraxisedges
        linestyle = as3.faraxisedgelinestyle
        # x axis is drawn at far y
        line(ctx, g.(Vec3[(xmin, yfar, zfar), (xmax, yfar, zfar)]); linestyle)
        # y
        line(ctx, g.(Vec3[(xfar, ymin, zfar), (xfar, ymax, zfar)]); linestyle)
        # z
        line(ctx, g.(Vec3[(xfar, yfar, zmin), (xfar, yfar, zmax)]); linestyle)
    end


    

    function exttickline(q::Vec3, offset)
        a = ctxfromaxis(axismap3, q)
        b = ctxfromcube(axismap3, cubefromaxis(axismap3, q) + offset)
        line(ctx, [a, b]; linestyle = LineStyle(Color(:black), 0.5))
    end


    if isnothing(as3.xtickoffset) 
        xtickoffset = Vec3(0, sign(ynear-yfar)*0.25, 0)
    else
        xtickoffset = as3.xtickoffset 
    end
    if isnothing(as3.ytickoffset)
        ytickoffset = Vec3(sign(xnear-xfar)*0.15, 0, 0)
    else
        ytickoffset = as3.ytickoffset
    end
    if isnothing(as3.ztickoffset)
        # now we choose ticks on z-axis to be in horizontal plane
        if xnearleft
            # offset in y direction
            ztickoffset = Vec3(0, -0.12*sign(ynear-yfar), 0)
        else
            # offset in x direction
            ztickoffset = Vec3(-0.12*sign(xnear-xfar), 0, 0)
        end
    else
        ztickoffset = as3.ztickoffset
    end

    

    # draw the font
    for i = 1:length(xticks)
        xt = xticks[i]
        if xmin < xt < xmax
            pos = cubefromaxis(axismap3, Vec3(xt, ynear, zfar)) + xtickoffset
            text(ctx, ctxfromcube(axismap3, pos),
                 as3.fontsize, as3.fontcolor, xtickstrings[i];
                 fname = as3.fontname,
                 horizontal = "center")
            exttickline(Vec3(xt, ynear, zfar), 0.4*xtickoffset)
        end
    end
    
    for i = 1:length(yticks)
        yt = yticks[i]
        if ymin < yt < ymax
            pos = cubefromaxis(axismap3, Vec3(xnear, yt, zfar)) + ytickoffset
            text(ctx, ctxfromcube(axismap3, pos),
                 as3.fontsize, as3.fontcolor, ytickstrings[i];
                 fname = as3.fontname,
                 horizontal = "center", vertical = "center")
            exttickline(Vec3(xnear, yt, zfar), 0.4*ytickoffset)
        end
    end
    
    for i = 1:length(zticks)
        zt = zticks[i]
        if zmin < zt < zmax
            pos = cubefromaxis(axismap3, Vec3(zaxis_x, zaxis_y, zt)) + ztickoffset
            text(ctx, ctxfromcube(axismap3, pos),
                 as3.fontsize, as3.fontcolor, ztickstrings[i];
                 fname = as3.fontname,
                 vertical = "center", horizontal = "right")
            exttickline(Vec3(zaxis_x, zaxis_y, zt), 0.4*ztickoffset)
        end
    end

    # draw the labels
    labelfontcolor = ifnotnothing(as3.labelfontcolor, as3.fontcolor)
    labelfontsize  = ifnotnothing(as3.labelfontsize, as3.fontsize)
    labelfontname  = ifnotnothing(as3.labelfontname, as3.fontname)

    xlabeloffset   = ifnotnothing(as3.xlabeloffset,  xtickoffset * 2)
    xt = (box.xmax + box.xmin)/2
    pos = cubefromaxis(axismap3, Vec3(xt, ynear, zfar)) + xlabeloffset
    text(ctx, ctxfromcube(axismap3, pos),
         labelfontsize, labelfontcolor, as3.xlabel;
         fname = labelfontname,
         horizontal = "center")
    
    ylabeloffset   = ifnotnothing(as3.ylabeloffset,  ytickoffset * 2)
    yt = (box.ymax + box.ymin)/2
    pos = cubefromaxis(axismap3, Vec3(xnear, yt, zfar)) + ylabeloffset
    text(ctx, ctxfromcube(axismap3, pos),
         labelfontsize, labelfontcolor, as3.ylabel;
         fname = labelfontname,
         horizontal = "center", vertical = "center")

    zlabeloffset   = ifnotnothing(as3.zlabeloffset,  ztickoffset * 3)
    zt = (box.zmax + box.zmin)/2
    pos = cubefromaxis(axismap3, Vec3(zaxis_x, zaxis_y, zt)) + zlabeloffset
    text(ctx, ctxfromcube(axismap3, pos),
         labelfontsize, labelfontcolor, as3.zlabel;
         fname = labelfontname,
         vertical = "center", horizontal = "right")

end





end

