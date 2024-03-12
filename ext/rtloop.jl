
module RTloop

# actual raytracer code lives here

using PlotKit

using ..Basic3D: Vec3, Array23, Array32, Array33
using ..Axes3D: Box3, Box3f
using ..CoreRT: Texture, Shape, CoreRT, Material


# this code does the actual raytracing
export traceray
export Light, Lighting



##############################################################################
struct Light
    lightpos::Vec3
    lightintensity::Color  # color
end

struct Lighting
    lights::Array{Light,1}
    ambientlight::Color    # color
end

# default lighting
function Lighting()
    light1 = Light(Vec3(-30,0,20), 0.8*Color(:white))
    light2 = Light(Vec3(-30,-20,-30), 0.4 * Color(:white))
    lighting = Lighting([light1,light2], 0.6 * Color(:white))
    return lighting
end


function inshadow(cube, light, point, surfacelist::Array{T,1}) where {T<:Shape}

    tmin = 0.01
    rayorigin = point
    raydirection = normalize(light.lightpos - point)
    nsurfaces = length(surfacelist)

    tenter, tleave, face = intersectcube(rayorigin, raydirection, cube)

    for i = 1:nsurfaces
        t, tl, hitdata = hitshape(surfacelist[i], rayorigin, raydirection, tmin, tleave)
        if t >= 0
            return true
        end
    end
    return false
end

function pixelcolor(cube, m::Material, lighting, viewer, normal, point, surfacelist::Array{T,1}, shadows) where {T<:Shape}

    # material color
    pointcolor = m.color

    
    # ambient

    rho = m.kambient * hadamard(pointcolor, lighting.ambientlight)



    for l in lighting.lights
        # check if in shadow
        if !shadows || !inshadow(cube, l, point, surfacelist)
            
            lightdir = normalize(l.lightpos - point)
            # diffuse reflection
            rho += m.kdiffuse * max(0, dot(normal, lightdir)) * hadamard(l.lightintensity, pointcolor)
            
            # specular reflection
            #S = m.metalness*pointcolor + (1-m.metalness) * Color(:white)
            S = interp(Color(:white),  pointcolor,  m.metalness)
            ss = viewer + lightdir
            h = normalize(ss)
            rho +=  m.kspecular*(max(0,dot(normal,h))^m.pspecular) * hadamard(l.lightintensity,  S)
        end
    end

    return rho
end


function getcolor(s::Shape, hitdata,
                            raydirection, cube, lighting, recdepth, 
                            surfacelist::Array{T,1}, shadows) where {T<:Shape}

    hitpoint = hitdata.point
    normal = hitnormal(s, hitdata)

    #@code_warntype hittexture(s, hitdata)
    #return
    
    material, needsrefinement = hittexture(s, hitdata)
    c = pixelcolor(cube, material, lighting, -raydirection, normal, hitpoint, surfacelist, shadows)
    kreflection = material.kreflection
    if kreflection > 0 && recdepth < 4
        reflectiondir = raydirection - 2*dot(normal, raydirection)*normal
        p = hitpoint + 0.01 * normal


        nr::Bool, c2::Color = traceray(p, reflectiondir, cube, surfacelist,  lighting, recdepth+1, true, shadows)
        if c2.r + c2.g + c2.b < 3
            c = c + kreflection * c2
            needsrefinement = needsrefinement | nr
        end
    end

    return needsrefinement, c
end

function intersectcube(rayorigin, raydirection, cube)
    sx = 2/(cube.xmax - cube.xmin)
    sy = 2/(cube.ymax - cube.ymin)
    sz = 2/(cube.zmax - cube.zmin)
    s = Vec3(sx, sy, sz)
    
    qx = (cube.xmax + cube.xmin)/2
    qy = (cube.ymax + cube.ymin)/2
    qz = (cube.zmax + cube.zmin)/2
    q = Vec3(qx, qy, qz)

    
    origin = vhadamard(rayorigin - q, s)
    direction = vhadamard(raydirection, s)

    te, tl, face = intersectunitcube(origin, direction)
    
    return te, tl, face
end


function traceray(rayorigin, raydirection, cube, surfacelist::Array{T,1},  lighting,
                            recdepth, assumeinside, shadows) where {T<:Shape}



    tenter, tleave, face = intersectcube(rayorigin, raydirection, cube)

    nsurfaces = length(surfacelist)
    tmin::Float64 = Inf

    local imin::Int64
    local hitdata::Hitdata

    if tenter >= 0 #|| assumeinside
        #if assumeinside
        #    tenter = 0.000
        #end
        for i = 1:nsurfaces
            # this ray enters the cube

            t, tl, hd = hitshape(surfacelist[i], rayorigin, raydirection, tenter, tleave)
            if t >= 0 &&  t < tmin
                tmin = t
                imin = i
                hitdata = hd
            end
        end

        if tmin < Inf
            needsrefinement, c = getcolor(surfacelist[imin], hitdata,
                                          raydirection, cube, lighting, 1, 
                                          surfacelist, shadows)
            return needsrefinement, c
        end


    end
    return  false, Color(:white)
end




end
