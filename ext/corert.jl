
module CoreRT

using LinearAlgebra
using PlotKit
using ..Basic3D: Vec3
using ..Axes3D: Box3


abstract type Shape end
abstract type Texture end

export Shape,  Texture, Material, Grid, Uniform, Testregion
export nohit, Hitdata, gettexture


# functions which which are public methods of Shape objects
export hitnormal,  hittexture,  hitshape,  init!,   getdistancefromstripe, transform
export addlimits




init!() = return
hitnormal() = return
hittexture() = return
hitshape() = return

getdistancefromstripe() = return
transform() = return


nohit(s) = (-1.0,-1.0,Hitdata(s))

transform(p::S, scale::Vec3, origin::Vec3) where S <: Shape = transform(p,
                                                                     diagm(vec(scale)),
                                                                     vec(origin))

transform(s::S, axis3) where S <: Shape = transform(s, axis3.axismap3.scale,
                                                    axis3.axismap3.origin)


addlimits(s::S, cube::Box3) where S <: Shape = return




mutable struct Hitdata{S<:Shape}
    shape::S
    point::Vec3
    face::Int64
    child::Hitdata   # this is an abstract type
    function Hitdata{S}(shape::S, point::Vec3, face::Int64, child::Hitdata) where {S<:Shape}
        new(shape, point, face, child)
    end
    function Hitdata{S}(shape::S, point::Vec3, face::Int64) where {S<:Shape}
        new(shape, point, face)
    end
    function Hitdata{S}(shape::S)  where {S<:Shape}
        new(shape)
    end
end

Hitdata(shape::S, point, face, hitdata) where {S<:Shape} = Hitdata{S}(shape, point, face, hitdata) 
Hitdata(shape::S, point, face) where {S<:Shape} = Hitdata{S}(shape, point, face) 
Hitdata(shape::S) where {S<:Shape} = Hitdata{S}(shape)

##############################################################################
# textures


struct Material
    color::Color
    kambient::Float64
    kdiffuse::Float64
    kspecular::Float64
    pspecular::Float64
    metalness::Float64
    kreflection::Float64
end

# defaults
Material(c::Color) = Material(c, 0.4, 0.45, 0.35, 8, 0.0, 0.0)
function Material(s::Symbol)
    if s == :black
        return Material(Color(:white), 0.0, 0.0, 0.0, 0, 0.0, 0.0)
    end
    if s == :mirror
        return  Material(0.9 * Color(:white), 0.2, 0.2, 0.0, 8, 0.0, 1.0)
    end
    return Material(Color(s))
end


struct Uniform <: Texture
    material::Material
end

# defaults
Uniform(c::Color) = Uniform(Material(c))
Uniform(s::Symbol) = Uniform(Material(s))

struct Grid <: Texture
    material::Material
    gridlines::Material
    width::Float64  # how thick the lines are
    spacingx::Float64  # how far apart the lines are in cube coords
    spacingy::Float64  # how far apart the lines are in cube coords
end

# defaults
Grid(smat::S, gmat::T) where {S <:Material, T<:Material} = Grid(smat, gmat, 0.01, 0.1, 0.1)
Grid(s, g) = Grid(Material(s), Material(g))
Grid(s) = Grid(Material(s), Material(:black))

struct Testregion <: Texture
    material1::Material
    material2::Material
end

function gettexture(s::Shape, t::Uniform, q::Vec3)
    return t.material, false
end




#
# returns
#    material
#    needsrefinement = true if resampling should occur here, (to avoid stripe jaggies)
#
# public
function gettexture(s::Shape, t::Grid, q::Vec3)
    dfs = getdistancefromstripe(s, q, t.spacingx, t.spacingy)
    needsrefinement = false
    if dfs < t.width*4
        needsrefinement = true
    end
            
    if dfs < t.width
        # stripe color
        return t.gridlines, needsrefinement
    end
    
    return t.material, needsrefinement
end



function gettexture(s::Shape, t::Testregion, q::Vec3)
    a = 0.5
    if 0.35<q[1]<0.4 && -0.6<q[2]<-0.4
        return t.material1, false
    end
    return t.material2, false
end
                  

end
