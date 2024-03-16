
module Basic3D

using PlotKit
using LinearAlgebra

const pk = PlotKit
export Vec3, Array33, Array32, Array23, Arrayn3, dot, normalize, vec3_hadamard, vec3_hadamarddiv

struct Vec3
    x::Float64
    y::Float64
    z::Float64
end

eval(makevector(Vec3))

pk.vec3(a...; kw...) = Vec3(a...; kw...)
# 3 by 3 matrix
struct Array33
    a::Vec3
    b::Vec3
    c::Vec3
end
Array33(A::Array{Float64,2}) = Array33([z[:] for z in eachcol(A)]...)
Base.:*(A::Array33, p::Vec3) = p.x*A.a + p.y*A.b + p.z*A.c
Base.convert(::Type{Array33}, x) = Array33(x)
Base.convert(::Type{Array33}, x::Array33) = x

# 2 by 3 matrix
struct Array23
    x::Point
    y::Point
    z::Point
end
Base.:*(P::Array23, q::Vec3) = P.x*q.x + P.y*q.y + P.z*q.z
Array23(A::Array{Float64,2}) = Array23([z[:] for z in eachcol(A)]...)
Base.convert(::Type{Array23}, x) = Array23(x)
Base.convert(::Type{Array23}, x::Array23) = x



# 3 by 2 matrix
struct Array32
    x::Vec3
    y::Vec3
end
Base.:*(P::Array32, q::Point) = P.x*q.x + P.y*q.y
Array32(A::Array{Float64,2}) = Array32([z[:] for z in eachcol(A)]...)
Base.convert(::Type{Array32}, x) = Array32(x)
Base.convert(::Type{Array32}, x::Array32) = x

end

