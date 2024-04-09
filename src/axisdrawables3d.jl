
module AxisDrawables3D

using PlotKit

using Cairo
using ..Basic3D: Vec3
using ..Axes3D: Axis3, drawaxis3, Axes3D, AxisOptions3, parse_axisoptions3


export AxisDrawable3


abstract type AxisDrawable3 <: Drawable end

Base.@kwdef mutable struct AxisImageDrawable3 <: AxisDrawable3
    surface
    ctx
    width
    height
    fname
    axis3::Axis3
end

Base.@kwdef mutable struct AxisPDFDrawable3 <: AxisDrawable3
    surface
    ctx
    width
    height
    fname
    axis3::Axis3
end

Base.@kwdef mutable struct AxisSVGDrawable3 <: AxisDrawable3
    surface
    ctx
    width
    height
    fname
    axis3::Axis3
end

Base.@kwdef mutable struct AxisRecorderDrawable3 <: AxisDrawable3
    surface
    ctx
    width
    height
    axis3::Axis3
end


PlotKitCairo.Drawable(ad::AxisImageDrawable3) = ImageDrawable(ad.surface, ad.ctx, ad.width, ad.height, ad.fname)
PlotKitCairo.Drawable(ad::AxisPDFDrawable3) = PDFDrawable(ad.surface, ad.ctx, ad.width, ad.height, ad.fname)
PlotKitCairo.Drawable(ad::AxisSVGDrawable3) = SVGDrawable(ad.surface, ad.ctx, ad.width, ad.height, ad.fname)
PlotKitCairo.Drawable(ad::AxisRecorderDrawable3) = RecorderDrawable(ad.surface, ad.ctx, ad.width, ad.height)



AxisDrawable3(axis3::Axis3, dw::ImageDrawable; kw...) = AxisImageDrawable3(; surface = dw.surface, ctx = dw.ctx, width = dw.width, height = dw.height, fname = dw.fname, axis3 = axis3, kw...)
AxisDrawable3(axis3::Axis3, dw::PDFDrawable; kw...) = AxisPDFDrawable3(; surface = dw.surface, ctx = dw.ctx, width = dw.width, height = dw.height, fname = dw.fname, axis3 = axis3, kw...)
AxisDrawable3(axis3::Axis3, dw::SVGDrawable; kw...) = AxisSVGDrawable3(; surface = dw.surface, ctx = dw.ctx, width = dw.width, height = dw.height, fname = dw.fname, axis3 = axis3, kw...)
AxisDrawable3(axis3::Axis3, dw::RecorderDrawable; kw...) = AxisRecorderDrawable3(; surface = dw.surface, ctx = dw.ctx, width = dw.width, height = dw.height, axis3 = axis3, kw...)


# close has different behavior for ImageDrawable vs other Drawables
ImageDrawable(dw::AxisImageDrawable3) = ImageDrawable(dw.surface, dw.ctx, dw.width, dw.height, dw.fname)
PlotKitCairo.close(dw::AxisImageDrawable3) = PlotKitCairo.close(ImageDrawable(dw))

# paint and save only apply to RecorderDrawables
RecorderDrawable(dw::AxisRecorderDrawable3) = RecorderDrawable(dw.surface, dw.ctx, dw.width, dw.height)
PlotKitCairo.paint(ctx::CairoContext, r::AxisRecorderDrawable3, args...) = PlotKitCairo.paint(ctx, RecorderDrawable(r), args...)
PlotKitCairo.save(r::AxisRecorderDrawable3, args...) = PlotKitCairo.save(RecorderDrawable(r), args...)


##############################################################################
# builder functions



function AxisDrawable3(axis3::Axis3; fname = nothing)
    # The call to Drawable starts the interaction with Cairo
    dw = Drawable(axis3.width, axis3.height; fname)
    return AxisDrawable3(axis3, dw)
end

                                    
function AxisDrawable3(p, ao3::AxisOptions3; fname = nothing)
    axis3 = Axis3(p, ao)
    return AxisDrawable3(axis3; fname)
end
                                    
function AxisDrawable3(ao::AxisOptions3; fname = nothing)
    axis3 = Axis3(ao)
    return AxisDrawable3(axis3; fname)
end
                
             



AxisDrawable3(p; fname = nothing, kw...) = AxisDrawable3(p, parse_axisoptions3(; kw...); fname)
AxisDrawable3( ; fname = nothing, kw...) = AxisDrawable3(   parse_axisoptions3(; kw...); fname)



##############################################################################
# drawing axes

function PlotKit.drawbackground(ad::AxisDrawable3)
    if ad.axis3.drawwindowbackground
        rect(ad.ctx, Point(0,0), Point(ad.width, ad.height); fillcolor = ad.axis3.windowbackgroundcolor)
    end
end


# also draw background
function Axes3D.drawaxis3(ad::AxisDrawable3)
    drawbackground(ad)
    drawaxis3(ad.ctx, ad.axis3)
end

##############################################################################
# drawing functions


PlotKitCairo.line(ad::AxisDrawable3, p::Array{Vec3}; kwargs...) =  line(ad.ctx, ad.axis3.axismap3(p); kwargs...)
PlotKitCairo.circle(ad::AxisDrawable3, p, r; kw...) =  circle(ad.ctx, ad.axis3.axismap3(p), r; kw...)
PlotKitCairo.text(ad::AxisDrawable3, p, fsize, color, txt; kw...) =  text(ad.ctx, ad.axis3.axismap3(p), fsize, color, txt; kw...)
PlotKitCairo.linear_pattern(ad::AxisDrawable3, p1::Vec3, p2::Vec3) = PlotKitCairo.linear_pattern(ad.axis3.axismap3(p1), ad.axis3.axismap3(p2))
PlotKitCairo.line(ad::AxisDrawable3, p::Vec3, q::Vec3, args...; kwargs...) = line(ad.ctx, ad.axis3.axismap3(p), ad.axis3.axismap3(q), args...; kwargs...)


##############################################################################
# this depends on PlotKitDiagrams


function PlotKit.draw(ad::AxisDrawable3, x::Vec3, dir::Vec3, arrow::TriangularArrow)
    f(q)  = ad.axis3.axismap3(q)
    dir2 = f(dir) - f(Vec3(0,0,0)) # direction in ctx coords
    draw(Drawable(ad), f(x), dir2, arrow)
end


end
