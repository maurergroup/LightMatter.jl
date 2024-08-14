#My edits to Dierckx
if isdefined(Base, :get_extension)
    using Dierckx
    using Dierckx: Spline1D
    using Symbolics
    using Symbolics: Num, unwrap, SymbolicUtils
else
    using ..Dierckx
    using ..Dierckx: Spline1D
    using ..Symbolics
    using ..Symbolics: Num, unwrap, SymbolicUtils
end
#import Base: length
(spl::Spline1D)(x::Num) = SymbolicUtils.term(spl, unwrap(x))
(spl::Spline1D)(x::SymbolicUtils.BasicSymbolic{Real}) = SymbolicUtils.term(spl, unwrap(x))
SymbolicUtils.promote_symtype(t::Spline1D, _...) = Real
Base.nameof(spl::Spline1D) = :Interpolation

Spline1D(x::Symbolics.Arr,y::Symbolics.Arr,bc="nearest") = SymbolicUtils.term(Spline1D,x,y,bc)
