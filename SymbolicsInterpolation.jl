#My edits to Dierckx
if isdefined(Base, :get_extension)
    using Dierckx
    using Dierckx: Spline1D
    using Symbolics
    using Symbolics: Num, unwrap, SymbolicUtils, has_symwrapper, wrapper_type, wrap, Symbolic
else
    using ..Dierckx
    using ..Dierckx: Spline1D
    using ..Symbolics
    using ..Symbolics: Num, unwrap, SymbolicUtils, has_symwrapper, wrapper_type, wrap, Symbolic
end

(spl::Spline1D)(x::Num) = SymbolicUtils.term(spl, unwrap(x))
(spl::Spline1D)(x::SymbolicUtils.BasicSymbolic{Real}) = SymbolicUtils.term(spl, unwrap(x))
SymbolicUtils.promote_symtype(t::Spline1D, _...) = Real
Base.nameof(spl::Spline1D) = :Interpolation

Spline1D(x::Symbolics.Arr,y::Symbolics.Arr,bc="nearest") = SymbolicUtils.term(Spline1D,x,y,bc)
#= function symbolic_t(T)
    return has_symwrapper(T) ? Union{wrapper_type(T), Symbolic{<:T}} : Symbolic{<:T}
end
R = symbolic_t(Real)
A = symbolic_t(Vector{Real}) =#
