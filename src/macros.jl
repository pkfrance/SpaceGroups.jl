# A place for macros used internally in the package

"""
    @delegate TypeName.field (f1, f2, ...)

Generate thin forwarding methods for each function name in `targets`, delegating
calls on `TypeName` instances to the specified field. 
"""
macro delegate(source, targets)
    # source is like SpaceGroupQuotient{N, T}.group
    # or just SpaceGroupQuotient.group
    
    type_expr = source.args[1] 
    field = source.args[2] isa QuoteNode ? source.args[2].value : source.args[2]
    
    # Extract functions to delegate
    funcs = targets isa Expr && targets.head == :tuple ? targets.args : [targets]
    
    defs = Any[]
    for f in funcs
        push!(defs, quote
            @inline function $(esc(f))(obj::$(esc(type_expr)), args...; kwargs...)
                return $(esc(f))(getfield(obj, $(QuoteNode(field))), args...; kwargs...)
            end
        end)
    end
    return Expr(:block, defs...)
end