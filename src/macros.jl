# A place for macros used internally in the package

"""
    @delegate TypeName.field (f1, f2, ...)

Generate thin forwarding methods for each function name in `targets`, delegating
calls on `TypeName` instances to the specified field. 
"""
macro delegate(source, targets)
    # source is an expression like SpaceGroupQuotient.group
    typename = source.args[1]
    # Extract the field name from the QuoteNode or Symbol
    fieldname = source.args[2] isa QuoteNode ? source.args[2].value : source.args[2]
    
    # Handle both single functions and tuples of functions
    funcs = targets isa Expr && targets.head == :tuple ? targets.args : [targets]
    
    defs = Any[]
    for f in funcs
        push!(defs, quote
            @inline $(esc(f))(obj::$(esc(typename)), args...; kwargs...) = 
                $(esc(f))(getfield(obj, $(QuoteNode(fieldname))), args...; kwargs...)
        end)
    end
    return Expr(:block, defs...)
end