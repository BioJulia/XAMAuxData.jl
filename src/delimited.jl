"""
    DelimitedIterator(data, x::T)

Return an iterator over memory-backed data `data` of eltype `T`.
Returns `MemoryView`s of the same elements as `data`, split by `x`.
"""
struct DelimitedIterator{T}
    v::ImmutableMemoryView{T}
    d::T
end

function DelimitedIterator(x, d::T) where {T}
    v = ImmutableMemoryView(x)
    eltype(v) == T || error("MemoryView element type is different from delimiter type")
    DelimitedIterator{T}(v, d)
end

Base.IteratorSize(::Type{<:DelimitedIterator}) = Base.SizeUnknown()
Base.eltype(::Type{DelimitedIterator{T}}) where {T} = ImmutableMemoryView{T}

function Base.iterate(d::DelimitedIterator, state::Int=1)
    len = length(d.v)
    if state > len
        return if state > len + 1 || iszero(len)
            nothing
        else
            v = @inbounds d.v[len : len - 1]
            (v, state + 1)
        end
    end
    next = findnext(==(d.d), d.v, state)
    if isnothing(next)
        (@inbounds d.v[state:len], len + 2)
    else
        (@inbounds d.v[state:(next - 1)], next + 1)
    end
end
