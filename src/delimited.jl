"""
    DelimitedIterator(data, x::T)

Return an iterator over memory-backed data `data` of eltype `T`.
Returns `MemView`s of the same elements as `data`, split by `x`.
"""
struct DelimitedIterator{T}
    v::ImmutableMemView{T}
    d::T
end

function DelimitedIterator(x, d::T) where {T}
    v = ImmutableMemView(x)
    eltype(v) == T || error("MemView(x) must be of eltype T") # TODO
    DelimitedIterator{T}(v, d)
end

Base.IteratorSize(::Type{<:DelimitedIterator}) = Base.SizeUnknown()
Base.eltype(::Type{DelimitedIterator{T}}) where {T} = ImmutableMemView{T}

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
