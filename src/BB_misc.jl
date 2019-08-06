"""
function find_all( A )
Return a tuple `(I, J, V)` where `I` and `J` are the row and column indexes of all the
values in matrix `A`, and `V` is a vector of all the values.
```jldoctest
julia> A = [ 1 2 ; 0 3 ; 0 4 ]
3×3 Array{Int64,2}:
 1  2
 0  3
 0  4
julia> find_all( A )
([1, 2, 3, 1, 2, 3], [1, 1, 1, 2, 2, 2], [1, 0, 0, 2, 3, 4])
```
"""
function find_all(A::AbstractMatrix{T}) where {T}
# adaptet from findnz()
    nnzA = length(A)
    I = zeros(Int, nnzA)
    J = zeros(Int, nnzA)
    NZs = Array{T,1}(nnzA)
    count = 1
    if nnzA > 0
        for j=indices(A,2), i=indices(A,1)
            Aij = A[i,j]
            I[count] = i
            J[count] = j
            NZs[count] = Aij
            count += 1
        end
    end
    return (I, J, NZs)
end

"""
    find_non_typemin( A )
Return a tuple `(I, J, V)` where `I` and `J` are the row and column indexes of all the
values in matrix `A`, and `V` is a vector of all the values, except for triples in which the value is typemin of the type of A.
```jldoctest
julia> A = [ typemin(Int64) 2 ; 0 3 ; 0 4 ]
3×3 Array{Int64,2}:
 -9223372036854775808  2
                    0  3
                    0  4
julia> find_non_typemin( A )
([2, 3, 1, 2, 3], [1, 1, 2, 2, 2], [0, 0, 2, 3, 4])
```
"""
function find_non_typemin(A::AbstractMatrix{T}) where {T}
# adaptet from findnz()
    nnzA = length(A)
    I = zeros(Int, nnzA)
    J = zeros(Int, nnzA)
    NZs = Array{T,1}(nnzA)
    if nnzA > 0
        count = 1
        for j=indices(A,2), i=indices(A,1)
            Aij = A[i,j]
            if Aij != typemin( T )
                I[count] = i
                J[count] = j
                NZs[count] = Aij
                count += 1
            end
        end
        resize!( I , count - 1 )
        resize!( J , count - 1 )
        resize!( NZs , count - 1 )
    end
    return (I, J, NZs)
end
"""
    find_non_typemax( A )
Return a tuple `(I, J, V)` where `I` and `J` are the row and column indexes of all the
values in matrix `A`, and `V` is a vector of all the values, except for triples in which the value is typemax of the type of A.
```jldoctest
julia> A = [ typemax(Int64) 2 ; 0 3 ; 0 4 ]
3×3 Array{Int64,2}:
  9223372036854775807  2
                    0  3
                    0  4
julia> find_non_typemax( A )
([2, 3, 1, 2, 3], [1, 1, 2, 2, 2], [0, 0, 2, 3, 4])
```
"""
function find_non_typemax(A::AbstractMatrix{T}) where {T}
# adaptet from findnz()
    nnzA = length(A)
    I = zeros(Int, nnzA)
    J = zeros(Int, nnzA)
    NZs = Array{T,1}(nnzA)
    if nnzA > 0
        count = 1
        for j=indices(A,2), i=indices(A,1)
            Aij = A[i,j]
            if Aij != typemax( T )
                I[count] = i
                J[count] = j
                NZs[count] = Aij
                count += 1
            end
        end
        resize!( I , count - 1 )
        resize!( J , count - 1 )
        resize!( NZs , count - 1 )
    end
    return (I, J, NZs)
end
