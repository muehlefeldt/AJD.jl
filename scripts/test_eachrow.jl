using BenchmarkTools
function get_diag_elements_eachrow1(A::AbstractArray)
    D = zeros(axes(A))
    #slower than 1:rows but "cleaner"
    row = 1
    for _ in eachrow(A[:,:,begin])
        D[row,row,:] = A[row,row,:]
        row += 1
    end
    return D
end
#slower than enumerate
function get_diag_elements_eachrow2(A::AbstractArray)
    D = zeros(axes(A))
    #slower than 1:rows but "cleaner"
    for row in eachrow(A[:,:,begin]).axes
        D[row,row,:] = A[row,row,:]
    end
    return D
end
#enumerate is slower than eachindex or eachrow1
function get_diag_elements_eachrow3(A::AbstractArray)
    D = zeros(axes(A))
    #slower than 1:rows but "cleaner"
    for (row,i) in enumerate(eachrow(A[:,:,begin]))
        D[row,row,:] = A[row,row,:]
    end
    return D
end
function get_diag_elements_eachindex(A::AbstractArray)
    D = zeros(axes(A))
    # bit faster than calling eachindex each iteration
    indexing = eachindex(A[:,begin,begin])
    #slower than 1:rows but "cleaner"
    for row in indexing
        D[row,row,:] = A[row,row,:]
    end
    return D
end
#pairs is the slowest way
function get_diag_elements_pairs(A::AbstractArray)
    D = zeros(axes(A))
    #slower than 1:rows but "cleaner"
    for (index,value) in pairs(IndexCartesian(),A[:,:,begin])
        if index[1] != index[2]
            D[index[1],index[1],:] = A[index[1],index[1],:]
        end
    end
    return D

end
#still the fastest way and takes the least allocations
function get_diag_elements(A::AbstractArray)
    rows,columns,k = size(A)
    D = zeros(rows, columns, k)
    #slower than 1:rows but "cleaner"
    for row in 1:rows
        D[row,row,:] = A[row,row,:]
    end
    return D
end
function get_diag_elements_axes(A::AbstractArray)
    
    D = zeros(axes(A))
    iterator = axes(A,1)
    for row in iterator
        D[row,row,:] = A[row,row,:]
    end
    
    return D
end
M = ones(3,3,3)

@btime get_diag_elements(M)
@btime get_diag_elements_eachindex(M)
@btime get_diag_elements_eachrow1(M)
@btime get_diag_elements_eachrow2(M)
@btime get_diag_elements_eachrow3(M)
@btime get_diag_elements_pairs(M)
