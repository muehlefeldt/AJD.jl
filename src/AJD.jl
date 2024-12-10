module AJD

function Is_Commuting(A::AbstractMatrix, B::AbstractMatrix)
    return A*B == B*A #check if the Matrices are commuting else the Diag won't work
end
function get_diag_elements(A::AbstractMatrix)
        
end
export Is_Commuting
end
