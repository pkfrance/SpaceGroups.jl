using LinearAlgebra

"""
    smith_normal_form(A::AbstractMatrix{T}) where T <: Integer

Returns (S, D, T) such that S*A*T = D, where D is the Smith Normal Form
of A, and S, T are unimodular matrices.
"""
function smith_normal_form(A::AbstractMatrix{T}) where T <: Integer
    m, n = size(A)
    D = copy(A)
    S = Matrix{T}(I, m, m)
    T_mat = Matrix{T}(I, n, n)
    
    # Work on the k-th diagonal element
    for k in 1:min(m, n)
        # 1. Pivot selection and zeroing out row/column
        while true
            # Find the smallest non-zero element in the remaining block
            idx = find_min_nonzero(D, k)
            idx === nothing && break # Entire block is zero
            
            # Swap to move pivot to (k, k)
            swap_rows!(D, S, k, idx[1])
            swap_cols!(D, T_mat, k, idx[2])
            
            # Eliminate other entries in column k
            changed = false
            for i in k+1:m
                if D[i, k] != 0
                    eliminate_row!(D, S, k, i)
                    changed = true
                end
            end
            
            # Eliminate other entries in row k
            for j in k+1:n
                if D[k, j] != 0
                    eliminate_col!(D, T_mat, k, j)
                    changed = true
                end
            end
            
            # If no changes were made in row/col elimination, 
            # we check if D[k,k] divides all other elements in the block.
            if !changed
                if check_and_fix_divisibility!(D, S, T_mat, k)
                    continue # Repeat reduction if divisibility was fixed
                end
                break
            end
        end
        
        # Ensure pivot is non-negative
        if D[k, k] < 0
            D[k, k] *= -1
            S[k, :] .*= -1
        end
    end
    
    return S, D, T_mat
end

function find_min_nonzero(D, k)
    m, n = size(D)
    best_val = nothing
    best_idx = nothing
    
    for i in k:m, j in k:n
        val = abs(D[i, j])
        if val > 0 && (best_val === nothing || val < best_val)
            best_val = val
            best_idx = (i, j)
        end
    end
    return best_idx
end

function swap_rows!(D, S, r1, r2)
    r1 == r2 && return
    D[r1, :], D[r2, :] = D[r2, :], D[r1, :]
    S[r1, :], S[r2, :] = S[r2, :], S[r1, :]
end

function swap_cols!(D, T, c1, c2)
    c1 == c2 && return
    D[:, c1], D[:, c2] = D[:, c2], D[:, c1]
    T[:, c1], T[:, c2] = T[:, c2], T[:, c1]
end

function eliminate_row!(D, S, k, i)
    a, b = D[k, k], D[i, k]
    g, x, y = gcdx(a, b)
    # Unimodular transformation matrix: [x y; -b/g a/g]
    u, v = div(a, g), div(b, g)
    
    row_k = D[k, :]
    row_i = D[i, :]
    D[k, :] = x * row_k + y * row_i
    D[i, :] = -v * row_k + u * row_i
    
    s_k = S[k, :]
    s_i = S[i, :]
    S[k, :] = x * s_k + y * s_i
    S[i, :] = -v * s_k + u * s_i
end

function eliminate_col!(D, T, k, j)
    a, b = D[k, k], D[k, j]
    g, x, y = gcdx(a, b)
    u, v = div(a, g), div(b, g)
    
    col_k = D[:, k]
    col_j = D[:, j]
    D[:, k] = x * col_k + y * col_j
    D[:, j] = -v * col_k + u * col_j
    
    t_k = T[:, k]
    t_j = T[:, j]
    T[:, k] = x * t_k + y * t_j
    T[:, j] = -v * t_k + u * t_j
end

function check_and_fix_divisibility!(D, S, T, k)
    m, n = size(D)
    for i in k+1:m, j in k+1:n
        if D[i, j] % D[k, k] != 0
            # Add row i to row k to introduce a non-divisible element in row k
            D[k, :] .+= D[i, :]
            S[k, :] .+= S[i, :]
            return true
        end
    end
    return false
end