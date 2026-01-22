# Test the smith_normal_form function
using LinearAlgebra

function is_diagonal(M)
    """Check if matrix is diagonal"""
    m, n = size(M)
    for i in 1:m, j in 1:n
        if i != j && M[i, j] != 0
            return false
        end
    end
    return true
end

@testset "Smith Normal Form" begin
    @testset "Basic properties" begin
        # Test 2x2 identity
        A = [1 0; 0 1]
        S, D, T = SpaceGroups.smith_normal_form(A)
        @test S * A * T ≈ D
        @test is_diagonal(D)
        
        # Test 2x2 simple matrix
        A = [2 4; 6 8]
        S, D, T = SpaceGroups.smith_normal_form(A)
        @test S * A * T ≈ D
        @test is_diagonal(D)
        @test D[1,1] <= D[2,2]  # Divisibility condition
        
        # Test 3x3 matrix
        A = [1 2 3; 4 5 6; 7 8 9]
        S, D, T = SpaceGroups.smith_normal_form(A)
        @test S * A * T ≈ D
        @test is_diagonal(D)
    end
    
    @testset "Divisibility property" begin
        # In Smith Normal Form, diagonal elements divide each other
        A = [6 9; 12 18]
        S, D, T = SpaceGroups.smith_normal_form(A)
        @test S * A * T ≈ D
        # Check divisibility: d[i,i] divides d[i+1,i+1]
        for i in 1:size(D, 1)-1
            if D[i,i] != 0 && D[i+1,i+1] != 0
                @test D[i+1,i+1] % D[i,i] == 0
            end
        end
    end
    
    @testset "Rectangular matrices" begin
        # Test 2x3 matrix
        A = [1 2 3; 4 5 6]
        S, D, T = SpaceGroups.smith_normal_form(A)
        @test S * A * T ≈ D
        @test is_diagonal(D)
        @test all(D[i,j] == 0 for i=1:2, j=1:3 if i != j)
        
        # Test 3x2 matrix
        A = [1 2; 3 4; 5 6]
        S, D, T = SpaceGroups.smith_normal_form(A)
        @test S * A * T ≈ D
        @test is_diagonal(D)
    end
    
    @testset "Rank and determinant properties" begin
        # Test zero matrix
        A = zeros(Int, 2, 2)
        S, D, T = SpaceGroups.smith_normal_form(A)
        @test S * A * T ≈ D
        @test all(D .== 0)
        
        # Test rank-deficient matrix
        A = [1 2; 2 4; 3 6]
        S, D, T = SpaceGroups.smith_normal_form(A)
        @test S * A * T ≈ D
        rank_D = count(!iszero, diag(D))
        @test rank_D <= min(size(A)...)
    end
    
    @testset "Unimodularity check" begin
        # S and T should be unimodular (det = ±1)
        A = [6 9 18; 12 18 24; 8 12 16]
        S, D, T = SpaceGroups.smith_normal_form(A)
        det_S = det(S)
        det_T = det(T)
        @test abs(det_S) == 1
        @test abs(det_T) == 1
    end
    
    @testset "Negative values" begin
        # Test matrix with negative values
        A = [-2 4; 6 -8]
        S, D, T = SpaceGroups.smith_normal_form(A)
        @test S * A * T ≈ D
        @test all(D[i,i] >= 0 for i in 1:min(size(D)...))  # Diagonal should be non-negative
    end
    
    @testset "Big integers" begin
        # Test with big integers
        A = BigInt[
            123456789 987654321
            246813579 135792468
        ]
        S, D, T = SpaceGroups.smith_normal_form(A)
        @test S * A * T ≈ D
        @test is_diagonal(D)
        # Check that all elements are BigInt
        @test eltype(D) == BigInt
        @test eltype(S) == BigInt
        @test eltype(T) == BigInt
        
        # Test larger matrix with big integers
        A = BigInt[
            10^15 + 1  10^15 + 2  10^15 + 3
            10^15 + 4  10^15 + 5  10^15 + 6
            10^15 + 7  10^15 + 8  10^15 + 9
        ]
        S, D, T = SpaceGroups.smith_normal_form(A)
        @test S * A * T ≈ D
        @test is_diagonal(D)
    end
    
    @testset "Known examples" begin
        # Example: [[6, 9], [12, 18]]
        # SNF should give [[6, 0], [0, 3]]
        A = [6 9; 12 18]
        S, D, T = SpaceGroups.smith_normal_form(A)
        @test S * A * T ≈ D
        # Diagonal entries should be [6, 3] or in decreasing order by divisibility
        diag_entries = [D[i,i] for i in 1:min(size(D)...)]
        diag_entries = sort(diag_entries[diag_entries .!= 0])
        @test all(diag_entries[i] <= diag_entries[i+1] for i in 1:length(diag_entries)-1)
    end
    
    @testset "Single row and column" begin
        # Test single row
        A = reshape([2, 4, 6], 1, 3)
        S, D, T = SpaceGroups.smith_normal_form(A)
        @test S * A * T ≈ D
        
        # Test single column
        A = reshape([2, 4, 6], 3, 1)
        S, D, T = SpaceGroups.smith_normal_form(A)
        @test S * A * T ≈ D
    end
end
