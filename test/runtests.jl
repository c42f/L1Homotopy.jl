using L1Homotopy
using Base.Test

# write your own tests here
@testset "Homotopy" begin
    N = 100
    A = randn(N,N)
    x_in = Float64.(rand(N) .> 0.9)
    x_out = L1_homotopy(A, A*x_in, 0.1)
    # Solution is close but approximate on the support
    @test x_in ≈ x_out
    # Solution should be exactly equal to zero off the support
    @test (x_out .== 0) == (x_in .== 0)
end
