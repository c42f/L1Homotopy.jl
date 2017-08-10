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


@testset "Random matrix special cases" begin
    A = [-1.04 -0.70 -0.82  0.58 -0.01;
          1.04  1.29 -0.57 -0.02 -1.02;
         -0.29 -0.02  0.35 -0.22  0.98;
         -1.56  1.13  0.49 -0.23 -0.26;
         -0.72  1.39 -0.28  0.75  0.07]

    x_in = [0.0, 0.0, 1.0, 1.0, 0.0]

    x_out = L1_homotopy(A, A*x_in, 0.001)

    @test x_in ≈ x_out
    # This test doesn't yet pass, because the maintenance of the active set
    # doesn't handle the case where γ⁺ ≈ γ⁻ to numerical precision.
    @test_broken (x_out .== 0) == (x_in .== 0)
end
