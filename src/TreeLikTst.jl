module TreeLikTst
include("TreeLik.jl")
using .TreeLik
using Test
using ApproxFun
using LinearAlgebra

normtol = 1e-10
layer1tests = @testset "test likelihood and gradient for single layer" begin
    uexact = Fun(x->exp(-r * x),  0..tau)
    exptest = @test norm(TreeLik.u - uexact, Inf)  ≈ 0 atol = normtol
    @test TreeLik.dudr(tau)  ≈ -tau * exp(-r * tau)
    u1exact = Fun(x->exp(-r * x) + (1 - exp( -r * tau)) * mu / r, 0..tau)
    @test norm(TreeLik.u1 - u1exact, Inf)  ≈ 0 atol = normtol
    @test TreeLik.du1dr(tau)  ≈ -tau * exp(-r * tau) + mu / r^2 * (r * (tau * exp(-r * tau)) - (1 - exp(-r * tau)))
    @test TreeLik.dudmu(tau)  ≈ (1 - exp(-r * tau)) / r
end;

end
