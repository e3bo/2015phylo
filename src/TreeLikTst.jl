module TreeLikTst
include("TreeLik.jl")
using .TreeLik
using Test
using ApproxFun
using LinearAlgebra

normtol = 1e-10
uexact = Fun(x-> exp(-r * x), 0..tau)
u1exact = Fun(x->exp(-r * x) + (1 - exp( -r * x)) * mu / r, 0..tau)
u1wrong = Fun(x->exp(-r * x) + (1 - exp( -r * tau)) * mu / r, 0..tau)

pinfnorm(x) = maximum(abs(x)) # norm(x, Inf) is not working as of 2019 March 17

layer1tests = @testset "test likelihood and gradient for single layer" begin
    @test pinfnorm(TreeLik.u - uexact)  ≈ 0 atol = normtol
    @test TreeLik.dudr(tau)  ≈ -tau * uexact(tau) atol = normtol
    @test pinfnorm(TreeLik.u1 - u1exact) ≈ 0 atol = normtol
    @test pinfnorm(TreeLik.u1 - u1wrong) ≉ 0 atol = normtol
    @test TreeLik.u1(1) ≉ u1wrong(1) atol = normtol
    @test TreeLik.du1dr(tau) ≈ -tau * exp(-r * tau) + mu / r^2 * (r * (tau * exp(-r * tau)) - (1 - exp(-r * tau))) atol = normtol
    @test TreeLik.dudmu(tau) ≈ (1 - exp(-r * tau)) / r atol = normtol
end;

end
