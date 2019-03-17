module TreeLikTst
include("TreeLik.jl")
using .TreeLik
using Test
using ApproxFun
using LinearAlgebra

uexact = Fun(x->exp(-2x),  0..tau)
exptest = @test norm(TreeLik.u - uexact, Inf) < 1e-14
dexptest = @test TreeLik.dudr(10)  ≈ -tau * exp(-2 * tau)

u1exact = Fun(x->exp(-2x) + 0.2, 0..10)
u1valtest = @test norm(TreeLik.u1 - u1exact, Inf) < 1e-14
du1drtest = @test TreeLik.du1dr(10)  ≈ -tau * exp(-r * tau) + mu / r^2 * (r * (tau * exp(-r * tau)) - (1 - exp(-r * tau))) 
dudmutest = @test TreeLik.dudmu(10)  ≈ (1 - exp(-r * tau)) / r

end
