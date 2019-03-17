module TreeLikTst
include("TreeLik.jl")
import .TreeLik
using Test
using ApproxFun
using LinearAlgebra

uexact = Fun(x->exp(-2x),  0..10)
exptest = @test norm(TreeLik.u - uexact, Inf) < 1e-14
dexptest = @test TreeLik.dudr(10)  ≈ -10 * exp(-2 * 10)

end
