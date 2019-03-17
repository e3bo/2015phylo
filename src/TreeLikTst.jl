module TreeLikTst
include("TreeLik.jl")
import .TreeLik
using Test
using ApproxFun
using LinearAlgebra

uexact = Fun(x->exp(-x),  0..10)
@test norm(TreeLik.u - uexact) < 1e-14





end
