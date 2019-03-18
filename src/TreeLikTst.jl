module TreeLikTst
include("TreeLik.jl")
using .TreeLik

using ApproxFun
using DifferentialEquations
using LinearAlgebra
using Test

normtol = 1e-10
uexact = Fun(x-> exp(-r * x), 0..τ)
u1exact = Fun(x->exp(-r * x) + (1 - exp( -r * x)) * μ / r, 0..τ)
u1wrong = Fun(x->exp(-r * x) + (1 - exp( -r * τ)) * μ / r, 0..τ)

pinfnorm(x) = maximum(abs(x)) # norm(x, Inf) is not working as of 2019 March 17

layer1tests = @testset "test likelihood and gradient for single layer" begin
    @test pinfnorm(TreeLik.u - uexact)  ≈ 0 atol = normtol
    @test TreeLik.dudr(τ)  ≈ -τ * uexact(τ) atol = normtol
    @test pinfnorm(TreeLik.u1 - u1exact) ≈ 0 atol = normtol
    @test pinfnorm(TreeLik.u1 - u1wrong) ≉ 0 atol = normtol
    @test TreeLik.u1(1) ≉ u1wrong(1) atol = normtol
    @test TreeLik.du1dr(τ) ≈ -τ * exp(-r * τ) +  μ / r^2 * (r * (τ * exp(-r * τ)) - (1 - exp(-r * τ))) atol = normtol
    @test TreeLik.dudmu(τ) ≈ (1 - exp(-r * τ)) / r atol = normtol
end;

using ParameterizedFunctions
using DifferentialEquations

f = @ode_def begin
  du1 = -(λ + ψ + μ) * u1 + μ
  du2 = -(λ + ψ + μ) * u2 + λ * u1 * u2 + μ 
end  λ ψ μ

u0 = [1.0;1.0]
tspan = (0.0,τ)
p = [λ, ψ, μ]
prob = ODEProblem(f,u0,tspan,p)

function test_f(p)
  _prob = remake(prob;u0=convert.(eltype(p),prob.u0),p=p)
  solve(_prob,Vern9(),abstol=1e-14,reltol=1e-14,save_everystep=false)[end]
end

using ForwardDiff, DiffEqSensitivity, Calculus
fd_res = ForwardDiff.jacobian(test_f,p)

calc_res = Calculus.finite_difference_jacobian(test_f, p)

#using DiffResults
#res = DiffResults.JacobianResult(u0,p) # Build the results object
#DiffResults.jacobian!(res,p) # Populate it with the results
#val = DiffResults.value(res) # This is the sol[end]
#jac = DiffResults.jacobian(res) 

sol = solve(prob)



layer2tests = @testset "test likelihood and gradient for two layers" begin
    @test TreeLikTst.TreeLik.u2.(sol.t) ≈ sol[2,:] rtol = 1e-4
end;

end
