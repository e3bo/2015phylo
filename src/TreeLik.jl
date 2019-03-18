module TreeLik
using ApproxFun
using LinearAlgebra

export r, λ, ψ, μ, τ

τ = 10
a,b = 0,τ
x = Fun(identity, a..b)

d = domain(x)
D = Derivative(d)
λ = 0.5
ψ = 1
μ = 0.2
r = λ + ψ + μ

L = D + r * I
u = [Evaluation(0); L] \ [1, 0]
dLdψ = Conversion(Chebyshev(0..τ), Ultraspherical(1,0..τ))
dudψ = [Evaluation(0); L] \ (-[0 * Evaluation(0); dLdψ] * u)
dudmu = [Evaluation(0); L] \ [0, 1]
dudμ = [Evaluation(0); L] \(-[0 * Evaluation(0); dLdψ] * u + [0; 1])

u1 = [Evaluation(0); L] \ [1, μ]
du1dμ = [Evaluation(0); L] \(-[0 * Evaluation(0); dLdψ] * u1 + [0; 1])
du1dψ = [Evaluation(0); L] \ (-[0 * Evaluation(0); dLdψ] * u1)

L1 = L - λ * u1
dL1dψ = dLdψ - λ * du1dψ
u2 = [Evaluation(0); L1] \ [1, μ]
du2dψ = [Evaluation(0); L1] \ (-[0 * Evaluation(0); dL1dψ] * u2)
dL1dμ = dLdψ - λ * du1dμ
du2dμ = [Evaluation(0); L1] \ (-[0 * Evaluation(0); dL1dμ] * u2 + [0, 1])
dL1dλ = dL1dψ - u1
du2dλ = [Evaluation(0); L1] \ (-[0 * Evaluation(0); dL1dλ] * u2)

end
