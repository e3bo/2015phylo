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
dLdr = Conversion(Chebyshev(0..τ), Ultraspherical(1,0..τ))
dudr = [Evaluation(0); L] \ (-[0 * Evaluation(0); dLdr] * u)
dudmu = [Evaluation(0); L] \ [0, 1]
dudmunew = [Evaluation(0); L] \(-[0 * Evaluation(0); dLdr] * u + [0; 1])

u1 = [Evaluation(0); L] \ [1, μ]
du1dr = [Evaluation(0); L] \ (-[0 * Evaluation(0); dLdr] * u1)

L1 = L - λ * u1
u2 = [Evaluation(0); L1] \ [1, μ]

du2dpsi = [Evaluation(0); L] \ (-[0 * Evaluation(0); dLdr] * u2)

du2dmu = [Evaluation(0); L1] \ [0, 1]



end
