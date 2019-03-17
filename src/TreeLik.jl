module TreeLik
using ApproxFun
using LinearAlgebra

export r
export mu
export tau

tau = 10
a,b = 0,tau
x = Fun(identity, a..b)



d = domain(x)
D = Derivative(d)
r = 2
L = D + r * I

u = [Evaluation(0); L] \ [1, 0]

dLdr = Conversion(Chebyshev(0..tau), Ultraspherical(1,0..tau))
dudr = [Evaluation(0); L] \ (-[0 * Evaluation(0); dLdr] * u)

mu = 0.2

u1 = [Evaluation(0); L] \ [1, mu]

du1dr = [Evaluation(0); L] \ (-[0 * Evaluation(0); dLdr] * u1)




end
