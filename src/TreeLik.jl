module TreeLik
using ApproxFun
using LinearAlgebra

a,b = 0,10
x = Fun(identity, a..b)

d = domain(x)
D = Derivative(d)
r = 2
L = D + r * I

u = [Evaluation(0); L] \ [1, 0]

dLdr = Conversion(Chebyshev(0..10), Ultraspherical(1,0..10))
dudr = [Evaluation(0); L] \ (-[0 * Evaluation(0); dLdr] * u)

end
