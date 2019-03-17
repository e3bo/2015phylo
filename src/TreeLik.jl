module TreeLik
using ApproxFun
using LinearAlgebra

export u

a,b = 0,10
x = Fun(identity, a..b)

d = domain(x)
D = Derivative(d)
L = D + I

u = [Evaluation(0); L] \ [1, 0]

end
