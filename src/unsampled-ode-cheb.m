A = chebop(0, 5);
x = chebfun('x', [0 5]);
phi = @(u) diff(u) + 2 * u;
u = x;
nrmv = 1;
y = 0;
while nrmv > 1e-11
  A.op = @(v) 2*v;
  A.lbc = @(v) u(0) + v(0) -1;
  v = A\(-phi(u));
  nrmv = norm(v); y = y + norm(v);
  u = u + v;
end

addpath('chebfun')
L = chebop(0, 10);
L.op = @(x, u1, u2, u3) [diff(u1) + 2 * u1 ; diff(u2) + 2 * u2 - 0.4 * u1; diff(u3) + 2 * u3 - 0.4 * u2];
L.lbc = @(u1, u2, u3) [u1-1; u2 - 1; u3 - 1];
rhs = [.2; .2; .2];
U = L\rhs;

addpath('chebfun')
L = chebop(0, 10);
L.op = @(u1, u2, u3, u4) [diff(u1) + 2 * u1 ; diff(u2) + 2 * u2 - 0.4 * u1; diff(u3) + 2 * u3 - 0.4 * u2; diff(u4) + 2 * u4 - 0.4 * u3];
L.lbc = @(u1, u2, u3, u4) [u1-1; u2 - 1; u3 - 1, u4 - 1];
rhs = [.2; .2; .2; .2];
U = solvebvp(L, rhs);
