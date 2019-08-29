syms cA k1A k2A I4 positive;

syms aj real;

M = [0,cos(aj),sin(aj)];

psi = k1A/(2*k2A)*(exp( k2A*(I4-1).^2 ) - 1)

dpsi = diff(psi,I4);

syms lr lt lz positive;

F = diag([lr lt lz]);
C = F'*F;

dpsi = subs(dpsi,I4,M*C*M');

s = simplify(2*F*(cA/2+dpsi*(M'*M))*F)