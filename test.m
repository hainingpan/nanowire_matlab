[L,U,p] = lu(A,'vector');
Afun = @(x) U\(L\(x(p)));
d = eigs(Afun,100,6,'smallestabs')