function [K, M, f, q, P, e, V] = ...
    AssembleUpdate(X, u, E, A, rho, q, ds, md, ...
                   domain_size, regularize, diagonal, lumped, update, constitutive)

if regularize == 1
  [q, qm] = Update(X, u, E, A, q, domain_size, update, constitutive);
  [K, M, f, q, P, e, V] = ...
      Assembly(X, u, E, A, rho, qm, ds, md, diagonal, lumped, 0, constitutive);
else
  [K, M, f, q, P, e, V] = ...
      Assembly(X, u, E, A, rho, q, ds, md, diagonal, lumped, update, constitutive);
end
