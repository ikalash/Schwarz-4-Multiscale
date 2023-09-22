function [K, M, C, f, q, P, e, V] = ...
    Assembly(X, x, E, A, rho, q, diagonal, lumped, update, constitutive, ...
    a0, a1)

num_nodes = length(X);
num_elements = num_nodes - 1;

K = zeros(num_nodes, num_nodes);

C = zeros(num_nodes, num_nodes);

if diagonal == 1
  M = zeros(num_nodes, 1);
else
  M = zeros(num_nodes, num_nodes);
end

f = zeros(num_nodes, 1);
P = zeros(num_elements, 1);
e = zeros(num_elements, 1);

V = 0;

for i = 1 : num_elements
  element.reference_positions = X(i : i + 1);
  element.current_positions = x(i : i + 1);
  element.internal = q(i);
  element.modulus = E(i);
  element.area = A(i);
  element.density = rho(i);
  element.diagonal= diagonal;
  element.lumped = lumped;
  element.update = update;
  element.constitutive_law = constitutive;

  element.mass_damping_coefficient = a0;
  element.stiffness_damping_coefficient = a1;

  [element] = StiffForce(element);

  Ki = element.stiffness;
  Mi = element.mass;
  Ci = element.damping;
  fi = element.internal_force;
  q(i) = element.internal;
  P(i) = element.stress;
  e(i) = element.stretch;

  K(i : i + 1, i : i + 1) = K(i : i + 1, i : i + 1) + Ki;
  C(i : i + 1, i : i + 1) = C(i : i + 1, i : i + 1) + Ci;

  if (diagonal == 1)
    M(i : i + 1) = M(i : i + 1) + Mi;
  else
    M(i : i + 1, i : i + 1) = M(i : i + 1, i : i + 1) + Mi;
  end

  f(i : i + 1) = f(i : i + 1) + fi;

  V = V + element.stored_energy;
end

K = sparse(K);
C = sparse(C);

if diagonal == 0
    M = sparse(M);
end
