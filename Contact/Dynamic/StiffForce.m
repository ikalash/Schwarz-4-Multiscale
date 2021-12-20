function [element] = StiffForce(element)

X = element.reference_positions;
x = element.current_positions;

L = X(2) - X(1);
l = x(2) - x(1);

lambda = l / L;
strain = lambda - 1;

%[element] = element.constitutive_law(element);

% Update energy, stress and Hessian.
E = element.modulus;
A = element.area;
W = 0.5 * E * strain * strain;
P = E * strain;
C = E;

rho = element.density;
lumped = element.lumped;
diagonal = element.diagonal;

% Use midpoint rule to integrate.
%[~, DNa] = Lagrange1D2(0.0);
B = [-1.0; 1.0] ./ L;
f = L * B * P * A;
K = L * B * C * A * B';

if diagonal == 1
  M = [0.5; 0.5] * L * A * rho;
else
  if lumped == 1
    M = L * A * rho * [0.5, 0; 0, 0.5];
  else
    M = L * A * rho * [2.0, 1.0; 1.0, 2.0] / 6.0;
  end
end

if diagonal == 1
  element_mass_matrix = full(diag(M));
else
  element_mass_matrix = M;
end

element.energy_density = W;
element.stress = P;
element.stretch = lambda;
element.hessian = C;
element.mass = M;
element.stiffness = K;
element.damping = ...
    element.stiffness * element.stiffness_damping_coefficient + ...
    element_mass_matrix * element.mass_damping_coefficient;

element.internal_force = f;
element.stored_energy = L * A * W;
