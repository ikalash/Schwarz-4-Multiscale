function [element] = Constitutive_Linear_Elastic(element)

E = element.modulus;
A = element.area;
l = element.stretch;
strain = l - 1;
q = element.internal;
update = element.update;

% Update energy, stress and Hessian.
W = 0.5 * E * strain * strain;
P = E * strain;
C = E;

element.energy_density = W;
element.stress = P;
element.hessian = C;
element.internal = q;
