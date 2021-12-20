function [element_properties, nodal_fields] = SetupSchwarz(top_level_params)

num_steps = top_level_params.num_steps;
num_domains = top_level_params.num_domains;

reference_positions = cell(num_domains);
current_positions = cell(num_domains);
velocities = cell(num_domains);
accelerations = cell(num_domains);
reactions = cell(num_domains);

limits_domains = top_level_params.limits_domains;
num_elements_domain = top_level_params.num_elements_domain;
initial_velocities = top_level_params.initial_velocities;

for domain = 1 : num_domains
  left = limits_domains(domain, 1);
  right = limits_domains(domain, 2);
  num_elements = num_elements_domain(domain);
  num_nodes = num_elements + 1;
  length = right - left;
  delta = length / num_elements;
  initial_velocity = initial_velocities(domain);
  reference_positions{domain} = (left : delta : right);
  current_positions{domain} = reference_positions{domain};
  velocities{domain} = initial_velocity * ones(1, num_nodes);
  accelerations{domain} = zeros(1, num_nodes);
  reactions{domain} = zeros(1, num_nodes);
  posn_histories{domain} = zeros(num_steps + 1, num_nodes);
  velo_histories{domain} = zeros(num_steps + 1, num_nodes);
  acce_histories{domain} = zeros(num_steps + 1, num_nodes);
  reac_histories{domain} = zeros(num_steps + 1, num_nodes);
  cont_histories{domain} = false(num_steps + 1, 2);
  schwarz_iter_history = zeros(num_steps + 1, 1);
end

nodal_fields.reference_positions = reference_positions;
nodal_fields.current_positions = current_positions;
nodal_fields.velocities = velocities;
nodal_fields.accelerations = accelerations;
nodal_fields.reactions = reactions;

nodal_fields.posn_histories = posn_histories;
nodal_fields.velo_histories = velo_histories;
nodal_fields.acce_histories = acce_histories;
nodal_fields.reac_histories = reac_histories;
nodal_fields.cont_histories = cont_histories;
nodal_fields.schwarz_iter_history = schwarz_iter_history;

element_properties.moduli = cell(num_domains);
element_properties.areas = cell(num_domains);
element_properties.densities = cell(num_domains);
element_properties.internals = cell(num_domains);
element_properties.forces = cell(num_domains);
element_properties.strains = cell(num_domains);

element_properties.mass_damping_coefficients = cell(num_domains);
element_properties.stiffness_damping_coefficients = cell(num_domains);

for domain = 1 : num_domains
  num_elements = num_elements_domain(domain);

  E = top_level_params.youngs_modulus * ones(1, num_elements);
  A = top_level_params.bar_area * ones(1, num_elements);
  rho = top_level_params.material_density * ones(1, num_elements);
  q = zeros(1, num_elements);

  element_properties.constitutive_law = top_level_params.constitutive_law;
  element_properties.moduli{domain} = E;
  element_properties.areas{domain} = A;
  element_properties.densities{domain} = rho;
  element_properties.internals{domain} = q;
  element_properties.forces{domain} = zeros(1, num_elements);
  element_properties.strains{domain} = zeros(1, num_elements);

  element_properties.mass_damping_coefficients{domain} = ...
    top_level_params.mass_damping_coefficient;

  element_properties.stiffness_damping_coefficients{domain} = ...
    top_level_params.stiffness_damping_coefficient;

end
