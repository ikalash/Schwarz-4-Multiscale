clear variables;
initial_time = -0.0e-03;
final_time = 2.0e-04;
time_step = 1.0e-06;

output_filename = 'schwarz_method.mat';

top_level_params.time_interval = [initial_time, final_time];
% newmark-beta method parameters
% -> explicit central difference when gamma = 0.5 & beta = 0.0
% -> Average constant acceleration when gamma = 0.5 & beta = 0.25
% ->-> Average constant acceleration is unconditionally stable
top_level_params.gamma = 0.5;
top_level_params.beta = 0.25;

%bar properties
top_level_params.bar_area = 0.0254*0.0254;

%material properties
top_level_params.constitutive_law = ...
@(element)Constitutive_Linear_Elastic(element);
top_level_params.youngs_modulus = 206.8e09;
top_level_params.material_density = 7844.07;


% These depend on the number of domains
num_domains = 2;
switch num_domains
case 1
  top_level_params.num_domains = num_domains;
  top_level_params.integration_schemes = 0;
  top_level_params.time_steps = 1e-7;
  top_level_params.num_elements_domain = 100;
  top_level_params.limits_domains = [0, 1.0];
  top_level_params.overlap_domains = [0, 0];
  top_level_params.theta = [0, 0];
  top_level_params.initial_velocities = 0.0;
case 2
  offset = 2.54e-04;
  offsets = [-1, -1; 1, 1] * offset;
  limits = [-0.254, 0; 0, 0.254];
  top_level_params.num_domains = num_domains;
  top_level_params.integration_schemes = [0, 0];
  top_level_params.time_steps = [1.0e-06, 1.0e-06];
  top_level_params.num_elements_domain = [20, 20]';
  top_level_params.limits_domains = limits + offsets;
  top_level_params.overlap_domains = [0, 2; 1, 0];
  top_level_params.theta = [0, 1.0; 1.0, 0];
  top_level_params.initial_velocities = [5.13588, -5.13588];
  A = top_level_params.bar_area;
  rho = top_level_params.material_density;

  L = 2*0.254;

  M = rho * A * L;

  % JonStart
  % proportional damping coefficients
  % [C] = a0 * M + a1 * K
  % a0 -> mass proportional damping coefficient
  % a1 -> stiffness proportional damping coefficient
  % damping ratio zeta_i = (a0 + a1 * omega_i ^ 2) / (2 * omega_i)
  % if a1 = 0 -> a0 =  zeta_i * 2 * omega_i
  % if a0 = 0  -> zeta_i * 2 / omega_i
  % omega_i -> ith natural frequency of vibration
  top_level_params.zeta_target = 0.0;
  top_level_params.stiffness_damping_contribution = 1;


  top_level_params.mass_damping_contribution = 1 - top_level_params.stiffness_damping_contribution;

  K = top_level_params.bar_area * top_level_params.youngs_modulus / L;
  omega_i_est = sqrt(K / M);
  top_level_params.mass_damping_coefficient = ...
    top_level_params.zeta_target * 2.0 * omega_i_est * top_level_params.mass_damping_contribution;

  top_level_params.stiffness_damping_coefficient = ...
    top_level_params.zeta_target * 2.0 / omega_i_est * top_level_params.stiffness_damping_contribution;
  % JonEnd

  V = 5.13588;
  reference_energy = 0.5 * M * V * V;
  reference_momentum = M * V;
  zero_reference_momentum = true;

case 3
  offset = 2.0e-02;
  offsets = [-1, -1; 0, 0; 1, 1] * offset;
  limits = [-3/8, -1/8; -1/8, 1/8; 1/8, 3/8];
  top_level_params.num_domains = num_domains;
  top_level_params.integration_schemes = [0, 0, 0];
  top_level_params.time_steps = [1.0e-07, 1.0e-07, 1.0e-07];
  top_level_params.num_elements_domain = [200, 200, 200]';
  top_level_params.limits_domains = limits + offsets;
  top_level_params.overlap_domains = [0, 2; 1, 3; 2, 0];
  top_level_params.theta = [0, 1; 1, 1; 1, 0];
  top_level_params.initial_velocities = [100.0, 0.0, -100.0];
  A = top_level_params.bar_area;
  rho = top_level_params.material_density;
  L = 0.50;
  M = rho * A * L;
  V = 100.0;
  reference_energy = 0.5 * M * V * V;
  reference_momentum = M * V;
  zero_reference_momentum = true;
otherwise


end
time_interval = top_level_params.time_interval;
time_steps = top_level_params.time_steps;
time_difference = final_time - initial_time;
num_steps = round(time_difference ./ time_step);
top_level_params.num_steps = num_steps;
top_level_params.time_step = time_step;

fprintf('Initial number of steps: %d\n', num_steps);

top_level_params.rel_tol_domain = 1e-9;
top_level_params.rel_tol_schwarz = 1e-9;
top_level_params.abs_tol_domain = 1e-12;
top_level_params.abs_tol_schwarz = 1e-12;

top_level_params.min_iter_schwarz = 1;
top_level_params.max_iter_schwarz = 8;
top_level_params.max_iter_domain = 100000;

top_level_params.step_interval = 10;
top_level_params.schwarz_interval = 10;
top_level_params.domain_interval = 1;

[element_properties, nodal_fields] =  SetupSchwarz(top_level_params);

[element_properties, nodal_fields] = ...
TopLevelSchwarz(top_level_params, element_properties, nodal_fields);

% Plotting of results
num_plots = 10;
num_domains = top_level_params.num_domains;
num_elements_domain = top_level_params.num_elements_domain;
num_steps = top_level_params.num_steps;
plot_str = ['r-';'g-';'b-';'k-'];
schemes = ['IMPLICIT'; 'EXPLICIT'];
integration_schemes = top_level_params.integration_schemes;
plot_interval = round(num_steps / num_plots);
title_str = schemes(integration_schemes(1) + 1, :);
for domain = 2 : num_domains
  integration_scheme = integration_schemes(domain);
  title_str = strcat(title_str, '-', schemes(integration_scheme + 1, :));
end

close all;
figure(1);
subplot(3, 4, 1);
title('CONFIGURATION');
hold on;
for step = 1 : num_steps + 1
  if mod(step - 1, plot_interval) == 0
    for domain = 1 : num_domains
      x = nodal_fields.posn_histories{domain}(step, :);
      time = initial_time + (step - 1) * time_step;
      times = time * ones(size(x));
      plot(x, times, plot_str(domain,:), 'LineWidth', 0.7);
    end
  end
end
xlabel('POSITION [m]');
ylabel('TIME [s]');
grid on;
hold off;

subplot(3, 4, 2);
title('VELOCITIES');
hold on;
vm = min(top_level_params.initial_velocities);
vp = max(top_level_params.initial_velocities);
for step = 1 : num_steps + 1
  if mod(step - 1, plot_interval) == 0
    for domain = 1 : num_domains
      x = nodal_fields.posn_histories{domain}(step, :);
      v = nodal_fields.velo_histories{domain}(step, :);
      plot(x, v, plot_str(domain,:), 'LineWidth', 0.7);
    end
  end
end
xlabel('POSITION [m]');
ylabel('VELOCITY [m/s]');
grid on;
hold off;

KE = zeros(num_steps + 1, 1);
PE = zeros(num_steps + 1, 1);
TE = zeros(num_steps + 1, 1);
LM = zeros(num_steps + 1, 1);
times = initial_time + (0 : num_steps)' * time_difference / num_steps;
for step = 1 : num_steps + 1
  for domain = 1 : num_domains
    KE(step) = KE(step) + nodal_fields.kinetic_energies(domain, step);
    PE(step) = PE(step) + nodal_fields.potential_energies(domain, step);
    TE(step) = TE(step) + nodal_fields.total_energies(domain, step);
    LM(step) = LM(step) + nodal_fields.linear_momenta(domain, step);
  end
end

subplot(3, 4, 3);
hold on;
title('TOTAL ENERGY');
plot(times, KE, 'r-', 'LineWidth', 0.7);
plot(times, PE, 'g-', 'LineWidth', 0.7);
plot(times, TE, 'b-', 'LineWidth', 0.7);
legend({'KINETIC', 'STRAIN', 'TOTAL'}, 'Location', 'Southwest');
xlabel('TIME [s]');
ylabel('ENERGY [J]');
grid on;
hold off;

subplot(3, 4, 4);
hold on;
title('ENERGY ERROR');
EE = 100 * (TE - reference_energy) / reference_energy;
plot(times, EE, 'r-', 'LineWidth', 0.7);
xlabel('TIME [s]');
ylabel('RELATIVE ERROR [%]');
grid on;
hold off;

%subplot(3, 4, 5);
%hold on;
%title('TOTAL LINEAR MOMENTUM');
%plot(times, LM, 'r-', 'LineWidth', 0.7);
%xlabel('TIME');
%ylabel('MOMENTUM [kg m/s]');
%grid on;
%hold off;

subplot(3, 4, 5);
hold on;
title('LINEAR MOMENTUM ERROR');
if zero_reference_momentum == true
  ME = 100 * LM / reference_momentum;
else
  ME = 100 * (LM - reference_momentum) / reference_momentum;
end
plot(times, ME, 'r-', 'LineWidth', 0.7);
xlabel('TIME');
ylabel('RELATIVE ERROR [%]');
grid on;
hold off;

X = zeros(num_steps + 1, num_domains, 2);
V = zeros(num_steps + 1, num_domains, 2);
A = zeros(num_steps + 1, num_domains, 2);
R = zeros(num_steps + 1, num_domains, 2);
times = initial_time + (0 : num_steps)' * time_difference / num_steps;

for step = 1 : num_steps + 1
  for domain = 1 : num_domains
    if top_level_params.overlap_domains(domain, 1) > 0
      X(step, domain, 1) = nodal_fields.posn_histories{domain}(step, 1);
      V(step, domain, 1) = nodal_fields.velo_histories{domain}(step, 1);
      A(step, domain, 1) = nodal_fields.acce_histories{domain}(step, 1);
      R(step, domain, 1) = nodal_fields.reac_histories{domain}(step, 1);
    end
    if top_level_params.overlap_domains(domain, 2) > 0
      X(step, domain, 2) = nodal_fields.posn_histories{domain}(step, end);
      V(step, domain, 2) = nodal_fields.velo_histories{domain}(step, end);
      A(step, domain, 2) = nodal_fields.acce_histories{domain}(step, end);
      R(step, domain, 2) = nodal_fields.reac_histories{domain}(step, end);
    end
  end
end

plot_str = ['r-';'g-';'b-';'k-'];

subplot(3, 4, 6);
hold on;
title('CONTACT POSITIONS');
for domain = 1 : num_domains
  if top_level_params.overlap_domains(domain, 1) > 0
    plot(times, X(:, domain, 1), plot_str(domain, :), 'LineWidth', 0.7)
  end
  if top_level_params.overlap_domains(domain, 2) > 0
    plot(times, X(:, domain, 2), plot_str(domain, :), 'LineWidth', 0.7)
  end
end
xlabel('TIME [s]');
ylabel('POSITION [m]');
grid on;
hold off;

subplot(3, 4, 7);
hold on;
title('CONTACT VELOCITIES');
for domain = 1 : num_domains
  if top_level_params.overlap_domains(domain, 1) > 0
    plot(times, V(:, domain, 1), plot_str(domain, :), 'LineWidth', 0.7)
  end
  if top_level_params.overlap_domains(domain, 2) > 0
    plot(times, V(:, domain, 2), plot_str(domain, :), 'LineWidth', 0.7)
  end
end
xlabel('TIME [s]');
ylabel('VELOCITY [m/s]');
grid on;
hold off;

subplot(3, 4, 8);
hold on;
title('CONTACT ACCELERATIONS');
for domain = 1 : num_domains
  if top_level_params.overlap_domains(domain, 1) > 0
    plot(times, A(:, domain, 1), plot_str(domain, :), 'LineWidth', 0.7)
  end
  if top_level_params.overlap_domains(domain, 2) > 0
    plot(times, A(:, domain, 2), plot_str(domain, :), 'LineWidth', 0.7)
  end
end
xlabel('TIME [s]');
ylabel('ACCELERATION [m/s^2]');
grid on;
hold off;

subplot(3, 4, 9);
hold on;
title('CONTACT FORCES');
for domain = 1 : num_domains
  if top_level_params.overlap_domains(domain, 1) > 0
    plot(times, R(:, domain, 1), plot_str(domain, :), 'LineWidth', 0.7)
  end
  if top_level_params.overlap_domains(domain, 2) > 0
    plot(times, R(:, domain, 2), plot_str(domain, :), 'LineWidth', 0.7)
  end
end
xlabel('TIME [s]');
ylabel('FORCE [N]');
grid on;
hold off;

CH = nodal_fields.cont_histories;

subplot(3, 4, 10);
hold on;
title('CONTACT LEFT');
for domain = 1 : num_domains
  plot(times, CH{domain}(:, 1), plot_str(domain, :), 'LineWidth', 0.7)
end
xlabel('TIME [s]');
ylabel('CONTACT');
grid on;
hold off;

subplot(3, 4, 11);
hold on;
title('CONTACT RIGHT');
for domain = 1 : num_domains
  plot(times, CH{domain}(:, 2), plot_str(domain, :), 'LineWidth', 0.7)
end
xlabel('TIME [s]');
ylabel('CONTACT');
grid on;
hold off;

subplot(3, 4, 12);
hold on;
title('SCHWARZ ITERATIONS');
plot(times, nodal_fields.schwarz_iter_history, 'r-', 'LineWidth', 0.7)
xlabel('TIME [s]');
ylabel('NUMBER OF ITERATIONS');
grid on;
hold off;

% Save the output
save(output_filename)
