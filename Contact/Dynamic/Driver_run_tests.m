function [nodal_fields, times, settings] = Driver_run_test(settings_filename)

%read settings from yaml file
addpath('YAMLMatlab_0.4.3/')
settings = ReadYaml(settings_filename);
disp(settings)

output_dir = "output/";
mkdir(output_dir);

% Plot comparison to conventional methods ?
plot_option = true;

initial_time = settings.t_start;
final_time = settings.t_end;
time_step = settings.dt;

output_filename = settings.output_filename;
output_figname = settings.output_figname;
output_figname2 = settings.output_figname2;

top_level_params.time_interval = [initial_time, final_time];
% newmark-beta method parameters
% -> explicit central difference when gamma = 0.5 & beta = 0.0
% -> Average constant acceleration when gamma = 0.5 & beta = 0.25
% ->-> Average constant acceleration is unconditionally stable
top_level_params.gamma = settings.gamma;
top_level_params.beta = settings.beta;

% Choose boundary conditons for Schwarz methood
% DN -> Dirichlet domain left, Neumann domain right, then alternate (swap)
% ND -> Neumann domain left, Dirichlet domain right, then alternate (swap)
% RR -> Robin-Robin 
top_level_params.boundary_cond = 'DN';

%bar properties
top_level_params.bar_area = 1.0e-06;

%material properties
top_level_params.constitutive_law = ...
@(element)Constitutive_Linear_Elastic(element);
top_level_params.youngs_modulus = 1.0e09;
top_level_params.material_density = 1000;


% These depend on the number of domains
num_domains = 2;
tol_zero_momentum = 1.0e-08;
switch num_domains
case 1
  top_level_params.num_domains = num_domains;
  top_level_params.integration_schemes = implicit;
  top_level_params.time_steps = 1e-7;
  top_level_params.num_elements_domain = 100;
  top_level_params.limits_domains = [0, 1.0];
  top_level_params.overlap_domains = [0, 0];
  top_level_params.initial_velocities = 0.0;
  A = top_level_params.bar_area;
  rho = top_level_params.material_density;
  L = top_level_params.limits_domains(2) - top_level_params.limits_domains(1);
  M = rho * L * A;
  V = top_level_params.initial_velocities;
  reference_energy = 0.5 * M * V * V;
  reference_momentum = M * V;
  zero_reference_momentum = abs(reference_momentum) <= tol_zero_momentum;
case 2
  offset = 2.0e-02;
  offsets = [-1, -1; 1, 1] * offset;
  top_level_params.num_domains = num_domains;
  %parameters from 'settings'
  limits = [settings.limits_domains{1,1}, settings.limits_domains{1,2}; 
            settings.limits_domains{2,1}, settings.limits_domains{2,2}];
  top_level_params.limits_domains = limits + offsets;
  top_level_params.integration_schemes = [settings.integration_schemes{1}, settings.integration_schemes{2}];
  top_level_params.time_steps = [settings.time_steps{1}, settings.time_steps{2}]; 
  top_level_params.num_elements_domain = [settings.number_elements_domain{1}, settings.number_elements_domain{2}]';
  top_level_params.overlap_domains = [settings.overlap_domains{1,1}, settings.overlap_domains{1,2}; 
                                      settings.overlap_domains{2,1}, settings.overlap_domains{2,2}];
  top_level_params.initial_velocities = [settings.initial_velocities{1}, settings.initial_velocities{2}];
  top_level_params.alpha1 = [0, 0.5; 0.5, 0]; %for dislacement (Robin-Robin BC)
  top_level_params.alpha2 = [0, 1.; 1., 0]; %for traction (Robin-Robin BC)
  top_level_params.mass_damping_coefficient = 0.0;
  top_level_params.stiffness_damping_coefficient = 0.0;

  A = top_level_params.bar_area;
  rho = top_level_params.material_density;

  L1 = top_level_params.limits_domains(1, 2) - top_level_params.limits_domains(1, 1);
  L2 = top_level_params.limits_domains(2, 2) - top_level_params.limits_domains(2, 1);
  M1 = rho * A * L1;
  M2 = rho * A * L2;
  V1 = top_level_params.initial_velocities(1);
  V2 = top_level_params.initial_velocities(2);
  reference_energy = 0.5 * M1 * V1 * V1 + 0.5 * M2 * V2 * V2;
  reference_momentum = M1 * V1 + M2 * V2;
  zero_reference_momentum = abs(reference_momentum) <= tol_zero_momentum;
case 3
  offset = 2.0e-02;
  offsets = [-1, -1; 0, 0; 1, 1] * offset;
  limits = [-3/8, -1/8; -1/8, 1/8; 1/8, 3/8];
  top_level_params.num_domains = num_domains;
  top_level_params.integration_schemes = [implicit, implicit, implicit];
  top_level_params.time_steps = [1.0e-07, 1.0e-07, 1.0e-07];
  top_level_params.num_elements_domain = [200, 200, 200]';
  top_level_params.limits_domains = limits + offsets;
  top_level_params.overlap_domains = [0, 2; 1, 3; 2, 0];
  top_level_params.initial_velocities = [100.0, 0.0, -100.0];
  top_level_params.mass_damping_coefficient = 0.0;
  top_level_params.stiffness_damping_coefficient = 0.0;
  A = top_level_params.bar_area;
  rho = top_level_params.material_density;
  L1 = top_level_params.limits_domains(1, 2) - top_level_params.limits_domains(1, 1);
  L2 = top_level_params.limits_domains(2, 2) - top_level_params.limits_domains(2, 1);
  L3 = top_level_params.limits_domains(3, 2) - top_level_params.limits_domains(3, 1);
  M1 = rho * A * L1;
  M2 = rho * A * L2;
  M3 = rho * A * L3;
  V1 = top_level_params.initial_velocities(1);
  V2 = top_level_params.initial_velocities(2);
  V3 = top_level_params.initial_velocities(3);
  reference_energy = 0.5 * M1 * V1 * V1 + 0.5 * M2 * V2 * V2 + 0.5 * M3 * V3 * V3;
  reference_momentum = M1 * V1 + M2 * V2 + M3 * V3;
  zero_reference_momentum = abs(reference_momentum) <= tol_zero_momentum;
otherwise
  error("There are only setups for 1 to 3 subdomains.");
end

time_interval = top_level_params.time_interval;
time_steps = top_level_params.time_steps;
time_difference = final_time - initial_time;
num_steps = round(time_difference ./ time_step);
top_level_params.num_steps = num_steps;
top_level_params.time_step = time_step;

fprintf('Initial number of steps: %d\n', num_steps);

top_level_params.rel_tol_domain = 1e-12;
top_level_params.rel_tol_schwarz = 1e-12;
top_level_params.abs_tol_domain = 1e-15;
top_level_params.abs_tol_schwarz = 1e-15;

top_level_params.min_iter_schwarz = 1;
top_level_params.max_iter_schwarz = 128;
top_level_params.max_iter_domain = 100000;

top_level_params.step_interval = 100;
top_level_params.schwarz_interval = 10;
top_level_params.domain_interval = 100;

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

KE  = zeros(num_steps + 1, 1);
KE1 = zeros(num_steps + 1, 1); %KE of the first domain
KE2 = zeros(num_steps + 1, 1); %KE of the second domain
PE  = zeros(num_steps + 1, 1);
PE1 = zeros(num_steps + 1, 1); %PE of the first domain
PE2 = zeros(num_steps + 1, 1); %PE of the second domain
TE  = zeros(num_steps + 1, 1); 
TE1 = zeros(num_steps + 1, 1); %TE of the first domain
TE2 = zeros(num_steps + 1, 1); %TE of the second domain
LM  = zeros(num_steps + 1, 1);

times = initial_time + (0 : num_steps)' * time_difference / num_steps;
for step = 1 : num_steps + 1
  for domain = 1 : num_domains
    KE(step) = KE(step) + nodal_fields.kinetic_energies(domain, step);
    PE(step) = PE(step) + nodal_fields.potential_energies(domain, step);
    TE(step) = TE(step) + nodal_fields.total_energies(domain, step);
    LM(step) = LM(step) + nodal_fields.linear_momenta(domain, step);
  end
  KE1(step) = KE1(step) + nodal_fields.kinetic_energies(1, step); %KE of the first domain
  PE1(step) = PE1(step) + nodal_fields.potential_energies(1, step); %PE of the first domain
  TE1(step) = TE1(step) + nodal_fields.total_energies(1, step); %TE of the first domain
  KE2(step) = KE2(step) + nodal_fields.kinetic_energies(2, step); 
  PE2(step) = PE2(step) + nodal_fields.potential_energies(2, step); 
  TE2(step) = TE2(step) + nodal_fields.total_energies(2, step); 
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

subplot(3, 4, 5);
hold on;
title('LINEAR MOMENTUM ERROR');
if zero_reference_momentum == true
  ME = LM;
  ylabel('DEVIATION FROM ZERO MOMENTUM');
else
  ME = 100 * (LM - reference_momentum) / reference_momentum;
  ylabel('RELATIVE ERROR [%]');
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
savefig(output_figname);
save(output_filename);

%%
h(3) = figure ;

num_plots = 6;
subplot(3, 2, 1);
hold on;
title('CONTACT VELOCITIES -- domain left');
for domain = 1 : 1
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

subplot(3, 2, 2);
hold on;
title('CONTACT VELOCITIES -- domain right');
for domain = 2 : 2
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

subplot(3, 2, 3);
hold on;
title('CONTACT ACCELERATIONS -- domain left');
for domain = 1 : 1
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

subplot(3, 2, 4);
hold on;
title('CONTACT ACCELERATIONS -- domain right');
for domain = 2 : 2
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

subplot(3, 2, 5);
hold on;
title('CONTACT FORCES -- domain left');
for domain = 1 : 1
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

subplot(3, 2, 6);
hold on;
title('CONTACT FORCES -- domain right');
for domain = 2 : 2
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

% save(output_figname2);
savefig(h(3),output_figname2);

%%
if plot_option == true

output     = open(settings.conventional);
%analytical = output.analytical_output; 
explicitLM = output.explicit_lagrange_output;
explicitP  = output.explicit_penalty_output;
implicitP  = output.implicit_penalty_output;
schwarz    = open(output_filename); 

num_plots = 8;
h(2) = figure;
 %figure(2);
  
subplot(2, 4, 1);
%hold on; 
%plot(analytical.tV, analytical.VV_1, 'color',[1, 0, 0], 'LineWidth', 1); 
hold on; 
plot(explicitLM.tV, explicitLM.VV{1}(end,:), 'k', 'LineWidth', 1); 
hold on; 
plot(explicitP.tV, explicitP.VV{1}(end,:),  'b', 'LineWidth', 1); 
hold on; 
plot(implicitP.tV, implicitP.VV{1}(end,:), 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);
hold on; 
plot(schwarz.times, schwarz.V(:,1,2), 'color', [0, 0.75, 0.75], 'LineWidth', 1, 'LineStyle', '--'); 
xlabel('time'); 
ylabel('Contact Point Velocity -- domain left'); 
title('Contact Point Velocity -- domain left');  
%legend('Analytical', 'Explicit LM', 'Explicit Penalty', 'Implicit Penalty', 'Schwarz', 'Location', 'Best'); 
legend('Explicit LM', 'Explicit Penalty', 'Implicit Penalty', 'Schwarz', 'Location', 'Best'); 
grid on;
hold off;

subplot(2, 4, 2);
%hold on; 
%plot(analytical.tV, analytical.VV_1, 'color',[1, 0, 0], 'LineWidth', 1); 
hold on; 
plot(explicitLM.tV, explicitLM.VV{2}(1,:), 'k', 'LineWidth', 1); 
hold on; 
plot(explicitP.tV, explicitP.VV{2}(1,:),  'b', 'LineWidth', 1); 
hold on; 
plot(implicitP.tV, implicitP.VV{2}(1,:), 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);
hold on; 
plot(schwarz.times, schwarz.V(:,2,1), 'color', [0, 0.75, 0.75], 'LineWidth', 1, 'LineStyle', '--'); 
xlabel('time'); 
ylabel('Contact Point Velocity -- domain right'); 
title('Contact Point Velocity -- domain right');  
%legend('Analytical', 'Explicit LM', 'Explicit Penalty', 'Implicit Penalty', 'Schwarz', 'Location', 'Best'); 
legend('Explicit LM', 'Explicit Penalty', 'Implicit Penalty', 'Schwarz', 'Location', 'Best'); 
grid on;
hold off;

subplot(2, 4, 3);
%hold on; 
%plot(analytical.tV, analytical.contact_force, 'color',[1, 0, 0], 'LineWidth', 1); 
hold on; 
plot(explicitLM.tV, explicitLM.contact_force, 'k', 'LineWidth', 1); 
hold on; 
plot(explicitP.tV, explicitP.contact_force,  'b', 'LineWidth', 1); 
hold on; 
plot(implicitP.tV, implicitP.contact_force, 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);
hold on; 
plot(schwarz.times, -1*schwarz.R(:, 1, 2), 'color', [0, 0.75, 0.75], 'LineWidth', 1, 'LineStyle', '--'); 
xlabel('time'); 
ylabel('Contact Point Force -- domain left'); 
title('Contact Point Force -- domain left');  
%legend('Analytical', 'Explicit LM', 'Explicit Penalty', 'Implicit Penalty', 'Schwarz', 'Location', 'Best'); 
legend('Explicit LM', 'Explicit Penalty', 'Implicit Penalty', 'Schwarz', 'Location', 'Best'); 
grid on;
hold off;

subplot(2, 4, 4);
%hold on; 
%plot(analytical.tV, analytical.contact_force, 'color',[1, 0, 0], 'LineWidth', 1); 
% hold on; 
% plot(explicitLM.tV, explicitLM.contact_force, 'k', 'LineWidth', 1); 
% hold on; 
% plot(explicitP.tV, explicitP.contact_force,  'b', 'LineWidth', 1); 
% hold on; 
% plot(implicitP.tV, implicitP.contact_force, 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);
hold on; 
plot(schwarz.times, schwarz.R(:, 2, 1), 'color', [0, 0.75, 0.75], 'LineWidth', 1, 'LineStyle', '--'); 
xlabel('time'); 
ylabel('Contact Point Force -- domain right'); 
title('Contact Point Force -- domain right');  
%legend('Analytical', 'Explicit LM', 'Explicit Penalty', 'Implicit Penalty', 'Schwarz', 'Location', 'Best'); 
legend( 'Schwarz', 'Location', 'Best'); 
grid on;
hold off;

subplot(2, 4, 5);
%hold on;
%plot(analytical.tV, analytical.total_energy, 'color',[1, 0, 0], 'LineWidth', 1); 
hold on;
plot(explicitLM.tV, (explicitLM.potential_energy{1} + explicitLM.kinetic_energy{1}),'k', 'LineWidth', 1); 
hold on; 
plot(explicitP.tV, (explicitP.potential_energy{1} + explicitP.kinetic_energy{1}), 'b', 'LineWidth', 1); 
hold on; 
plot(implicitP.tV, (implicitP.potential_energy{1} + implicitP.kinetic_energy{1}), 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);
hold on; 
plot(schwarz.times, schwarz.PE1 + schwarz.KE1, 'color', [0, 0.75, 0.75], 'LineWidth', 1, 'LineStyle', '--'); 
xlabel('time'); 
ylabel('Total energy -- domain left '); 
title('Total energy -- domain left');  
%legend('Analytical', 'Explicit LM', 'Explicit Penalty', 'Implicit Penalty', 'Schwarz', 'Location', 'Best'); 
legend('Explicit LM', 'Explicit Penalty', 'Implicit Penalty', 'Schwarz', 'Location', 'Best'); 
grid on;
hold off;

subplot(2, 4, 6);
%hold on;
%plot(analytical.tV, analytical.total_energy, 'color',[1, 0, 0], 'LineWidth', 1); 
hold on;
plot(explicitLM.tV, (explicitLM.potential_energy{2} + explicitLM.kinetic_energy{2}),'k', 'LineWidth', 1); 
hold on; 
plot(explicitP.tV, (explicitP.potential_energy{2} + explicitP.kinetic_energy{2}), 'b', 'LineWidth', 1); 
hold on; 
plot(implicitP.tV, (implicitP.potential_energy{2} + implicitP.kinetic_energy{2}), 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);
hold on; 
plot(schwarz.times, schwarz.PE2 + schwarz.KE2, 'color', [0, 0.75, 0.75], 'LineWidth', 1, 'LineStyle', '--'); 
xlabel('time'); 
ylabel('Total energy -- domain right'); 
title('Total energy -- domain right');  
%legend('Analytical', 'Explicit LM', 'Explicit Penalty', 'Implicit Penalty', 'Schwarz', 'Location', 'Best'); 
legend('Explicit LM', 'Explicit Penalty', 'Implicit Penalty', 'Schwarz', 'Location', 'Best'); 
grid on;
hold off;

subplot(2, 4, 7);
hold on; 
plot(explicitLM.tV, explicitLM.kinetic_energy{1}, 'k', 'LineWidth', 1); 
hold on; 
plot(explicitP.tV, explicitP.kinetic_energy{1},  'b', 'LineWidth', 1); 
hold on; 
plot(implicitP.tV, implicitP.kinetic_energy{1}, 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);
hold on; 
plot(schwarz.times, schwarz.KE1, 'color', [0, 0.75, 0.75], 'LineWidth', 1, 'LineStyle', '--'); 
xlabel('time'); 
ylabel('Kinetic energy -- domain left'); 
title('Kinetic energy -- domain left');  
legend('Explicit LM', 'Explicit Penalty', 'Implicit Penalty', 'Schwarz', 'Location', 'Best'); 
grid on;
hold off;

subplot(2, 4, 8);
hold on; 
plot(explicitLM.tV, explicitLM.kinetic_energy{2}, 'k', 'LineWidth', 1); 
hold on; 
plot(explicitP.tV, explicitP.kinetic_energy{2},  'b', 'LineWidth', 1); 
hold on; 
plot(implicitP.tV, implicitP.kinetic_energy{2}, 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);
hold on; 
plot(schwarz.times, schwarz.KE2, 'color', [0, 0.75, 0.75], 'LineWidth', 1, 'LineStyle', '--'); 
xlabel('time'); 
ylabel('Kinetic energy -- domain right'); 
title('Kinetic energy -- domain right');  
legend('Explicit LM', 'Explicit Penalty', 'Implicit Penalty', 'Schwarz', 'Location', 'Best'); 
grid on;
hold off;

savefig(h(2),settings.comp_fig);
%close(h);
end

