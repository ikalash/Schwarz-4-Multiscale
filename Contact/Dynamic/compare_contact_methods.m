close all
clear all
%% Read settings file
% Add the yaml library
addpath('../../YAMLMatlab_0.4.3/')

%set the config filename
%settings_filename = 'contact_settings.yaml';
settings_filename = 'contact_settings_carpenter.yaml';

output_dir = "output_low_speed/";
mkdir(output_dir);
output_filename = output_dir + "output.mat";

settings = ReadYaml(settings_filename);

%% run implicit penalty method
implicit_penalty_output = penalty_implicit(settings);

%% run lagrange multiplier method
%settings.augmented_lagrange_mode = 0;
%lagrange_output = TopLevelLagrange_Disp(settings);  

%% run augmented lagrange method
%settings.augmented_lagrange_mode = 1;
%augmented_lagrange_output = TopLevelLagrange_Disp(settings);

%% run explicit penalty method
settings.penalty_mode = 1;
explicit_penalty_output = penalty_explicit(settings);

%% run explicit lagrange method
settings.penalty_mode = 0;
explicit_lagrange_output = lagrange_explicit(settings);


%% run implicit implicit Schwarz method
settings.dt_schwarz = 1.0e-7;
settings.dt_domains = {1.0e-7, 1.0e-7};
settings.integration_schemes = {0, 0};
implicit_implicit_schwarz_output = SchwarzMethod(settings);

%% run explicit implicit schwarz method
settings.dt_schwarz = 2.0e-7;
settings.dt_domains = {(2.0e-7) / 3, (2.0e-7) / 2};
settings.integration_schemes = {1, 0};
explicit_implicit_schwarz_output = SchwarzMethod(settings); 

%% compare to analytical solution
g = (settings.x0_start{2} -  settings.x1_start{1}) / 2;
L = settings.x1_start{1} -  settings.x0_start{1};
v0 = settings.v0{1};
t_impact = settings.t_start + g / v0;
t_release = t_impact + 2 * L * sqrt(settings.rho / settings.E);

t_start = settings.t_start;
t_end = settings.t_end;
dt = settings.dt;

nt = int32(round((t_end - t_start) / dt + 1)); 
t_end_actual = round(double(nt - 1)) * dt + t_start;
tV = linspace(t_start, t_end_actual, nt);

contact_force = zeros(nt, 1);
contact = false(nt, 1);
linear_momenta = zeros(nt, 1);
linear_momenta(:) = v0 * settings.rho * settings.A * L;
contact_force(tV > t_impact & tV < t_release) = v0 * sqrt(settings.E * settings.rho) * settings.A;
contact(tV > t_impact & tV < t_release) = true;

XV_1 = zeros(nt, 1);
VV_1 = zeros(nt, 1);

VV_1(tV < t_impact) = v0;
VV_1(tV > t_impact & tV < t_release) = 0;
VV_1(tV > t_release) = -v0;

XV_1(tV < t_impact) = v0 * (tV(tV < t_impact) - t_start) - g;
XV_1(tV > t_impact & tV < t_release) = 0;
XV_1(tV > t_release) = -v0 * (tV(tV > t_release) - t_release);

mass_averaged_velocity = zeros(nt, 1);
mass_averaged_velocity(tV < t_impact) = v0;
mass_averaged_velocity(tV >= t_impact & tV < t_release) = v0 - v0 * sqrt(settings.E) / (L * sqrt(settings.rho)) * (tV(tV >= t_impact & tV < t_release) - t_impact);
mass_averaged_velocity(tV >= t_release) = -v0;

kinetic_energy = zeros(nt, 1);
total_energy = zeros(nt, 1);
t1 = t_impact + L * sqrt(settings.rho / settings.E);
ke_start = 0.5 * settings.rho * settings.A * L * v0 * v0;
kinetic_energy(tV <= t_impact) = ke_start;

total_energy(:) = ke_start;

%kinetic_energy(tV > t_impact & tV <= t1) = ke_start - 0.5 * sqrt(settings.rho * settings.E) * L * settings.A * v0 * v0 * (tV(tV > t_impact & tV <= t1) - t_impact); 
%kinetic_energy(tV > t1 & tV <= t_release) = 0.5 * sqrt(settings.rho * settings.E) * L * settings.A * v0 * v0 * (tV(tV > t1 & tV <= t_release) - t1);

kinetic_energy(tV > t_impact & tV <= t1) = ke_start - 0.5 * sqrt(settings.rho * settings.E) * settings.A * v0 * v0 * (tV(tV > t_impact & tV <= t1) - t_impact); 
kinetic_energy(tV > t1 & tV <= t_release) = 0.5 * sqrt(settings.rho * settings.E) * settings.A * v0 * v0 * (tV(tV > t1 & tV <= t_release) - t1);

kinetic_energy(tV >= t_release) = ke_start;

potential_energy = total_energy - kinetic_energy;

disp(['impact time = ', num2str(t_impact), ' (sec)']);
disp(['release time = ', num2str(t_release), ' (sec)']);

analytical_output.tV = tV;
analytical_output.contact = contact;
analytical_output.contact_force = contact_force;
analytical_output.kinetic_energy = kinetic_energy;
analytical_output.total_energy = total_energy;
analytical_output.potential_energy = potential_energy;
analytical_output.mass_averaged_velocity = mass_averaged_velocity;
analytical_output.XV_1 = XV_1;
analytical_output.VV_1 = VV_1;

%% plot the quantities of interest

%{

figure();
hold on 
plot(tV, contact, '-k');
plot(lagrange_output.tV, lagrange_output.contact, '-r');
plot(augmented_lagrange_output.tV, augmented_lagrange_output.contact, '-y');
plot(penalty_output.tV, penalty_output.contact, '-g');

plot(explicit_lagrange_output.tV, explicit_lagrange_output.contact, '-m');
plot(explicit_penalty_output.tV, explicit_penalty_output.contact, '-c');

plot(schwarz_output.tV, schwarz_output.contact, '-b');

ylabel('contact condition')
xlabel('time, t, sec')
legend('analytical','lagrange', 'augmented lagrange','penalty', 'explicit lagrange', 'explicit penalty', 'schwarz');
title('Contact');
savefig('output/contact');
close

figure();
hold on
plot(tV, contact_force, '-k');
plot(lagrange_output.tV, lagrange_output.contact_force, '-r');
plot(augmented_lagrange_output.tV, augmented_lagrange_output.contact_force, '-y');
plot(penalty_output.tV, penalty_output.contact_force, '-g');

plot(explicit_lagrange_output.tV, explicit_lagrange_output.contact_force, '-m');
plot(explicit_penalty_output.tV, explicit_penalty_output.contact_force, '-c');

plot(schwarz_output.tV, schwarz_output.contact_force, '-b');
xlabel('time, t, sec')
ylabel('contact force, \lambda, (N)')
legend('analytical','lagrange', 'augmented lagrange','penalty', 'explicit lagrange', 'explicit penalty', 'schwarz');
savefig('output/contact_force');
close

figure();
hold on
plot(tV, XV_1, '-k');
plot(lagrange_output.tV, lagrange_output.XV{1}(end, :), '-r');
plot(augmented_lagrange_output.tV, augmented_lagrange_output.XV{1}(end, :), '-y');
plot(penalty_output.tV, penalty_output.XV{1}(end, :), '-g');

plot(explicit_lagrange_output.tV, explicit_lagrange_output.XV{1}(end, :), '-m');
plot(explicit_penalty_output.tV, explicit_penalty_output.XV{1}(end, :), '-c');

plot(schwarz_output.tV, schwarz_output.XV{1}(end, :), '-b');
xlabel('time, t, sec');
ylabel('bar 1 right position, x, m');
legend('analytical','lagrange', 'augmented lagrange','penalty', 'explicit lagrange', 'explicit penalty', 'schwarz');
savefig('output/x1_right');
close

figure();
hold on
plot(tV, -XV_1, '-k');
plot(lagrange_output.tV, lagrange_output.XV{2}(1, :), '-r');
plot(augmented_lagrange_output.tV, augmented_lagrange_output.XV{2}(1, :), '-y');
plot(penalty_output.tV, penalty_output.XV{2}(1, :), '-g');

plot(explicit_lagrange_output.tV, explicit_lagrange_output.XV{2}(1, :), '-m');
plot(explicit_penalty_output.tV, explicit_penalty_output.XV{2}(1, :), '-c');

plot(schwarz_output.tV, schwarz_output.XV{1}(end, :), '-b');
xlabel('time, t, sec');
ylabel('bar 2 left position, x, m');
legend('analytical','lagrange', 'augmented lagrange','penalty', 'explicit lagrange', 'explicit penalty', 'schwarz');
savefig('output/x2_left');
close

figure();
hold on
plot(tV, VV_1, '-k');
plot(lagrange_output.tV, lagrange_output.VV{1}(end, :), '-r');
plot(augmented_lagrange_output.tV, augmented_lagrange_output.VV{1}(end, :), '-y');
plot(penalty_output.tV, penalty_output.VV{1}(end, :), '-g');

plot(explicit_lagrange_output.tV, explicit_lagrange_output.VV{1}(end, :), '-m');
plot(explicit_penalty_output.tV, explicit_penalty_output.VV{1}(end, :), '-c');

plot(schwarz_output.tV, schwarz_output.VV{1}(end, :), '-b');
xlabel('time, t, sec');
ylabel('bar 1 right velocity, v, m/s');
legend('analytical','lagrange', 'augmented lagrange','penalty', 'explicit lagrange', 'explicit penalty', 'schwarz');
savefig('output/v1_right');
close

figure();
hold on
plot(tV, -VV_1, '-k');
plot(lagrange_output.tV, lagrange_output.VV{2}(1, :), '-r');
plot(augmented_lagrange_output.tV, augmented_lagrange_output.VV{2}(1, :), '-y');
plot(penalty_output.tV, penalty_output.VV{2}(1, :), '-g');

plot(explicit_lagrange_output.tV, explicit_lagrange_output.VV{2}(1, :), '-m');
plot(explicit_penalty_output.tV, explicit_penalty_output.VV{2}(1, :), '-c');

plot(schwarz_output.tV, schwarz_output.VV{2}(1, :), '-b');
xlabel('time, t, sec');
ylabel('bar 2 left velocity, v, m/s');
legend('analytical','lagrange', 'augmented lagrange','penalty', 'explicit lagrange', 'explicit penalty', 'schwarz');
savefig('output/v2_left');
close

figure();
hold on
plot(lagrange_output.tV, lagrange_output.AV{1}(end, :), '-r');
plot(augmented_lagrange_output.tV, augmented_lagrange_output.AV{1}(end, :), '-y');
plot(penalty_output.tV, penalty_output.AV{1}(end, :), '-g');

plot(explicit_lagrange_output.tV, explicit_lagrange_output.AV{1}(end, :), '-m');
plot(explicit_penalty_output.tV, explicit_penalty_output.AV{1}(end, :), '-c');

plot(schwarz_output.tV, schwarz_output.AV{1}(end, :), '-b');
xlabel('time, t, sec');
ylabel('bar 1 right right acceleration, a, m/s^2');
legend('lagrange', 'augmented lagrange','penalty', 'explicit lagrange', 'explicit penalty', 'schwarz');
savefig('output/a1_right');
close

figure();
hold on
plot(lagrange_output.tV, lagrange_output.AV{2}(1, :), '-r');
plot(augmented_lagrange_output.tV, lagrange_output.AV{2}(1, :), '-y');
plot(penalty_output.tV, penalty_output.AV{2}(1, :), '-g');

plot(explicit_lagrange_output.tV, explicit_lagrange_output.AV{2}(1, :), '-m');
plot(explicit_penalty_output.tV, explicit_penalty_output.AV{2}(1, :), '-c');

plot(schwarz_output.tV, schwarz_output.AV{2}(1, :), '-b');
xlabel('time, t, sec');
ylabel('bar 2 left acceleration, a, m/s^2');
legend('lagrange', 'augmented lagrange','penalty', 'explicit lagrange', 'explicit penalty', 'schwarz');
savefig('output/a2_left');
close

figure();
hold on
plot(tV, kinetic_energy, '-k')
plot(lagrange_output.tV, lagrange_output.kinetic_energy{1}, '-r');
plot(augmented_lagrange_output.tV, augmented_lagrange_output.kinetic_energy{1}, '-y');
plot(penalty_output.tV, penalty_output.kinetic_energy{1}, '-g');

plot(explicit_lagrange_output.tV, explicit_lagrange_output.kinetic_energy{1}, '-m');
plot(explicit_penalty_output.tV, explicit_penalty_output.kinetic_energy{1}, '-c');

plot(schwarz_output.tV, schwarz_output.kinetic_energy{1}, '-b');
xlabel('time, t, sec');
ylabel('bar 1 potential energy, PE, J');
legend('analytical','lagrange', 'augmented lagrange','penalty', 'explicit lagrange', 'explicit penalty', 'schwarz');
savefig('output/kinetic_energy');
close

figure();
hold on
plot(tV, potential_energy, '-k')
plot(lagrange_output.tV, lagrange_output.potential_energy{1}, '-r');
plot(augmented_lagrange_output.tV, augmented_lagrange_output.potential_energy{1}, '-y');
plot(penalty_output.tV, penalty_output.potential_energy{1}, '-g');

plot(explicit_lagrange_output.tV, explicit_lagrange_output.potential_energy{1}, '-m');
plot(explicit_penalty_output.tV, explicit_penalty_output.potential_energy{1}, '-c');

plot(schwarz_output.tV, schwarz_output.potential_energy{1}, '-b');
xlabel('time, t, sec');
ylabel('bar 1 potential energy, PE, J');
legend('analytical','lagrange', 'augmented lagrange','penalty', 'explicit lagrange', 'explicit penalty', 'schwarz');
savefig('output/potential_energy');
close

figure();
hold on
plot(tV, total_energy, '-k')
plot(lagrange_output.tV, lagrange_output.total_energy{1}, '-r');
plot(augmented_lagrange_output.tV, augmented_lagrange_output.total_energy{1}, '-y');
plot(penalty_output.tV, penalty_output.total_energy{1}, '-g');

plot(explicit_lagrange_output.tV, explicit_lagrange_output.total_energy{1}, '-m');
plot(explicit_penalty_output.tV, explicit_penalty_output.total_energy{1}, '-c');

plot(schwarz_output.tV, schwarz_output.total_energy{1}, '-b');
xlabel('time, t, sec');
ylabel('bar 1 total energy, KE, J');
legend('analytical','lagrange', 'augmented lagrange','penalty', 'explicit lagrange', 'explicit penalty', 'schwarz');
savefig('output/total_energy');
close

%figure();
%hold on
%lagrange_momenta_error = lagrange_output.linear_momenta / linear_momenta - 1.0;
%penalty_momenta_error = penalty_output.linear_momenta / linear_momenta - 1.0;
%plot(lagrange_output.tV, lagrange_momenta_error);
%plot(penalty_output.tV, penalty_momenta_error);
%xlabel('time, t, sec');
%ylabel('momentum error');
%legend('lagrange','penalty');
%savefig('output/bar1_momentum_error');
%close

figure();
hold on
lagrange_energy_error = lagrange_output.total_energy{1} ./ total_energy' - 1.0;
augmented_lagrange_energy_error = augmented_lagrange_output.total_energy{1} ./ total_energy' - 1.0;
penalty_energy_error = penalty_output.total_energy{1} ./ total_energy' - 1.0;

explicit_lagrange_energy_error = explicit_lagrange_output.total_energy{1} ./ total_energy' - 1.0;
explicit_penalty_energy_error = explicit_penalty_output.total_energy{1} ./ total_energy' - 1.0;

schwarz_energy_error = schwarz_output.total_energy{1} ./ total_energy' - 1.0;
plot(lagrange_output.tV, lagrange_energy_error, '-r');
plot(augmented_lagrange_output.tV, augmented_lagrange_energy_error, '-y');
plot(penalty_output.tV, penalty_energy_error, '-g');

plot(explicit_lagrange_output.tV, explicit_lagrange_energy_error, '-m');
plot(explicit_penalty_output.tV, explicit_penalty_energy_error, '-c');

plot(schwarz_output.tV, schwarz_energy_error, '-b');
xlabel('time, t, sec');
ylabel('energy error');
legend('lagrange', 'augmented_lagrange',  'penalty', 'explicit lagrange', 'explicit penalty', 'schwarz');
savefig('output/bar1_energy_error');
close


figure();
hold on
plot(tV, mass_averaged_velocity, '-k')
plot(lagrange_output.tV, lagrange_output.mass_averaged_velocity{1}, '-r');
plot(augmented_lagrange_output.tV, augmented_lagrange_output.mass_averaged_velocity{1}, '-y');
plot(penalty_output.tV, penalty_output.mass_averaged_velocity{1}, '-g');

plot(explicit_lagrange_output.tV, explicit_lagrange_output.mass_averaged_velocity{1}, '-m');
plot(explicit_penalty_output.tV, explicit_penalty_output.mass_averaged_velocity{1}, '-c');

plot(schwarz_output.tV, schwarz_output.mass_averaged_velocity{1}, '-b');
xlabel('time, t, sec');
ylabel('bar 1 mass_averaged_velocity, v, m/s');
legend('analytical', 'lagrange', 'augmented lagrange', 'penalty', 'explicit lagrange', 'explicit penalty', 'schwarz');
savefig('output/bar1_mass_averaged_velocity');
close

% open the figures and place into subfigure panel
%}
%% save the output

save(output_filename,'settings', 'analytical_output', 'implicit_penalty_output', 'explicit_lagrange_output','explicit_penalty_output','implicit_implicit_schwarz_output', 'explicit_implicit_schwarz_output'); 

