function [output] = SchwarzMethod(settings)

top_level_params.overlap_domains = [0, 2; 1, 0];
top_level_params.theta = [0, 1.0; 1.0, 0];

top_level_params.num_domains = 2;
top_level_params.bar_area = settings.A;

dt = settings.dt;

if isfield(settings, 'dt_domains')
   % check for separate time step sizes for each domain and the global time step
   top_level_params.time_step = settings.dt_schwarz;
   top_level_params.time_steps = cell2mat(settings.dt_domains);
   dt = settings.dt_schwarz;
else
   top_level_params.time_step = settings.dt;
   top_level_params.time_steps = [settings.dt, settings.dt];
end

%top_level_params.time_step = settings.dt;
%top_level_params.time_steps = [settings.dt, settings.dt];

top_level_params.time_interval = [settings.t_start, settings.t_end];

time_difference = settings.t_end - settings.t_start;


num_steps = round(time_difference / dt);
top_level_params.num_steps = num_steps;

top_level_params.num_elements_domain = cell2mat(settings.nx);

top_level_params.initial_velocities = cell2mat(settings.v0);

top_level_params.zeta_target = settings.zeta_target;
top_level_params.stiffness_damping_contribution = settings.stiffness_damping_contribution;
top_level_params.mass_damping_contribution = 1 - settings.stiffness_damping_contribution;

top_level_params.mass_damping_coefficient = settings.a0;
top_level_params.stiffness_damping_coefficient = settings.a1;

top_level_params.limits_domains = zeros(2, 2);
top_level_params.limits_domains(:, 1) = cell2mat(settings.x0_start);
top_level_params.limits_domains(:, 2) = cell2mat(settings.x1_start);

top_level_params.gamma = settings.gamma;
top_level_params.beta = settings.beta;

top_level_params.bar_area = settings.A;
top_level_params.material_density = settings.rho;
top_level_params.youngs_modulus = settings.E;

top_level_params.constitutive_law = @(element)Constitutive_Linear_Elastic(element);


% convert settings to top_level_params
top_level_params.rel_tol_domain = settings.rel_tol_domain;
top_level_params.rel_tol_schwarz = settings.rel_tol_schwarz;
top_level_params.abs_tol_domain = settings.abs_tol_domain;
top_level_params.abs_tol_schwarz = settings.abs_tol_schwarz;
top_level_params.min_iter_schwarz = settings.min_iter_schwarz;
top_level_params.max_iter_schwarz = settings.max_iter_schwarz;
top_level_params.max_iter_domain = settings.max_iter_domain;
top_level_params.step_interval = settings.step_interval;
top_level_params.schwarz_interval = settings.schwarz_interval;
top_level_params.domain_interval = settings.domain_interval;

top_level_params.integration_schemes = cell2mat(settings.integration_schemes);


[element_properties, nodal_fields] = SetupSchwarz(top_level_params);

[element_properties, nodal_fields] = TopLevelSchwarz(top_level_params, element_properties, nodal_fields);

% convert nodal fields to
output.XV = nodal_fields.posn_histories;
output.VV = nodal_fields.velo_histories;
output.AV = nodal_fields.acce_histories;

for idomain = 1:top_level_params.num_domains
   output.XV{idomain} = output.XV{idomain}';
   output.VV{idomain} = output.VV{idomain}';
   output.AV{idomain} = output.AV{idomain}';
end

output.contact_force = -nodal_fields.reac_histories{1}(:, end);
contact_1 = nodal_fields.cont_histories{1};
output.contact = contact_1(:, 1) | contact_1(:, 2);
output.num_schwarz_iterations = nodal_fields.schwarz_iter_history;

output.tV = nodal_fields.tV;

% calculate energies
[output] = compute_kinetic_energy(output, settings);

end
