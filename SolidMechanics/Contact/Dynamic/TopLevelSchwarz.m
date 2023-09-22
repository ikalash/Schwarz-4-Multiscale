function [element_properties, nodal_fields] = ...
        TopLevelSchwarz(top_level_params, element_properties, nodal_fields)


num_domains = top_level_params.num_domains;
num_steps = top_level_params.num_steps;
overlap_domains = top_level_params.overlap_domains;
num_elements_domain = top_level_params.num_elements_domain;
time_interval = top_level_params.time_interval;
rel_tol_domain = top_level_params.rel_tol_domain;
rel_tol_schwarz = top_level_params.rel_tol_schwarz;
abs_tol_domain = top_level_params.abs_tol_domain;
abs_tol_schwarz = top_level_params.abs_tol_schwarz;
max_iter_domain = top_level_params.max_iter_domain;
max_iter_schwarz = top_level_params.max_iter_schwarz;
min_iter_schwarz = top_level_params.min_iter_schwarz;
domain_interval = top_level_params.domain_interval;
schwarz_interval = top_level_params.schwarz_interval;
step_interval = top_level_params.step_interval;
time_steps = top_level_params.time_steps;
time_step = top_level_params.time_step;
integration_schemes = top_level_params.integration_schemes;


alpha1 = top_level_params.alpha1; %for Robin-Robin
alpha2 = top_level_params.alpha2; %for Robin-Robin

reference_positions = nodal_fields.reference_positions;
current_positions = nodal_fields.current_positions;
velocities = nodal_fields.velocities;
accelerations = nodal_fields.accelerations;
% a_int = nodal_fields.a_int;
reactions = nodal_fields.reactions;
posn_histories = nodal_fields.posn_histories;
velo_histories = nodal_fields.velo_histories;
acce_histories = nodal_fields.acce_histories;
reac_histories = nodal_fields.reac_histories;
cont_histories = nodal_fields.cont_histories;
schwarz_iter_history = nodal_fields.schwarz_iter_history;
end_masses = nodal_fields.end_masses;

moduli = element_properties.moduli;
areas = element_properties.areas;
densities = element_properties.densities;
internals = element_properties.internals;
forces = element_properties.forces;
strains = element_properties.strains;

initial_time = time_interval(1);
final_time = time_interval(2);

beta = top_level_params.beta;
gamma = top_level_params.gamma;
implicit = 0;
explicit = 1;

% Choose boundary conditons for Schwarz methood
% DN -> Dirichlet domain left, Neumann domain right; then, after each time step swap them
% ND -> Neumann domain left, Dirichlet domain right; then, after each time step swap them
% RR -> Robin-Robin 
boundary_cond = top_level_params.boundary_cond;
if boundary_cond == 'DN'
    dirichlet_bc = false(num_domains, 2);
    dirichlet_bc(:, 1) = true;
    neumann_bc = ~dirichlet_bc;
    robinrobin_bc = false;
end
if boundary_cond == 'ND'
    dirichlet_bc = true(num_domains, 2);
    dirichlet_bc(:, 1) = false;
    neumann_bc = ~dirichlet_bc;
    robinrobin_bc = false;
end
if boundary_cond == 'RR'
    robinrobin_bc = true;
    dirichlet_bc  = false(num_domains, 2);
    neumann_bc    = false(num_domains, 2);
end

% Estimate stable time step for each domain in case of explicit time
% integration
fmt_str = 'Warning domain %d, reduce time step from %0.8e to at most %0.8e.\n';
for domain = 1 : num_domains
  if integration_schemes(domain) == 1
    num_elements = num_elements_domain(domain);
    element_sizes = zeros(1, num_elements);
    X = reference_positions{domain};
    for element = 1 : num_elements
      element_sizes(element) = X(element + 1) - X(element);
    end
    wave_speeds = sqrt(moduli{domain} ./ densities{domain});
    stable_steps = element_sizes ./ wave_speeds;
    stable_step = 1.0 * min(stable_steps);
    tolerance = 1e-15;
    if time_steps(domain) - stable_step > tolerance
      fprintf(fmt_str, domain, time_steps(domain), stable_step);
      error('Program terminated.');
    end
  end
end
top_level_params.time_steps = time_steps;
proposed_time_step = max(time_steps);
if proposed_time_step > time_step
  fmt_str = 'Maximum domain step is: %0.8e, larger than global step %0.8e.\n';
  fprintf(fmt_str, proposed_time_step, time_step);
  error('Program terminated.');
end

potential_energies = zeros(num_domains, num_steps);
kinetic_energies = zeros(num_domains, num_steps);
total_energies = zeros(num_domains, num_steps);
linear_momenta = zeros(num_domains, num_steps);

 
% Before the first time step, compute the kinetic and potential energy as
% well as the initial acceleration in each domain
for domain = 1 : num_domains
  X = reference_positions{domain}';
  x = current_positions{domain}';
  v = velocities{domain}';
  E = moduli{domain}';
  A = areas{domain}';
  rho = densities{domain}';
  Q = internals{domain}';

  [~, M, ~, f, ~, ~, ~, PE] = Assembly(X, x, E, A, ...
    rho, Q, 0, 0, 1, ...
    element_properties.constitutive_law, ...
    top_level_params.mass_damping_coefficient, ...
    top_level_params.stiffness_damping_coefficient);

  a = (- M \ f)';
  posn_histories{domain}(1, :) = x';
  velo_histories{domain}(1, :) = v';
  acce_histories{domain}(1, :) = a';
  accelerations{domain} = acce_histories{domain}(1,:);
  p = M * v;
  KE = 0.5 * v' * p;
  LM = sum(p);

  potential_energies(domain, 1) = PE;
  kinetic_energies(domain, 1) = KE;
  total_energies(domain, 1) = KE + PE;
  linear_momenta(domain, 1) = LM;

  % Fill next history slot so that there is initial data for interpolation
  % of BCs if needed
  posn_histories{domain}(2, :) = x';
  velo_histories{domain}(2, :) = v';
  acce_histories{domain}(2, :) = a';
  end_masses{domain}(1) = rho(1) * A(1) * 0.5 * (X(2) - X(1));
  end_masses{domain}(2) = rho(end) * A(end) * 0.5 * (X(end) - X(end - 1));
end

contact_left = false(num_domains, 1);
contact_right = false(num_domains, 1);

%
% Main loop
%
for step = 1 : num_steps

  time = initial_time + (step - 1) * time_step;
  display_step = mod(step - 1, step_interval) == 0 || step == num_steps;

  if display_step == 1
    fprintf('step: %d, time: %0.8e\n', step, time);
  end

  prev_internals = internals;

  % Determine if contact has occurred
  for domain = 1 : num_domains
    domain_left = overlap_domains(domain, 1);
    domain_right = overlap_domains(domain, 2);

    if domain_right > 0
      curr_posn_right = current_positions{domain_right}(1);
      curr_posn = current_positions{domain}(end);
      reac = reactions{domain}(end);
      overlap = curr_posn > curr_posn_right;
      compress = reac < 0.0; 
      persist = compress && contact_right(domain) == true;
      contact_right(domain) = overlap || persist;
      contact_left(domain_right) = contact_right(domain);
    end

    if domain_left > 0
      curr_posn_left = current_positions{domain_left}(end);
      curr_posn = current_positions{domain}(1);
      reac = reactions{domain}(1);
      overlap = curr_posn_left > curr_posn ;
      compress = reac > 0.0;
      persist = compress && contact_left(domain) == true;
      contact_left(domain) = overlap || persist;
      contact_right(domain_left) = contact_left(domain);
    end

    cont_histories{domain}(step, 1) = contact_left(domain);
    cont_histories{domain}(step, 2) = contact_right(domain);

  end % Contact

  % Schwarz iteration
  iter_schwarz = 0;
  while true


      

    prev_schwarz_positions = current_positions;
    internals = prev_internals;

    display_schwarz = ...
        mod(iter_schwarz + 1, schwarz_interval) == 0 && display_step;

    baseline_norms = zeros(num_domains, 1);
    difference_norms = zeros(num_domains, 1);


    % Go one domain at a time.
    for domain = 1 : num_domains

      x_prev = posn_histories{domain}(step, :)';
      v_prev = velo_histories{domain}(step, :)';
      a_prev = acce_histories{domain}(step, :)';

      % Adjust the domain time step to be an integer multiple
      % of the global time step.
      num_domain_per_global = round(time_step / time_steps(domain));
      if num_domain_per_global < 1
        num_domain_per_global = 1;
      end
      domain_time_step = time_step / num_domain_per_global;

      num_elements = num_elements_domain(domain);
      num_nodes = num_elements + 1;

      domain_left = overlap_domains(domain, 1);
      domain_right = overlap_domains(domain, 2);

      dirichlet_bc_left = dirichlet_bc(domain, 1);
      dirichlet_bc_right = dirichlet_bc(domain, 2);

      neumann_bc_left = neumann_bc(domain, 1);
      neumann_bc_right = neumann_bc(domain, 2);

      alpha1_left  = alpha1(domain, 1);
      alpha1_right = alpha1(domain, 2);
      alpha2_left  = alpha2(domain, 1);
      alpha2_right = alpha2(domain, 2);

      fext = zeros(num_nodes, 1);

      for domain_step = 1 : num_domain_per_global

        domain_time = time + domain_time_step * domain_step;
        free_dof = true(num_nodes, 1);

        % If contact has occurred determine the appropriate boundary conditions.
        if contact_right(domain) == true
          % This results in interpolation with previous history data or
          % pure extrapolation in case there is enough data to do so.
          T = initial_time + (step - 1 : 1 : step) * time_step;
          X = posn_histories{domain_right}(step : step + 1, 1);
          V = velo_histories{domain_right}(step : step + 1, 1);
          A = acce_histories{domain_right}(step : step + 1, 1);
          R = reac_histories{domain_right}(step : step + 1, 1);
          t = domain_time;
          xn = interp1(T, X, t, 'linear', 'extrap');
          vn = interp1(T, V, t, 'linear', 'extrap');
          an = interp1(T, A, t, 'linear', 'extrap');
          rn = interp1(T, R, t, 'linear', 'extrap');
          if dirichlet_bc_right == true
            current_positions{domain}(end) = xn;
            velocities{domain}(end) = vn;
            accelerations{domain}(end) = an;
            free_dof(end) = false;
          end
          if neumann_bc_right == true
            fext(end) = -rn;
          end
          if robinrobin_bc == true
            current_positions{domain}(end) = alpha1_right * xn;
            velocities{domain}(end) = alpha1_right * vn;
            accelerations{domain}(end) = alpha1_right * an;
            free_dof(end) = false;
            fext(end) = alpha2_right * (-rn);
          end  
        end

        if contact_left(domain) == true
          % This results in interpolation with previous history data or
          % pure extrapolation in case there is enough data to do so.
          T = initial_time + (step - 1 : 1 : step) * time_step;
          X = posn_histories{domain_left}(step : step + 1, end);
          V = velo_histories{domain_left}(step : step + 1, end);
          A = acce_histories{domain_left}(step : step + 1, end);
          R = reac_histories{domain_left}(step : step + 1, end);
          t = domain_time;
          xn = interp1(T, X, t, 'linear', 'extrap');
          vn = interp1(T, V, t, 'linear', 'extrap');
          an = interp1(T, A, t, 'linear', 'extrap');
          rn = interp1(T, R, t, 'linear', 'extrap');
          if dirichlet_bc_left == true
            current_positions{domain}(1) = xn;
            velocities{domain}(1) = vn;
            accelerations{domain}(1) = an;
            free_dof(1) = false;
          end
          if neumann_bc_left == true
            fext(1) = -rn;
          end
          if robinrobin_bc == true
            current_positions{domain}(1) = alpha1_left * xn;
            velocities{domain}(1) = alpha1_left * vn;
            accelerations{domain}(1) = alpha1_left * an;
            free_dof(1) = false;
            fext(1) = alpha2_left * (-rn);
         end  
        end

        E = moduli{domain}';
        A = areas{domain}';
        rho = densities{domain}';
        Q = internals{domain}';

        X = reference_positions{domain}';
        x = current_positions{domain}';
        v = velocities{domain}';
        a = accelerations{domain}';

        display_domain = mod(domain_step, domain_interval) == 0 ...
            && display_schwarz;

        fmt_str = 'domain: %d, domain step: %d\n';

        if display_domain == 1
          fprintf(fmt_str, domain, domain_step);
        end

        % Integrate implicitly or explicitly
        dt = domain_time_step;
        if integration_schemes(domain) == implicit
          % Implicit integration

          % Predictor     
          % Classical Newmark         

          x_pre = x;
          x_pre(free_dof) = x_prev(free_dof) + ...
                            v_prev(free_dof) * dt + ...
                            a_prev(free_dof) * (0.5 - beta) * dt * dt;

          v_pre = v;
          v_pre(free_dof) = v_prev(free_dof) + ...
                            a_prev(free_dof) * (1 - gamma) * dt;

          a(free_dof) = (x(free_dof) - x_pre(free_dof)) / beta / dt / dt;

          iter_domain = 1;

          % Nonlinear solver iteration
          while true

            q = Q;

            [K, M, C, fres, q, P, e, PE] = Assembly(X, x, ...
              E, A, rho, q, 0, 0, 1, ...
              element_properties.constitutive_law, ...
              top_level_params.mass_damping_coefficient, ...
              top_level_params.stiffness_damping_coefficient);

            TM = M / beta / dt / dt + C * gamma / beta / dt + K;
            fvis = C * v;
            fine = M * a;
            fint = fres + fvis;
            r = fext - fint - fine ;
            delta = TM(free_dof, free_dof) \ r(free_dof);
            x(free_dof) = x(free_dof) + delta;
            a(free_dof) = (x(free_dof) - x_pre(free_dof)) / beta / dt / dt;
            v(free_dof) = v_pre(free_dof) + gamma * dt * a(free_dof); 

            p = M * v;
            KE = 0.5 * v' * p;
            LM = sum(p);

            abs_error = norm(delta);
            norm_posn = norm(x);

            if norm_posn > 0.0
              rel_error = abs_error / norm_posn;
            else
              if abs_error > 0
                rel_error = 1.0;
              else
                rel_error = 0.0;
              end
            end

            converged_rel = rel_error < rel_tol_domain;
            converged_abs = abs_error < abs_tol_domain;
            converged = (converged_rel || converged_abs);

            fmt_str = 'IMPL scheme nl iter: %d, rel error: %0.8e, abs error: %0.8e\n';

            if display_domain == 1
              fprintf(fmt_str, iter_domain, rel_error, abs_error);
            end

            if converged == 1 || iter_domain == max_iter_domain
              if display_domain == 0 && display_schwarz == 1
                fprintf(fmt_str, iter_domain, rel_error, abs_error);
              end
              break;
            end

            iter_domain = iter_domain + 1;

          end % iterate nonlinear solver

        elseif integration_schemes(domain) == explicit
          % Explicit integration

          % Predictor       
          % Central finite difference scheme
          x(free_dof) = x_prev(free_dof) + ...
                        v_prev(free_dof) * dt + ...
                        a_prev(free_dof) * dt * dt * 0.5 ;

          v(free_dof) = v_prev(free_dof) + ...
                        a_prev(free_dof) * (1 - gamma) * dt;

          [~, M, C, fres, q, P, e, PE] = Assembly(X, x, E, A, ...
          rho, Q, 1, 1, 1, ...
          element_properties.constitutive_law, ...
          top_level_params.mass_damping_coefficient, ...
          top_level_params.stiffness_damping_coefficient);
   
          fvis = diag(C) .* v;
          fine = M .* a;
          fint = fres + fvis;
          r = fext - fint - fine;
          delta = r(free_dof) ./ M(free_dof);
          a(free_dof) = a(free_dof) + delta;     

          % Corrector
          v(free_dof) = v(free_dof) + gamma * dt * a(free_dof);

          p = M .* v;
          KE = 0.5 * v' * p;
          LM = sum(p);

        else

          error('Unknown integration scheme');

        end

        % Update previous state
        x_prev = x;
        v_prev = v;
        a_prev = a;

      end % domain steps

      current_positions{domain} = x';
      velocities{domain} = v';
      accelerations{domain} = a';
      reactions{domain} = fint';

      internals{domain} = q';
      forces{domain} = P';
      strains{domain} = e';

      % Advanced one global step, record histories
      posn_histories{domain}(step + 1, :) = current_positions{domain};
      velo_histories{domain}(step + 1, :) = velocities{domain};
      acce_histories{domain}(step + 1, :) = accelerations{domain};
      reac_histories{domain}(step + 1, :) = reactions{domain};

      % Fill next history slot so that there is initial data for interpolation
      % of BCs if needed
      if (step < num_steps)
        posn_histories{domain}(step + 2, :) = current_positions{domain};
        velo_histories{domain}(step + 2, :) = velocities{domain};
        acce_histories{domain}(step + 2, :) = accelerations{domain};
        reac_histories{domain}(step + 2, :) = reactions{domain};
      end

      potential_energies(domain, step + 1) = PE;
      kinetic_energies(domain, step + 1) = KE;
      total_energies(domain, step + 1) = KE + PE;
      linear_momenta(domain, step + 1) = LM;

      baseline_norms(domain) = norm(x);
      dx = x' - prev_schwarz_positions{domain};
      difference_norms(domain) = norm(dx);

    end % for each domain

    if all(~contact_left) && all(~contact_right)
      if display_schwarz == 1
        fprintf('No contact, no Schwarz iterations\n');
      end
      break;
    end

    iter_schwarz = iter_schwarz + 1;

    baseline_norm = norm(baseline_norms);
    abs_schwarz_error = norm(difference_norms);

    if baseline_norm > 0.0
      rel_schwarz_error = abs_schwarz_error / baseline_norm;
    else
      if abs_schwarz_error > 0
        rel_schwarz_error = 1.0;
      else
        rel_schwarz_error = 0.0;
      end
    end

    converged_rel = rel_schwarz_error < rel_tol_schwarz;
    converged_abs = abs_schwarz_error < abs_tol_schwarz;
    converged = (converged_rel || converged_abs) && ...
                (iter_schwarz >= min_iter_schwarz);

    fmt_str = 'Schwarz iter: %d, rel error: %0.8e, abs error: %0.8e\n';
    if display_schwarz == 1
      fprintf(fmt_str, iter_schwarz, rel_schwarz_error, abs_schwarz_error);
    end

    if converged == 1 || iter_schwarz >= max_iter_schwarz
      if display_schwarz == 0 && display_step == 1
        fprintf(fmt_str, iter_schwarz, rel_schwarz_error, abs_schwarz_error);
      end
      break;
    end

  end % iterate domains (Schwarz iteration)

  if robinrobin_bc == false
    swap = dirichlet_bc;
    dirichlet_bc = neumann_bc;
    neumann_bc = swap;
  end  

  schwarz_iter_history(step) = iter_schwarz;

end % global time step

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

nodal_fields.tV = linspace(initial_time, final_time, num_steps + 1);

nodal_fields.potential_energies = potential_energies;
nodal_fields.kinetic_energies = kinetic_energies;
nodal_fields.total_energies = total_energies;
nodal_fields.linear_momenta = linear_momenta;

element_properties.internals = internals;
element_properties.forces = forces;
element_properties.strains = strains;
