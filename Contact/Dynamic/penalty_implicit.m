function [output] = penalty_implicit(settings)
  % unpack the struct

  dt = settings.dt;
  t_end = settings.t_end;
  t_start = settings.t_start;
  nx = settings.nx;
  v0 = settings.v0;
  A = settings.A;
  E = settings.E;
  rho = settings.rho;
  a0 = settings.a0;
  a1 = settings.a1;
  x0_start = settings.x0_start;
  x1_start = settings.x1_start;

  beta = settings.beta;
  gamma = settings.gamma;

  max_iterations = settings.max_iterations;
  newmark_tolerance = settings.newmark_tolerance;
  lambda = settings.lambda;

  small_double = 1.0e-64;

  num_domains = length(nx);

  num_nodes = zeros(num_domains,1);
  for idomain = 1:num_domains
    num_nodes(idomain) = nx{idomain} + 1;
  end
  num_nodes_total = sum(num_nodes);
  left_domain_indices = zeros(num_domains,1);
  right_domain_indices = zeros(num_domains,1);

  for idomain = 1:num_domains
    if idomain == 1
      left_domain_indices(idomain) = 1;
      right_domain_indices(idomain) = nx{idomain} + 1;
    else
      left_domain_indices(idomain) = 1 + right_domain_indices(idomain - 1);
      right_domain_indices(idomain) = left_domain_indices(idomain) + nx{idomain};
    end
  end

  nt = int32(round((t_end - t_start) / dt + 1));
  t_end_actual = round(double(nt - 1)) * dt + t_start;
  tV = linspace(t_start, t_end_actual, nt);

  XV = cell(num_domains, 1);
  VV = cell(num_domains, 1);
  AV = cell(num_domains, 1);

  %F_extV = cell(num_domains, 1);
  %F_intV = cell(num_domains, 1);

  it_per_time_step = zeros(nt,1);
  error_per_time_step = zeros(nt,1);

  contact_force = zeros(nt,1);
  contact = false(nt, 1);

  for idomain = 1:num_domains
    %TV{idomain} = tV;

    xv_initial = ...
      linspace(x0_start{idomain}, x1_start{idomain}, nx{idomain} + 1);
    XV{idomain} = zeros(nx{idomain} + 1, nt);
    VV{idomain} = zeros(nx{idomain} + 1, nt);
    AV{idomain} = zeros(nx{idomain} + 1, nt);

    %F_intV{idomain} = zeros(nx(idomain), nt);
    %F_extV{idomain} = zeros(nx(idomain), nt);

    % set initial conditions
    XV{idomain}(:, 1) = xv_initial;
    VV{idomain}(:, 1) = v0{idomain};

  end

  f_int = zeros(num_nodes_total, 1);

  for idt = 1:(nt-1)

    M = zeros(num_nodes_total, num_nodes_total);
    C = zeros(num_nodes_total, num_nodes_total);
    K = zeros(num_nodes_total, num_nodes_total);
    % assemble the block matrices
    for idomain = 1:num_domains
      idx_start = left_domain_indices(idomain);
      idx_end = right_domain_indices(idomain);
      xV = XV{idomain}(:, idt);
      [Ki, Mi, Ci] = FEMAssembly(xV, A, E, rho, a0, a1);
      K(idx_start:idx_end, idx_start:idx_end) = Ki;
      M(idx_start:idx_end, idx_start:idx_end) = Mi;
      C(idx_start:idx_end, idx_start:idx_end) = Ci;
    end

    % make them sparse
    K = sparse(K);
    M = sparse(M);
    C = sparse(C);

    x = zeros(num_nodes_total, 1);
    a = zeros(num_nodes_total, 1);
    v = zeros(num_nodes_total, 1);

    for idomain = 1:num_domains
      idx_start = left_domain_indices(idomain);
      idx_end = right_domain_indices(idomain);
      a(idx_start:idx_end) = AV{idomain}(:, idt);
      v(idx_start:idx_end) = VV{idomain}(:, idt);
    end

    abs_error_tol = settings.newmark_tolerance;

    [a_new, v_new, x_new, it_count, abs_error, f_contact] = ...
      newmark_iterative(M, C, K, f_int, a, v, x, XV, abs_error_tol, max_iterations, beta, gamma, dt, lambda, idt);

    contact_force(idt) = f_contact;
    if abs(f_contact) - small_double > 0
      contact(idt) = true;
    end

    it_per_time_step(idt) = it_count;
    error_per_time_step(idt) = abs_error;

    % update internal force vector
    f_int = f_int + K * x_new;

    for idomain = 1:num_domains
      idx_start = left_domain_indices(idomain);
      idx_end = right_domain_indices(idomain);
      AV{idomain}(:, idt + 1) = a_new(idx_start:idx_end);
      VV{idomain}(:, idt + 1) = v_new(idx_start:idx_end);
      XV{idomain}(:, idt + 1) = XV{idomain}(:, idt) + x_new(idx_start:idx_end);
    end

  end

  output.XV = XV;
  output.VV = VV;
  output.AV = AV;
  output.tV = tV;

  output.contact = contact;
  output.contact_force = contact_force;

  output = compute_kinetic_energy(output, settings);

end