function [output] = TopLevelLagrange_Vel(settings)

    % This follows the approach of
    % LAGRANGE CONSTRAINTS FOR TRANSIENT FINITE ELEMENT SURFACE CONTACT
    % by N.J. Carpenter, R. L. Taylor, and M. G. Katona
    disp('TopLevelLagrange_Velocity method');

    beta_0 = 0.5;
    beta_1 = 0.5;

    % equations 4, 5, and 6

    % unpack the struct
    dt = settings.dt;
    t_end = settings.t_end;
    t_start = settings.t_start;
    nx = cell2mat(settings.nx);
    v0 = cell2mat(settings.v0);

    A = settings.A;
    E = settings.E;
    rho = settings.rho;
    zeta_target = settings.zeta_target;

    K_damping_contrib = settings.stiffness_damping_contribution;
    M_damping_contrib = 1 - K_damping_contrib;

    x0_start = cell2mat(settings.x0_start);
    x1_start = cell2mat(settings.x1_start);

    num_domains = length(nx);

    bar_length = abs(x1_start(1) - x0_start(1));
    bar_stiffness = A * E / bar_length;
    bar_mass = bar_length * A * rho;

    omega = sqrt(bar_stiffness / bar_mass);

    a0 = 2.0 * omega * zeta_target * M_damping_contrib;
    a1 = 2.0 / omega * zeta_target * K_damping_contrib;

    num_nodes = nx + 1;
    num_nodes_total = sum(num_nodes);

    left_domain_indices = zeros(num_domains,1);
    right_domain_indices = zeros(num_domains,1);

    left_domain_indices(1) = 1;
    right_domain_indices(1) = num_nodes(1);

    for idomain = 2:num_domains
      left_domain_indices(idomain) = right_domain_indices(idomain - 1) + 1;
      right_domain_indices(idomain) = left_domain_indices(idomain) + nx(idomain);
    end

    nt = int32(round((t_end - t_start) / dt + 1));
    t_end_actual = round(double(nt - 1)) * dt + t_start;
    tV = linspace(t_start, t_end_actual, nt);

    K = zeros(num_nodes_total);
    M = zeros(num_nodes_total);
    C = zeros(num_nodes_total);

    f_int = zeros(num_nodes_total, 1);
    f_ext = zeros(num_nodes_total, 1);


    u = zeros(num_nodes_total, 1);
    dudt = zeros(num_nodes_total, 1);
    d2udt2 = zeros(num_nodes_total, 1);

    num_contact_points = 1;

    contact_force = zeros(nt, 1);
    contact = zeros(nt,1);

    K_domains = cell(num_domains, 1);

    XV = cell(num_domains, 1);
    VV = cell(num_domains, 1);
    AV = cell(num_domains, 1);
    for idomain = 1:num_domains
      XV{idomain} = zeros(num_nodes(idomain), nt);
      VV{idomain} = zeros(num_nodes(idomain), nt);
      AV{idomain} = zeros(num_nodes(idomain), nt);

      % set initial conditions
      XV{idomain}(:, 1) = linspace(x0_start(idomain), x1_start(idomain), num_nodes(idomain));
      VV{idomain}(:, 1) = v0(idomain);

      idx_start = left_domain_indices(idomain);
      idx_end = right_domain_indices(idomain);
      dudt(idx_start:idx_end) = v0(idomain);

    end

    % Allocate room for contact constraints matrix
    G = zeros(1, num_nodes_total);
    G(num_nodes(idomain)) = -1;
    G(num_nodes(idomain) + 1) = 1;

    % calculate velocity at time step n- 1/2
    v_old = zeros(num_nodes_total, 1);

    for idt = 1:(nt-1)
      disp(idt);
      % assemble mass and stiffness for each domain
      for idomain = 1:num_domains
        [K_ii, M_ii, C_ii] = FEMAssembly(XV{idomain}(:, idt), A, E, rho, a0, a1);
        K_domains{idomain} = K_ii;
        idx_start = left_domain_indices(idomain);
        idx_end = right_domain_indices(idomain);
        K(idx_start:idx_end, idx_start:idx_end) = K_ii;
        M(idx_start:idx_end, idx_start:idx_end) = M_ii;
        C(idx_start:idx_end, idx_start:idx_end) = C_ii;
      end

      % determine if there is contact at time step n+1 using
      %domain1_right = XV{1}(end, idt) + VV{1}(end, idt) * dt;
      %domain2_left = XV{2}(1, idt) + VV{2}(1, idt) * dt;

      domain1_right = XV{1}(end, idt) + VV{1}(end, idt) * dt + 0.5 * AV{1}(end, idt) * dt^2;
      domain2_left = XV{2}(1, idt) + VV{2}(1, idt) * dt + 0.5 * AV{2}(1, idt) * dt^2;

      if domain1_right > domain2_left
        contact(idt + 1) = true;
      end

      % 5d, 5e, 5f
      q0 = u + dt * dudt + 0.5 * dt ^ 2 * d2udt2;
      q1 = dudt + dt * d2udt2;
      q2 = d2udt2;

      % Eq 6
      % assemble fully coupled simulations
      A_mat = zeros(num_nodes_total + num_contact_points);
      b = zeros(num_nodes_total + num_contact_points, 1);

      A_mat(1:num_nodes_total, 1:num_nodes_total) = M
      A_mat(1:num_nodes_total, end) = -dt * G';
      A_mat(end, 1:num_nodes_total) = -dt * G;
      A_mat(end, end) = 0

      m_inv_f = M \ (f_int + f_ext)

      vstar = v_m1 + dt * m_inv_f


      b(1:num_nodes_total) = M * vstar
      b(end) = 0;

      if contact(idt + 1)

      else
        v_new = v_start
      end


      for idomain = 1:num_domains
        idx_start = left_domain_indices(idomain);
        idx_end = right_domain_indices(idomain);
        XV{idomain}(:, idt + 1) = XV{idomain}(:, idt) + u_new(idx_start:idx_end);
        VV{idomain}(:, idt + 1) = dudt_new(idx_start:idx_end);
        AV{idomain}(:, idt + 1) = d2udt2_new(idx_start:idx_end);
      end

      % update the internal forces for each domain
      f_int(left_domain_indices(idomain):right_domain_indices(idomain)) = ...
        f_int(left_domain_indices(idomain):right_domain_indices(idomain)) + ...
        K_domains{idomain} * u_new(idx_start:idx_end);

    end
    output.tV = tV;
    output.contact = contact;
    output.contact_force = contact_force;
    output.XV = XV;
    output.VV = VV;
    output.AV = AV;

end

