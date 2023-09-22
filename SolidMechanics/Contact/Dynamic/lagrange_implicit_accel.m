function [output] = lagrange_implicit_accel(settings)

    % small double
    tol = 1.0e-15;
    
    % unpack the struct
    dt = settings.dt;
    t_end = settings.t_end;
    t_start = settings.t_start;
    nx = cell2mat(settings.nx);
    v0 = cell2mat(settings.v0);
    A = settings.A;
    E = settings.E;
    rho = settings.rho;
  
    max_iterations = settings.max_iterations;

    augmented_lagrange_mode = settings.augmented_lagrange_mode;
    relaxation_factor = settings.relaxation_factor;
    lambda = settings.lambda;
    
    % newmark
%     beta = settings.beta;
%     gamma = settings.gamma;
    beta0 = 0.5;
    beta1 = 0.5;

    zeta_target = settings.zeta_target;
    
    K_damping_contrib = settings.stiffness_damping_contribution;
    M_damping_contrib = 1 - K_damping_contrib;
    
    x0_start = cell2mat(settings.x0_start);
    x1_start = cell2mat(settings.x1_start);    
  
    num_domains = length(nx);
  
    %bar_stiffness = zeros(num_domains, 2);
    %bar_mass = zeros(num_domains, 2);
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
  
    %jonStart
    %num_contact_points = (num_domains * num_domains - num_domains) / 2;
    %jonEnd
  
    %jonEasyStart
    num_contact_points = 1;
    K = zeros(num_nodes_total + num_contact_points);
    M = zeros(num_nodes_total + num_contact_points);
    C = zeros(num_nodes_total + num_contact_points);
%   
    contact_node_pair = [right_domain_indices(1), left_domain_indices(2)];
    contact = false(nt,1);
    contact_force = zeros(nt, num_contact_points);
    f_ext = zeros(num_nodes_total +1, 1);
    f_int = zeros(num_nodes_total +1, 1);
    R = zeros(num_nodes_total +1, 1);
  
    delta_u = cell(num_domains,1);
    u = zeros(num_nodes_total + 1, 1);
    u_new = zeros(num_nodes_total + 1, 1);
    xx = zeros(num_nodes_total + 1, 1);
    dudt = zeros(num_nodes_total + 1, 1);
    d2udt2 = zeros(num_nodes_total + 1, 1);
  
    %jonEasyEnd
  
    XV = cell(num_domains,1);
    VV = cell(num_domains,1);
    AV = cell(num_domains,1);
  
    K_domains = cell(num_domains, 1);
    C_domains = cell(num_domains, 1);
    M_domains = cell(num_domains, 1);
    A_mat = zeros(num_nodes_total + num_contact_points);
    b = zeros(num_nodes_total + num_contact_points, 1);
%     A_domains = cell(num_domains, 1);

    b0 = 0.5 * dt ^ 2 * beta0;
    b1 = dt * beta1;
    b2 = 1;

    for idomain = 1:num_domains
  
      xv_initial = linspace(x0_start(idomain), x1_start(idomain), num_nodes(idomain));
  
      XV{idomain} = zeros(nx(idomain) + 1, nt);
      VV{idomain} = zeros(nx(idomain) + 1, nt);
      XV{idomain}(:, 1) = xv_initial;
      VV{idomain}(:, 1) = v0(idomain);
      AV{idomain} = zeros(nx(idomain) + 1, nt);
    end

    G = zeros(1, num_nodes_total);
    G(num_nodes(idomain)) = 1;
    G(num_nodes(idomain) + 1) = -1;

    for idomain = 1:num_domains
%         [K_ii, M_ii, C_ii] = FEMAssembly(XV{idomain}(:, idt), A, E, rho, a0, a1);
%         K_domains{idomain} = K_ii;
%         M_domains{idomain} = M_ii;
%         C_domains{idomain} = C_ii;
        
        start_idx = left_domain_indices(idomain);
        end_idx = right_domain_indices(idomain);
%         K(start_idx:end_idx, start_idx:end_idx) = K_ii;
%         M(start_idx:end_idx, start_idx:end_idx) = M_ii;
%         C(start_idx:end_idx, start_idx:end_idx) = C_ii;
        u(start_idx:end_idx) = 0;
        u_new(start_idx:end_idx) = 0;
%         u(end) = contact_force(idt);
        xx(start_idx:end_idx) = XV{idomain}(:, 1);
        dudt(start_idx:end_idx) = VV{idomain}(:, 1);
        d2udt2(start_idx:end_idx) = AV{idomain}(:, 1);
    end

    for idt = 1:(nt-1)
      disp(idt);

      domain_1_pt = XV{1}(end, idt);
      domain_2_pt = XV{2}(1, idt);
      
      if domain_1_pt > domain_2_pt
        contact(idt) = true;
      else
        contact(idt) = false;
      end
      
      % assemble the element matrices [M, C, K]
       for idomain = 1:num_domains
        [K_ii, M_ii, C_ii] = FEMAssembly(XV{idomain}(:, idt), A, E, rho, a0, a1);
        K_domains{idomain} = K_ii;
        M_domains{idomain} = M_ii;
        C_domains{idomain} = C_ii;
        
        start_idx = left_domain_indices(idomain);
        end_idx = right_domain_indices(idomain);
        K(start_idx:end_idx, start_idx:end_idx) = K_ii;
        M(start_idx:end_idx, start_idx:end_idx) = M_ii;
        C(start_idx:end_idx, start_idx:end_idx) = C_ii;
%         u(start_idx:end_idx) = 0;
%         u_new(start_idx:end_idx) = 0;
% %         u(end) = contact_force(idt);
%         xx(start_idx:end_idx) = XV{idomain}(:, 1);
%         dudt(start_idx:end_idx) = VV{idomain}(:, 1);
%         d2udt2(start_idx:end_idx) = AV{idomain}(:, 1);
      end
      
      q0 = u + dt * dudt + 0.5 * dt ^ 2 * d2udt2;
      q1 = dudt + dt * d2udt2;
      q2 = d2udt2;
 
      A_mat(1:num_nodes_total, 1:num_nodes_total) = b2 * M(1:num_nodes_total, 1:num_nodes_total) + b1 * C(1:num_nodes_total, 1:num_nodes_total) + b0 * K(1:num_nodes_total, 1:num_nodes_total);
      A_mat(1:num_nodes_total, end) = G';
%       A_mat(end, 1:num_nodes_total) = G;
%       A_mat(end, end) = G * q0(1:num_nodes_total, 1);
      A_mat(end, 1:num_nodes_total) = b0 * G;
%       A_mat(end, end) = 0;

      b = R - f_int - (M * q2 + C * q1 + K * q0);
      b(end) = -G * (q0(1:num_nodes_total,1) + xx(1:num_nodes_total,1));  
%       b(end) = 0;

      if contact(idt) 
        A_mat = sparse(A_mat);
        x_new = A_mat \ b;
        delta_d2udt2_new = x_new(1:num_nodes_total);
        contact_force(idt + 1, 1) = x_new(end);
      else % contact 
        A_mat_no_contact = sparse(A_mat(1:num_nodes_total, 1:num_nodes_total));
        x_new = A_mat_no_contact \ b(1:num_nodes_total);
        delta_d2udt2_new = x_new;
         
      end % contact

      u_new(1:num_nodes_total,1) = q0(1:num_nodes_total,1) + b0 * delta_d2udt2_new;
%       u(1:num_nodes_total,1) = zeros(num_nodes_total,1);
      dudt_new(1:num_nodes_total,1) = q1(1:num_nodes_total,1) + b1 * delta_d2udt2_new;
      d2udt2_new(1:num_nodes_total,1) = q2(1:num_nodes_total,1) + b2 * delta_d2udt2_new;

%        u(1:num_nodes_total,1) = q0(1:num_nodes_total,1) + b0 * delta_d2udt2_new;
      u(1:num_nodes_total,1) = zeros(num_nodes_total,1);
%       u(1:num_nodes_total,1) = u_new(1:num_nodes_total,1);
      dudt(1:num_nodes_total,1) = dudt_new(1:num_nodes_total,1);
      d2udt2(1:num_nodes_total,1) = d2udt2_new(1:num_nodes_total,1);
        
      for idomain = 1:num_domains
        idx_start = left_domain_indices(idomain);
        idx_end = right_domain_indices(idomain);
%         XV{idomain}(:, idt + 1) = u(idx_start:idx_end);
%         XV{idomain}(:, idt + 1) = XV{idomain}(:, idt) + u_new(idx_start:idx_end);
        XV{idomain}(:, idt + 1) =  u_new(idx_start:idx_end);
        VV{idomain}(:, idt + 1) = dudt_new(idx_start:idx_end);
        AV{idomain}(:, idt + 1) = d2udt2_new(idx_start:idx_end);
      end

        f_int(left_domain_indices(idomain):right_domain_indices(idomain)) = ...
        f_int(left_domain_indices(idomain):right_domain_indices(idomain)) + ...
        K_domains{idomain} * u_new(idx_start:idx_end);

    end

    output.tV = tV;
    output.XV = XV;
    output.VV = VV;
    output.AV = AV;
  
    output.contact_force = contact_force;
    output.contact = contact;
  
    output = compute_kinetic_energy(output, settings);
    
  end
