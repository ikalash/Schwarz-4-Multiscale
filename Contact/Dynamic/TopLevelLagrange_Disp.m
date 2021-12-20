function [output] = TopLevelLagrange_Disp(settings)

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
  
    %a0 = settings.a0;
    %a1 = settings.a1;
  
    augmented_lagrange_tolerance = settings.augmented_lagrange_tolerance;
    max_lagrange_iterations = settings.max_lagrange_iterations;
    augmented_lagrange_mode = settings.augmented_lagrange_mode;
    relaxation_factor = settings.relaxation_factor;
    lambda = settings.lambda;
    
    % newmark
    beta = settings.beta;
    gamma = settings.gamma;

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
    K_aug = zeros(num_nodes_total + num_contact_points);
    M_aug = zeros(num_nodes_total + num_contact_points);
    C_aug = zeros(num_nodes_total + num_contact_points);
  
    contact_node_pair = [right_domain_indices(1), left_domain_indices(2)];
    contact = false(nt,1);
    contact_force = zeros(nt, num_contact_points);
    f_ext = zeros(num_nodes_total +1, 1);
    f_int = zeros(num_nodes_total +1, 1);
  
    delta_u = cell(num_domains,1);
    u = zeros(num_nodes_total + 1, 1);
    dudt = zeros(num_nodes_total + 1, 1);
    d2udt2 = zeros(num_nodes_total + 1, 1);
  
    %jonEasyEnd
  
    XV = cell(num_domains,1);
    VV = cell(num_domains,1);
    AV = cell(num_domains,1);
  
    K_domains = cell(num_domains, 1);
    C_domains = cell(num_domains, 1);
    M_domains = cell(num_domains, 1);
    A_domains = cell(num_domains, 1);
    
    for idomain = 1:num_domains
  
      xv_initial = linspace(x0_start(idomain), x1_start(idomain), num_nodes(idomain));
  
      XV{idomain} = zeros(nx(idomain) + 1, nt);
      VV{idomain} = zeros(nx(idomain) + 1, nt);
      XV{idomain}(:, 1) = xv_initial;
      VV{idomain}(:, 1) = v0(idomain);
      AV{idomain} = zeros(nx(idomain) + 1, nt);
  
      %K_ii = zeros(num_nodes(idomain));
      %C_ii = zeros(num_nodes(idomain));
      %M_ii = zeros(num_nodes(idomain));
    end
  
    prev_f_contact = 0.0;
    
    for idt = 1:(nt-1)
      disp(idt);
      % determine if contact occurs (2 conditions)
      % first check if the points are geometrically in contact
      
      domain_1_pt = XV{1}(end, idt);
      domain_2_pt = XV{2}(1, idt);
      
      domain_1_pt_future_est = domain_1_pt;
      domain_2_pt_future_est = domain_2_pt;
      
      %domain_1_pt_future_est = XV{1}(end, idt) + VV{1}(end, idt) * dt; 
      %domain_2_pt_future_est = XV{2}(1, idt) + VV{2}(1, idt) * dt;
      
      %domain_1_pt_future_est = XV{1}(end, idt) + VV{1}(end, idt) * dt + 0.5 * AV{1}(end, idt) * dt^2;  
      %domain_2_pt_future_est = XV{2}(1, idt) + VV{2}(1, idt) * dt + 0.5 * AV{2}(1, idt) * dt^2;

      if domain_1_pt > domain_2_pt
        contact(idt) = true;
      else
        contact(idt) = false;
      end
      % second check that there isnt a negative contact force from the previous time step
      % that is holding the two domains together
  
      %if idt ~= 1 && contact_force(idt - 1, 1) < 0
      %  contact(idt) = false;
      %end
      
  
      %if contact(idt)
      %  K_aug(num_nodes_total +1, contact_node_pair(1)) = 1;
      %  K_aug(num_nodes_total +1, contact_node_pair(2)) = -1;
      %  K_aug(contact_node_pair(1), num_nodes_total + 1) = -1;
      %  K_aug(contact_node_pair(2), num_nodes_total + 1) = 1;
      %  
      %  K_aug(num_nodes_total + 1, num_nodes_total + 1) =  0;
      %else
      %  K_aug(num_nodes_total +1, contact_node_pair(1)) = 0;
      %  K_aug(num_nodes_total +1, contact_node_pair(2)) = 0;
      %  K_aug(contact_node_pair(1), num_nodes_total + 1) = 0;
      %  K_aug(contact_node_pair(2), num_nodes_total + 1) = 0;
      %  K_aug(num_nodes_total + 1, num_nodes_total + 1) = 1;
      %end
      
      % assemble the element matrices [M, C, K]
      for idomain = 1:num_domains
        [K_ii, M_ii, C_ii] = FEMAssembly(XV{idomain}(:, idt), A, E, rho, a0, a1);
        K_domains{idomain} = K_ii;
        M_domains{idomain} = M_ii;
        C_domains{idomain} = C_ii;
        
        % backward euler
        A_ii = M_ii / dt^2 + C_ii / dt + K_ii;
        % midpoint
        %A_ii = K_ii + 2 / dt * C_ii + 4 / dt^2 * M_ii;
        % newmark beta
        %A_ii = M_ii / (beta * dt^2) + gamma * C_ii / (beta * dt) + K_ii;  

        A_domains{idomain} = A_ii;
        
        start_idx = left_domain_indices(idomain);
        end_idx = right_domain_indices(idomain);
        K_aug(start_idx:end_idx, start_idx:end_idx) = K_ii;
        M_aug(start_idx:end_idx, start_idx:end_idx) = M_ii;
        C_aug(start_idx:end_idx, start_idx:end_idx) = C_ii;
        u(start_idx:end_idx) = 0;
        u(end) = contact_force(idt);
        dudt(start_idx:end_idx) = VV{idomain}(:, idt);
        d2udt2(start_idx:end_idx) = AV{idomain}(:, idt);
      end
  
      if augmented_lagrange_mode
        %prev_f_contact = 0.0;
        f_ext(right_domain_indices(1)) = 0.0;
        f_ext(left_domain_indices(2)) = 0.0;
      end
      
      % midpoint rule
      %A_mat = K_aug + 2 / dt * C_aug + 4 / dt^2 * M_aug;  
      %b = f_ext - f_int + M_aug * (4 / dt^2 * u + 4 / dt * dudt + d2udt2) + C_aug * (2 / dt * u + dudt);
      % backward euler    
      A_mat = M_aug / dt^2 + C_aug / dt + K_aug;
      b = f_ext - f_int + M_aug * (u / dt^2 + dudt / dt) + C_aug * u / dt; 

      % use of newmark beta
      %A_mat = M_aug / (beta * dt^2) + gamma * C_aug / (beta * dt) + K_aug;
      %b = f_ext - f_int + ...
      %M_aug * (u / (beta * dt^2) + dudt / (beta * dt) +  (1 / (2*beta) - 1) * d2udt2) + ...
      %C_aug * (gamma * u / (beta * dt) - dudt* (1 - gamma / beta) - dt*d2udt2 * (1 - gamma / (2*beta)) );

      A_mat = sparse(A_mat);
      
      if augmented_lagrange_mode == false
        if contact(idt)
          A_mat(num_nodes_total +1, contact_node_pair(1)) = 1;
          A_mat(num_nodes_total +1, contact_node_pair(2)) = -1;
          A_mat(contact_node_pair(1), num_nodes_total + 1) = 1;
          A_mat(contact_node_pair(2), num_nodes_total + 1) = -1;

          A_mat(num_nodes_total + 1, num_nodes_total + 1) =  0;
          
          
          b(num_nodes_total + 1) = (domain_2_pt - domain_1_pt) * relaxation_factor;
        else
          A_mat(num_nodes_total +1, contact_node_pair(1)) = 0;
          A_mat(num_nodes_total +1, contact_node_pair(2)) = 0;
          A_mat(contact_node_pair(1), num_nodes_total + 1) = 0;
          A_mat(contact_node_pair(2), num_nodes_total + 1) = 0;
          A_mat(num_nodes_total + 1, num_nodes_total + 1) = 1;
          b(num_nodes_total + 1) = 0.0;
        end
      end
      
      if contact(idt)
        if augmented_lagrange_mode
          
          it_count = 0;
          while it_count < max_lagrange_iterations
            gap = domain_1_pt_future_est - domain_2_pt_future_est;
          
            f_contact = prev_f_contact + gap * lambda;
            
            % solve the system with the current estimate
            if f_contact > 0.0
              f_ext(right_domain_indices(1)) = -f_contact;
              f_ext(left_domain_indices(2)) = f_contact;
            else
              f_contact = 0.0;
              f_ext(right_domain_indices(1)) = 0.0;
              f_ext(left_domain_indices(2)) = 0.0;
            end
            % recalculate
            
            b = f_ext - f_int + M_aug * (u / dt^2 + dudt / dt) + C_aug * u / dt;

            
            u_new = zeros(num_nodes_total+1, 1);
            for idomain = 1:num_domains
              A_ii = A_domains{idomain};
              b_ii = b(left_domain_indices(idomain):right_domain_indices(idomain));
              u_new(left_domain_indices(idomain):right_domain_indices(idomain)) = A_ii \ b_ii;
            end
            u_new(num_nodes_total+1) = f_contact;
            domain_1_pt_future_est = XV{1}(end, idt) + u_new(right_domain_indices(1));
            domain_2_pt_future_est = XV{2}(1, idt) + u_new(left_domain_indices(2));
            
            if abs(domain_1_pt_future_est - domain_2_pt_future_est) * relaxation_factor < augmented_lagrange_tolerance
              break;
            end
            
            prev_f_contact = f_contact;
            it_count = it_count + 1;
          end
        else
          u_new = A_mat \ b;
        end
        
      else
                  
        A_mat_no_contact = sparse(A_mat(1:num_nodes_total, 1:num_nodes_total));
        b_no_contact = b(1:num_nodes_total);
        u_new = A_mat_no_contact \ b_no_contact;
      end
      
      if contact(idt)
        contact_force(idt, :) = u_new(num_nodes_total+1:end); 
      end
      
      for idomain = 1:num_domains

        delta_u{idomain} = u_new(left_domain_indices(idomain):right_domain_indices(idomain));
        
        XV{idomain}(:, idt + 1) = XV{idomain}(:, idt) + delta_u{idomain};
        % midpoint rule
        %VV{idomain}(:, idt + 1) = 2 / dt * delta_u{idomain} - VV{idomain}(:, idt);
        %AV{idomain}(:, idt + 1) = 4 / dt^2 * delta_u{idomain} - 4 / dt * VV{idomain}(:, idt) - AV{idomain}(:, idt);
        % backward euler 
        VV{idomain}(:, idt + 1) = 1 / dt * delta_u{idomain}; 
        AV{idomain}(:, idt + 1) = 1 / dt^2 * delta_u{idomain} - 1 / dt * VV{idomain}(:, idt);
        % newmark beta
        %AV{idomain}(:, idt + 1) = delta_u{idomain} / (beta * dt^2) - VV{idomain}(:, idt) / (beta * dt) - (1 / (2 * beta) - 1) * AV{idomain}(:, idt);
        %VV{idomain}(:, idt + 1) = gamma / (beta * dt) * delta_u{idomain} + VV{idomain}(:, idt) * (1 - gamma / beta) + dt * AV{idomain}(:, idt) * (1 - gamma / (2*beta));
        
        % update the internal forces for each domain
        f_int(left_domain_indices(idomain):right_domain_indices(idomain)) = ...
          f_int(left_domain_indices(idomain):right_domain_indices(idomain)) + ... 
          K_domains{idomain} * delta_u{idomain};
      end
    end
  
    output.tV = tV;
    output.XV = XV;
    output.VV = VV;
    output.AV = AV;
  
    output.contact_force = contact_force;
    output.contact = contact;
  
    output = compute_kinetic_energy(output, settings);
    
  end
