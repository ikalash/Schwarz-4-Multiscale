function [output] = lagrange_implicit(settings)
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

%   augmented_lagrange_tolerance = settings.augmented_lagrange_tolerance;
%   max_lagrange_iterations = settings.max_lagrange_iterations;
  augmented_lagrange_mode = settings.augmented_lagrange_mode;
%   relaxation_factor = settings.relaxation_factor;
  lambda = settings.lambda;

  max_iterations = settings.max_iterations;
  newmark_tolerance = settings.newmark_tolerance;

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

  num_contact_points = 1;
  K_aug = zeros(num_nodes_total + num_contact_points);
  M_aug = zeros(num_nodes_total + num_contact_points);
  C_aug = zeros(num_nodes_total + num_contact_points);
%   G = zeros(num_nodes_total + num_contact_points);
  G = zeros(1, num_nodes_total);
  G(num_nodes(idomain)) = 1;
  G(num_nodes(idomain) + 1) = -1;


  K_domains = cell(num_domains, 1);
  C_domains = cell(num_domains, 1);
  M_domains = cell(num_domains, 1);

  contact_node_pair = [right_domain_indices(1), left_domain_indices(2)];
  contact = false(nt,1);
  contact_force = zeros(nt, num_contact_points);
  f_ext = zeros(num_nodes_total, 1);
  f_int = zeros(num_nodes_total, 1);
  
  delta_u = cell(num_domains,1);
  u = zeros(num_nodes_total + 1, 1);
  dudt = zeros(num_nodes_total + 1, 1);
  d2udt2 = zeros(num_nodes_total + 1, 1);

  r = zeros(num_nodes_total +1, 1);

  beta0 = 0.5;
  beta1 = 0.5;

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
    xv_initial = linspace(x0_start{idomain}, x1_start{idomain}, nx{idomain} + 1);
    XV{idomain} = zeros(nx{idomain} + 1, nt);
    VV{idomain} = zeros(nx{idomain} + 1, nt);
    AV{idomain} = zeros(nx{idomain} + 1, nt);
    % set initial conditions
    XV{idomain}(:, 1) = xv_initial;
    VV{idomain}(:, 1) = v0{idomain};
  end

  prev_f_contact = 0.0;
%   f_int = zeros(num_nodes_total, 1);

  for idt = 1:(nt-1)
%     if idt == 1  
        domain_1_pt = XV{1}(end, idt);
        domain_2_pt = XV{2}(1, idt);
%     else
%         domain_1_pt = XV{1}(end, idt) + XV{1}(end, idt-1);
%         domain_2_pt = XV{2}(1, idt) - XV{2}(1, idt-1);
%     end    
    domain_1_pt_future_est = domain_1_pt;
    domain_2_pt_future_est = domain_2_pt;
    if domain_1_pt > domain_2_pt
      contact(idt) = true; 
    else
      contact(idt) = false;
    end

    % assemble the block matrices
    for idomain = 1:num_domains
        xV = XV{idomain}(:, idt);
        idx_start = left_domain_indices(idomain);
        idx_end = right_domain_indices(idomain);
        [Ki, Mi, Ci] = FEMAssembly(xV, A, E, rho, a0, a1);
        %K_domains{idomain} = K_ii;
        %M_domains{idomain} = M_ii;
        %C_domains{idomain} = C_ii;
        
        K(idx_start:idx_end, idx_start:idx_end) = Ki;
        M(idx_start:idx_end, idx_start:idx_end) = Mi;
        C(idx_start:idx_end, idx_start:idx_end) = Ci;
    end
    
    K = sparse(K);
    M = sparse(M);
    C = sparse(C);

    x = zeros(num_nodes_total, 1);
    xx = zeros(num_nodes_total, 1);
    a = zeros(num_nodes_total, 1);
    v = zeros(num_nodes_total, 1);
%     lambda = zeros(num_nodes_total, 1);

    for idomain = 1:num_domains
      idx_start = left_domain_indices(idomain);
      idx_end = right_domain_indices(idomain);
      xx(idx_start:idx_end) = XV{idomain}(:, idt);
      a(idx_start:idx_end)  = AV{idomain}(:, idt);    
      v(idx_start:idx_end)  = VV{idomain}(:, idt);
%       lambda(idx_start:idx_end)  = LM{idomain}(:, idt);
        %u(start_idx:end_idx) = 0;
        %u(end) = contact_force(idt);
        %dudt(start_idx:end_idx) = VV{idomain}(:, idt);
       % d2udt2(start_idx:end_idx) = AV{idomain}(:, idt);
    end

%     b0 = 0.5 * dt ^ 2 * beta0;
%     b1 = dt * beta1;
%     b2 = 1;
%     G = zeros(1, num_nodes_total);
%     G(num_nodes(idomain)) = 1;
%     G(num_nodes(idomain) + 1) = -1;

    % backward euler    
    %A_mat = M_aug / dt^2 + C_aug / dt + K_aug;
    %b = f_ext - f_int + M_aug * (u / dt^2 + dudt / dt) + C_aug * u / dt; 

    % Newmark beta
    %A_mat = M_aug / beta / dt / dt + C_aug * gamma / beta / dt + K_aug;
    %b = f_ext - f_int + M_aug * (u / dt^2 + dudt / dt) + C_aug * u / dt;   
    
    %A_mat = sparse(A_mat);

    abs_error_tol = settings.newmark_tolerance;

    %%% iterative newmark implicit solve
    x_new = xx;
    v_new = v;
    x_prev = xx;
    v_prev = v;
    a_prev = a;
    lambda_new = contact_force;

    rel_tol_domain = 1e-12;
    abs_tol_domain = 1e-15;
    
    it_count = 1;
%     x_pre = x_prev + v_prev * dt + a_prev * (0.5 - beta) * dt * dt; 
%     v_pre = v_prev + a_prev * (1 - gamma) * dt;
%     a_new = (x_new - x_pre) / beta / dt / dt;

%     x_pre = x_prev + v_prev * dt + a_prev * 0.5 * dt * dt; 
%     v_pre = v_prev + a_prev * dt;
% %     a_new = (x_new - x_pre) / dt / dt;
%     a_new = a_prev;


%     q0 = x_prev + dt * v_prev + 0.5 * dt ^ 2 * a_prev;
%     q1 = v_prev + dt * a_prev;
%     q2 = a_prev;

%      while it_count < 2
%        [f_ext_new, f_contact] = load_vector_assembly((x_new - x_prev), XV, lambda, idt);
      
%        A_mat = M / beta / dt / dt + C * gamma / beta / dt + K;
%        A_mat = b2 * M(1:num_nodes_total, 1:num_nodes_total) + b1 * C(1:num_nodes_total, 1:num_nodes_total) + b0 * K(1:num_nodes_total, 1:num_nodes_total);
%        A_mat = M + dt * beta1 * C + dt * dt * 0.5 * beta0 * K;

       A_mat = M / dt^2 + C / dt + K;
       b = f_ext - f_int + M * (xx / dt^2 + v / dt) + C * xx / dt; 
%        fvis = C * v_new;
%        fine = M * a_new;
%        f_int_new = f_int + K * (x_new - x_prev) ;
%        b = f_ext_new - (f_int_new + fvis) - fine; 
%        b = f_ext_new - f_int_new - (fine + fvis);
      

       %A_mat = M / dt^2 + C / dt + K;
       %f_int_new = f_int + K * (x_new - x_prev) ;
       %b = f_ext_new - f_int_new + M * (x_new / dt^2 + x_new / dt) + C * x_new / dt; 

% 
       if contact(idt)
          A_mat(num_nodes_total +1, contact_node_pair(1)) = 1;
          A_mat(num_nodes_total +1, contact_node_pair(2)) = -1;
          A_mat(contact_node_pair(1), num_nodes_total + 1) = 1;
          A_mat(contact_node_pair(2), num_nodes_total + 1) = -1;
          A_mat(num_nodes_total + 1, num_nodes_total + 1) =  0;
          
          b(num_nodes_total + 1) = (domain_2_pt - domain_1_pt) ;
%           b(num_nodes_total + 1) = -G * (x_pre(1:num_nodes_total,1)); 
%         else
%           A_mat(num_nodes_total +1, contact_node_pair(1)) = 0;
%           A_mat(num_nodes_total +1, contact_node_pair(2)) = 0;
%           A_mat(contact_node_pair(1), num_nodes_total + 1) = 0;
%           A_mat(contact_node_pair(2), num_nodes_total + 1) = 0;
%           A_mat(num_nodes_total + 1, num_nodes_total + 1) = 0;
%           b(num_nodes_total + 1) = 0.0;
%           A_mat = sparse(A_mat(1:num_nodes_total, 1:num_nodes_total));
%           b = b(1:num_nodes_total);
       end

%        if contact(idt)
%          A_mat(num_nodes_total +1, contact_node_pair(1)) = 1 ;
%          A_mat(num_nodes_total +1, contact_node_pair(2)) = -1 ;
%          A_mat(contact_node_pair(1), num_nodes_total + 1) = 1 ;
%          A_mat(contact_node_pair(2), num_nodes_total + 1) = -1 ;
%          A_mat(num_nodes_total + 1, num_nodes_total + 1) =  0 ;
%          
%          b(num_nodes_total + 1, 1) = (domain_2_pt - domain_1_pt) ;
% %        else  
% %          A_mat = sparse(A_mat(1:num_nodes_total, 1:num_nodes_total));
% %          b = b(1:num_nodes_total);
%        end % contact
 
       %fvis = C * v_new;
       %fine = M * a_new;
       %f_int_new = f_int + K * (x_new - x_prev) ;
       %r = f_ext_new - (f_int_new + fvis) - fine;
       
       delta = A_mat \ b;
       xx = xx + delta(1:num_nodes_total, 1);
       a = 1 / dt^2 * delta(1:num_nodes_total, 1) - 1 / dt * v;
       v = 1 / dt * delta(1:num_nodes_total, 1);  
        

%        x_new = x_new + delta(1:num_nodes_total, 1);
%        a_new = (x_new - x_pre) / beta / dt / dt;
%        v_new = v_pre + gamma * dt * a_new; 
%        x_new = x_pre + b0 * delta(1:num_nodes_total, 1);
%        a_new = a_new + b1 * delta(1:num_nodes_total, 1);
%        v_new = v_pre + b2 * delta(1:num_nodes_total, 1);

%        if contact(idt)
%            x_new
% %          lambda_new = lambda_new + delta(num_nodes_total+1:end);
%        end

       %VV{idomain}(:, idt + 1) = 1 / dt * delta_u{idomain}; 
       %v_new = 1 / dt * delta; 
       %AV{idomain}(:, idt + 1) = 1 / dt^2 * delta_u{idomain} - 1 / dt * VV{idomain}(:, idt);
       %a_new = 1 / dt^2 * delta - 1 / dt * v_new;

%        if it_count ~= 1
%          abs_error = norm(delta(1:num_nodes_total, 1));
%          norm_posn = norm(x_new);
%          if norm_posn > 0.0
%            rel_error = abs_error / norm_posn;
%          else
%            if abs_error > 0
%              rel_error = 1.0;
%            else
%             rel_error = 0.0;
%            end
%          end
%          converged_rel = rel_error < rel_tol_domain;
%          converged_abs = abs_error < abs_tol_domain;
%          converged = (converged_rel || converged_abs);
%          if converged == 1
%            break;
%          end
%        end
%        it_count = it_count + 1;

%    for idomain = 1:num_domains
%     idx_start = left_domain_indices(idomain);
%     idx_end   = right_domain_indices(idomain);   
%     XV{idomain}(:, idt + 1) = xx(idx_start:idx_end);
% %     if idt == 1  
%         domain_1_pt = XV{1}(end, idt);
%         domain_2_pt = XV{2}(1, idt);
% %     else
% %         domain_1_pt = XV{1}(end, idt) + XV{1}(end, idt-1);
% %         domain_2_pt = XV{2}(1, idt) - XV{2}(1, idt-1);
% %     end
%    end

%    end

   if contact(idt)
     contact_force(idt, :) = delta(num_nodes_total+1:end);
   end

   for idomain = 1:num_domains
      idx_start = left_domain_indices(idomain);
      idx_end = right_domain_indices(idomain);
      AV{idomain}(:, idt + 1) = a(idx_start:idx_end);
      VV{idomain}(:, idt + 1) = v(idx_start:idx_end);
      XV{idomain}(:, idt + 1) = xx(idx_start:idx_end);

     %delta_u{idomain} = u_new(left_domain_indices(idomain):right_domain_indices(idomain));
     %XV{idomain}(:, idt + 1) = XV{idomain}(:, idt) + delta_u{idomain};

     % backward euler 
     %         VV{idomain}(:, idt + 1) = 1 / dt * delta_u{idomain}; 
     %         AV{idomain}(:, idt + 1) = 1 / dt^2 * delta_u{idomain} - 1 / dt * VV{idomain}(:, idt);
      % newmark beta
     %         VV{idomain}(:, idt + 1) = gamma / (beta * dt) * delta_u{idomain} + VV{idomain}(:, idt) * (1 - gamma / beta) + dt * AV{idomain}(:, idt) * (1 - gamma / (2*beta));
     %         AV{idomain}(:, idt + 1) = delta_u{idomain} / (beta * dt^2) - VV{idomain}(:, idt) / (beta * dt) - (1 / (2 * beta) - 1) * AV{idomain}(:, idt);
 
     %    AV{idomain}(:, idt + 1) = (XV{idomain}(:, idt + 1) - XV{idomain}(:, idt) - VV{idomain}(:, idt) * dt - AV{idomain}(:, idt) * (0.5 - beta) * dt * dt) / beta / dt / dt;
     %         AV{idomain}(:, idt + 1) = delta_u{idomain} / beta / dt / dt;
     %   VV{idomain}(:, idt + 1) = VV{idomain}(:, idt) + dt * (1 - gamma) * AV{idomain}(:, idt) + ... 
     %        gamma * dt * AV{idomain}(:, idt + 1);

     % update the internal forces for each domain
%      f_int = f_int + K * (x_new - xx);
        f_int = f_int + K * (delta(1:num_nodes_total, 1));
        
     %      f_int(left_domain_indices(idomain):right_domain_indices(idomain)) = ...
     %      f_int(left_domain_indices(idomain):right_domain_indices(idomain)) + ... 
     %      K_domains{idomain} * delta{idomain}; 
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