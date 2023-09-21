function [output] = lagrange_explicit(settings)

  % FORWARD INCREMENT LAGRANGE MULTIPLIER METHOD

  % unpack the struct
  dt = settings.dt;
  t_end = settings.t_end;
  t_start = settings.t_start;
  nx = cell2mat(settings.nx);
  v0 = cell2mat(settings.v0);
  
  A = settings.A;
  E = settings.E;

  rho = settings.rho;
  lambda = settings.lambda;
  
  a0 = 0.0;
  a1 = 0.0;



  x0_start = cell2mat(settings.x0_start);
  x1_start = cell2mat(settings.x1_start);    

  num_domains = length(nx);

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

  %%%% CHECK STABLE TIME STEP
moduli    = E * ones(1, nx(1));
densities = rho * ones(1, nx(1));
%for domain = 1 : num_domains
 
    num_elements = nx(1);
    element_sizes = zeros(1, num_elements);
    X = linspace(x0_start(1), x1_start(1), num_nodes(1));

    for element = 1 : num_elements
      element_sizes(element) = X(element + 1) - X(element);
    end
    wave_speeds = sqrt(moduli ./ densities);
    stable_steps = element_sizes ./ wave_speeds;
    stable_step = 1.0 * min(stable_steps);
    tolerance = 1e-15;
%     if time_steps(domain) - stable_step > tolerance
%       fprintf(fmt_str, domain, time_steps(domain), stable_step);
%       error('Program terminated.');
%     end
  
% end
%%%%%%%

  nt = int32(round((t_end - t_start) / dt + 1)); 
  t_end_actual = round(double(nt - 1)) * dt + t_start;
  tV = linspace(t_start, t_end_actual, nt);

  contact = false(nt,1);
  contact_force = zeros(nt, 1);

  f_int = zeros(num_nodes_total, 1);

  XV = cell(num_domains,1);
  VV = cell(num_domains,1);
  AV = cell(num_domains,1);
  
  K_domains =  cell(num_domains, 1);
  M_domains =  cell(num_domains, 1);

  G_new = zeros(1, num_nodes_total);
  G_new(1, num_nodes(1)) = -1;
  G_new(1, num_nodes(2) + 1) = 1;

  for idomain = 1:num_domains

    xv_initial = linspace(x0_start(idomain), x1_start(idomain), num_nodes(idomain));

    XV{idomain} = zeros(nx(idomain) + 1, nt);
    VV{idomain} = zeros(nx(idomain) + 1, nt);
    XV{idomain}(:, 1) = xv_initial;
    VV{idomain}(:, 1) = v0(idomain);
    AV{idomain} = zeros(nx(idomain) + 1, nt);

  end

  dudt = zeros(num_nodes_total, 1);
  u = zeros(num_nodes_total, 1);
%   u_prev = zeros(num_nodes_total, 1);

  m_vec = zeros(num_nodes_total,1);
  
  for idt = 1:(nt-1)
    
     disp(idt)

    for idomain = 1:num_domains
      idx_start = left_domain_indices(idomain);
      idx_end = right_domain_indices(idomain);
      [K_ii, M_ii, C_ii] = FEMAssembly(XV{idomain}(:, idt), A, E, rho, a0, a1);
      m_vec(idx_start:idx_end) = sum(M_ii,2);
      K_domains{idomain} = K_ii;
      M_domains{idomain} = M_ii;
    end 

    M_inv = diag(m_vec.^-1);
    
    domain_1_pt = XV{1}(end, idt);
    domain_2_pt = XV{2}(1, idt);

    %if domain_1_pt > domain_2_pt
    %  contact(idt) = true;
    %else
    %  contact(idt) = false;
    %end

    %if contact(idt)
    %  G_new(1, num_nodes(1)) = -1;
    %  G_new(1, num_nodes(2) + 1) = 1;
    %else
    %  G_new(1, num_nodes(1)) = 0;
    %  G_new(1, num_nodes(2) + 1) = 0;
    %end  

    for idomain = 1:num_domains
      idx_start = left_domain_indices(idomain);
      idx_end = right_domain_indices(idomain);
    
      dudt(idx_start:idx_end) = VV{idomain}(:, idt);
      u(idx_start:idx_end) = XV{idomain}(:, idt);
%       if idt > 1
%        u_prev(idx_start:idx_end) = XV{idomain}(:, idt-1);
%       else
%        u_prev(idx_start:idx_end) = XV{idomain}(:, idt);
%       end    
    end

    u_new_star = dt^2 * M_inv * -f_int + u + dt * dudt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % Dissipative Tchamwa–Wielgosz scheme
%     phi = 1.1;
%     u_new_star1 = phi * (dt^2 * M_inv * -f_int) + u + dt * dudt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    gap = -(u_new_star(num_nodes(1) + 1) - u_new_star(num_nodes(1)));
    
%     if gap > -1.0e-6
%       contact(idt) = true;
%     else
%       contact(idt) = false;
%     end

%%% Introducing new contact conditions - the same as for Schwarz
%idt
      m1 = m_vec(right_domain_indices(1));
      m2 = m_vec(left_domain_indices(2));
      
      if settings.penalty_mode == false
        % contact force for lagrange explicit 
        f_contact = gap / (dt^2 * (1/m1 + 1/m2));
      else
        % results for penalty explicit
        f_contact = gap * lambda;
      end
      
      curr_posn_right = XV{2}(1, idt);
      curr_posn = XV{1}(end, idt);
      reac = -f_contact;
      overlap = curr_posn > curr_posn_right;
      compress = reac < 0.0; 
      persist = compress && contact(idt-1) == true;
      contact(idt) = overlap || persist;
%%%%%%%%

    u_new_c = zeros(num_nodes_total, 1);

    if contact(idt)
      %f_contact = dt^2 * inv(G_new * M_inv * G_new') * gap;
      %u_new_c = -dt^2 * M_inv * G_new' * f_contact; 
      m1 = m_vec(right_domain_indices(1));
      m2 = m_vec(left_domain_indices(2));
      

      if settings.penalty_mode == false
        % contact force for lagrange explicit 
        f_contact = gap / (dt^2 * (1/m1 + 1/m2));
      else
        % results for penalty explicit
        f_contact = gap * lambda;
      end
      
      if f_contact < 0.0
        contact(idt) = false;
        f_contact = 0.0;
      end

      u_new_c(right_domain_indices(1)) = -dt^2 * f_contact * (1 / m1);
      u_new_c(left_domain_indices(2)) = dt^2 * f_contact * (1 / m2);
    
    else
      f_contact = 0.0;
    end

    
    u_new = u_new_star + u_new_c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % Dissipative Tchamwa–Wielgosz scheme
%     u_new = u_new_star1 + phi*u_new_c;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     a_int_new = M_inv * -f_int ;
%     a_con_new = 2*(u_new_c / dt / dt - 0.5 * a_int_new);
%     a = a_int_new + a_con_new;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%


    contact_force(idt) = f_contact;

    for idomain = 1:num_domains
      idx_start = left_domain_indices(idomain);
      idx_end = right_domain_indices(idomain);
      
      XV{idomain}(:, idt+1) = u_new(idx_start:idx_end);
      VV{idomain}(:, idt+1) = (u_new(idx_start:idx_end) - u(idx_start:idx_end)) / dt;
      if idt > 1
        VV{idomain}(:, idt) = (XV{idomain}(:, idt + 1) - XV{idomain}(:, idt - 1)) / (2*dt);
        AV{idomain}(:, idt) = (XV{idomain}(:, idt+1) - 2*XV{idomain}(:, idt) +  XV{idomain}(:, idt-1)) / dt^2;
%         AV{idomain}(:, idt+1) = a(idx_start:idx_end);
      end
      f_int(idx_start:idx_end) = f_int(idx_start:idx_end) + K_domains{idomain} * (u_new(idx_start:idx_end) - u(idx_start:idx_end));

%       if idomain == 1
%         AV{idomain}(idx_end, idt) = 0.;
%       else
%         AV{idomain}(1, idt) = 0.;  
%       end    
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
