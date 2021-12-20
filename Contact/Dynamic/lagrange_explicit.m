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
    end

    u_new_star = dt^2 * M_inv * -f_int + u + dt * dudt;
    gap = -(u_new_star(num_nodes(1) + 1) - u_new_star(num_nodes(1)));
    
    if gap > 0.0
      contact(idt) = true;
    else
      contact(idt) = false;
    end
      
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

    contact_force(idt) = f_contact;

    for idomain = 1:num_domains
      idx_start = left_domain_indices(idomain);
      idx_end = right_domain_indices(idomain);
      
      XV{idomain}(:, idt+1) = u_new(idx_start:idx_end);
      VV{idomain}(:, idt+1) = (u_new(idx_start:idx_end) - u(idx_start:idx_end)) / dt;
      if idt > 1
        VV{idomain}(:, idt) = (XV{idomain}(:, idt + 1) - XV{idomain}(:, idt - 1)) / (2*dt);
        AV{idomain}(:, idt) = (XV{idomain}(:, idt+1) - 2*XV{idomain}(:, idt) +  XV{idomain}(:, idt-1)) / dt^2;
      end
      f_int(idx_start:idx_end) = f_int(idx_start:idx_end) + K_domains{idomain} * (u_new(idx_start:idx_end) - u(idx_start:idx_end));

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
