function [f_ext, f_contact] = load_vector_assembly(x_new, XV, lambda, idt)
  num_domains = length(XV);
  nx = zeros(num_domains,1);
  for idomain = 1:num_domains
    nx(idomain) = size(XV{idomain}, 1) - 1;
  end
  num_nodes_total = sum(nx) + num_domains;
  left_domain_indices = ones(num_domains,1);
  right_domain_indices = ones(num_domains,1);

  for idomain = 1:num_domains
    if idomain == 1
      left_domain_indices(idomain) = 1;
      right_domain_indices(idomain) = nx(idomain) + 1;
    else
      left_domain_indices(idomain) = 1 + right_domain_indices(idomain - 1);
      right_domain_indices(idomain) = left_domain_indices(idomain) + nx(idomain);
    end
  end

  domain_spans = zeros(num_domains, 2);
  for idomain = 1:num_domains
    domain_spans(idomain, 1) = XV{idomain}(1, idt) + x_new(left_domain_indices(idomain));
    domain_spans(idomain, 2) = XV{idomain}(end, idt) +  x_new(right_domain_indices(idomain));
  end

  % compute contact forces
  [contact_forces, contact_sides] = compute_contact_forces_matrix(domain_spans, lambda);

  % reset external load vector
  f_ext = zeros(num_nodes_total,1);

  % recalculate contact forces
  for idomain = 1:num_domains
    for jdomain = 1:num_domains
      contact_side = contact_sides(idomain, jdomain);
      contact_force = contact_forces(idomain, jdomain);
      if contact_side == 0
        f_ext(left_domain_indices(idomain)) = f_ext(left_domain_indices(idomain)) + contact_force;
      elseif contact_side == 1
        f_ext(right_domain_indices(idomain)) = f_ext(right_domain_indices(idomain)) - contact_force;
      end
    end
  end

  f_contact = contact_forces(1,2);

end