% using a lumped mass approximation
function [output] = compute_kinetic_energy(output, settings)

rho = settings.rho;
A = settings.A;
E = settings.E;
XV = output.XV;
VV = output.VV;

num_domains = length(XV);
mass_averaged_velocity = cell(num_domains,1);

kinetic_energy = cell(num_domains,1);
linear_momenta = cell(num_domains,1);
potential_energy = cell(num_domains, 1);
total_energy = cell(num_domains, 1);

lumped_kinetic_energy = cell(num_domains,1);
lumped_potential_energy = cell(num_domains, 1);
lumped_total_energy = cell(num_domains, 1);

orig_elem_lengths = cell(num_domains, 1);

for idomain = 1:num_domains
    num_elements = size(XV{idomain},1) - 1;
    nt = size(XV{idomain}, 2);
    kinetic_energy{idomain} = zeros(1, nt);
    potential_energy{idomain} = zeros(1, nt);
    
    for idt = 1:nt
        domain_kinetic_energy = 0.0;
        %domain_potential_energy = 0.0;
        for jelem = 1:num_elements
            Le = abs(XV{idomain}(jelem,idt) - XV{idomain}(jelem + 1, idt));
            [Ke, Me, Ce] = Element(Le, A, E, rho, 0.0, 0.0);
            ve = VV{idomain}(jelem:jelem + 1, idt);
            %xe = XV{idomain}(jelem:jelem + 1, idt) - XV{idomain}(jelem:jelem + 1, 1);
            %elem_potential_energy = 0.5 * xe' * Ke * xe;
            elem_kinetic_energy = 0.5 * ve' * Me * ve;
            domain_kinetic_energy = domain_kinetic_energy + elem_kinetic_energy;
            %domain_potential_energy = domain_potential_energy + elem_potential_energy;
        end
        kinetic_energy{idomain}(1, idt) = domain_kinetic_energy;
        %potential_energy{idomain}(1, idt) = domain_potential_energy;
    end
    %total_energy{idomain} = kinetic_energy{idomain} + potential_energy{idomain};
end
    
for idomain = 1:num_domains
  num_elements = size(XV{idomain},1) - 1;
  nodal_mass = zeros(size(XV{idomain}));
  nt = size(XV{idomain}, 2);
  element_potential_energy = zeros(num_elements, nt);

  orig_elem_lengths{idomain} = XV{idomain}(2:end,1) - XV{idomain}(1:end-1,1);

  for jelem = 1:num_elements
    L_elem = abs(XV{idomain}(jelem,:) - XV{idomain}(jelem + 1,:));
    elem_delta = L_elem - orig_elem_lengths{idomain}(jelem);
    elem_stiffness = settings.A * settings.E / orig_elem_lengths{idomain}(jelem);

    element_potential_energy(jelem, :) = 0.5 * elem_stiffness * elem_delta.^2;
    mass_elem = rho * A * L_elem;
    nodal_mass(jelem, :) = nodal_mass(jelem, :) + 0.5 * mass_elem;
    nodal_mass(jelem + 1, :) = nodal_mass(jelem + 1, :) + 0.5 * mass_elem;
  end
  nodal_kinetic_energy = 0.5 * nodal_mass.*VV{idomain}.^2;
  lumped_kinetic_energy{idomain} = sum(nodal_kinetic_energy, 1);
  linear_momenta{idomain} = sum(nodal_mass.*VV{idomain},1);
  mass_averaged_velocity{idomain} = sum(nodal_mass.*VV{idomain},1)./sum(nodal_mass,1);
  potential_energy{idomain} = sum(element_potential_energy,1);

  total_energy{idomain} = kinetic_energy{idomain} + potential_energy{idomain};
end

% do some more calcs
gap = XV{2}(1, :) - XV{1}(end, :);
gap_rate = VV{2}(1, :) - VV{1}(end, :);

output.mass_averaged_velocity = mass_averaged_velocity;
output.kinetic_energy = kinetic_energy;
output.potential_energy = potential_energy;
output.total_energy = total_energy;
output.linear_momenta = linear_momenta;

output.lumped_kinetic_energy = lumped_kinetic_energy;
%output.lumped_potential_energy = lumped_potential_energy;
%output.lumped_total_energy = lumped_total_energy;

output.gap = gap;
output.gap_rate = gap_rate;

end

