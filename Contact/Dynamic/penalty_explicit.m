function [output] = penalty_explicit(settings)

settings.penalty_mode = true;
output = lagrange_explicit(settings);

end

