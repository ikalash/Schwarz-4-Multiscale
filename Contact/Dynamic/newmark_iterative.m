function [a_new, v_new, x_new, it_count, abs_error, f_contact] = newmark_iterative(M, C, K, f_int, a, v, x, XV, xx, abs_error_tol, max_iterations, beta, gamma, dt, lambda, idt, num_domains, left_domain_indices, right_domain_indices)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % Jonathan's code (+ uncomment changes in penalty_implicit.m)
%   a_new = a;
%   it_count = 1;
%   while it_count < max_iterations
%     v_new = v + (1 - gamma) * dt * a + gamma * dt * a_new;
%     x_new = x + dt * v + 0.5 * dt * dt * ((1 - 2 * beta) * a + 2 * beta * a_new);
%     f_int_new = f_int + K * (x_new - x);
%     [f_ext_new, f_contact] = load_vector_assembly(x_new, XV, lambda, idt);
%     rhs = f_ext_new - f_int_new - C * v_new;
%     a_new_old = a_new;
%     a_new = M \ rhs;
% 
%     % now compute the error from last iteration
%     if it_count ~= 1
%       abs_error = max(abs(a_new - a_new_old));
%       if abs_error < abs_error_tol
%         break;
%       end
%     end
%     it_count = it_count + 1;
%   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Daria's code - solve in terms of delta displacement (delta_x)

% Classic Newmark
    x_new = xx;
    v_new = v;
    x_prev = xx;
    v_prev = v;
    a_prev = a;

    rel_tol_domain = 1e-12;
    abs_tol_domain = 1e-15;
    
    it_count = 1;
    x_pre = x_prev + v_prev * dt + a_prev * (0.5 - beta) * dt * dt; 
    v_pre = v_prev + a_prev * (1 - gamma) * dt;
    a_new = (x_new - x_pre) / beta / dt / dt;

    while it_count < max_iterations

       [f_ext_new, f_contact] = load_vector_assembly((x_new - x_prev), XV, lambda, idt);

       TM = M / beta / dt / dt + C * gamma / beta / dt + K;
       fvis = C * v_new;
       fine = M * a_new;
       f_int_new = f_int + K * (x_new - x_prev) ;
       r = f_ext_new - (f_int_new + fvis) - fine;
       delta = TM \ r;
       x_new = x_new + delta;
% 
% % Different calculation of accelerations:
% % Classic Newmark      
       a_new = (x_new - x_pre) / beta / dt / dt;  
%%%        
% % Contact-implicit Newmark from Kane (modified) 
% % aims to set acceleration on the contact boundary to be equal to 0
% %      a_int = M \ (- f_int_new); 
% %      a_con = 2*((x_new - x_pre) / dt / dt - 2 * beta * a_int);
% %      a_new = a_int + a_con; 
%%%
       v_new = v_pre + gamma * dt * a_new; 
%        
% % now compute the error from last iteration
        if it_count ~= 1
          abs_error = norm(delta);
          norm_posn = norm(x_new);
          if norm_posn > 0.0
            rel_error = abs_error / norm_posn;
          else
            if abs_error > 0
              rel_error = 1.0;
            else
              rel_error = 0.0;
            end
          end
          converged_rel = rel_error < rel_tol_domain;
          converged_abs = abs_error < abs_tol_domain;
          converged = (converged_rel || converged_abs);
          if converged == 1
            break;
          end
        end
        it_count = it_count + 1;
    end

%     for idomain = 1:num_domains 
%           idx_start = left_domain_indices(idomain);
%           idx_end = right_domain_indices(idomain);
%           if idomain == 1
%               a_new(idx_end) = 0.;
%           else
%               a_new(idx_start) = 0.;
%           end  
%        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % Contact-implicit Newmark from Kane (original)
% % Classic Newmark
%     x_new = xx;
%     v_new = v;
%     x_prev = xx;
%     v_prev = v;
%     a_prev = a;
% 
%     rel_tol_domain = 1e-12;
%     abs_tol_domain = 1e-15;
%     
%     it_count = 1;
%     x_pre = x_prev + v_prev * dt ; 
%     v_pre = v_prev + a_prev * (1 - gamma) * dt;
%     a_prev_int = M \ (- f_int);
%     [f_ext_prev, f_contact] = load_vector_assembly(x_prev, XV, lambda, idt);
%     a_prev_con = - M \ f_ext_prev; 
%     a_new = a_prev_int + a_prev_con;
% 
%     while it_count < max_iterations
% 
%        [f_ext_new, f_contact] = load_vector_assembly((x_new - x_prev), XV, lambda, idt);
% 
%        TM = M / beta / dt / dt + C * gamma / beta / dt + K;
%        fvis = C * v_new;
%        fine = M * a_new;
%        x_tilde = x_pre + dt * dt * 0.5 * (1 - 2*beta) * a_prev_int;
%        f_int_new = f_int + K * (x_new - x_prev) ;
%        r = ((-2/dt/dt) * (M * x_new - M * x_tilde) + 2*beta*(f_ext_new - (f_int_new + fvis)) - fine);
%        delta = TM \ r;
%        x_new = x_new + delta;
% 
%        a_int = M \ (- f_int_new); 
%        a_con = 2*((x_new - x_pre) / dt / dt - 2 * beta * a_int);
%        a_new = a_int + a_con; 
% 
%        v_new = v_pre + dt * (1 - gamma) * a_prev_int + dt * gamma * a_int + dt * a_con; 
%        
% % now compute the error from last iteration
%         if it_count ~= 1
%           abs_error = norm(delta);
%           norm_posn = norm(x_new);
%           if norm_posn > 0.0
%             rel_error = abs_error / norm_posn;
%           else
%             if abs_error > 0
%               rel_error = 1.0;
%             else
%               rel_error = 0.0;
%             end
%           end
%           converged_rel = rel_error < rel_tol_domain;
%           converged_abs = abs_error < abs_tol_domain;
%           converged = (converged_rel || converged_abs);
%           if converged == 1
%             break;
%           end
%         end
%         it_count = it_count + 1;
%     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

end