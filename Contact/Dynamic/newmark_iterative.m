function [a_new, v_new, x_new, it_count, abs_error, f_contact] = newmark_iterative(M, C, K, f_int, a, v, x, XV, abs_error_tol, max_iterations, beta, gamma, dt, lambda, idt)
  a_new = a;
  it_count = 1;
  while it_count < max_iterations
    v_new = v + (1 - gamma) * dt * a + gamma * dt * a_new;
    x_new = x + dt * v + 0.5 * dt * dt * ((1 - 2 * beta) * a + 2 * beta * a_new);
    f_int_new = f_int + K * (x_new - x);
    [f_ext_new, f_contact] = load_vector_assembly(x_new, XV, lambda, idt);
    rhs = f_ext_new - f_int_new - C * v_new;
    a_new_old = a_new;
    a_new = M \ rhs;

    % now compute the error from last iteration
    if it_count ~= 1
      abs_error = max(abs(a_new - a_new_old));
      if abs_error < abs_error_tol
        break;
      end
    end
    it_count = it_count + 1;
  end
  %disp(it_count)
  %disp(abs_error)

end