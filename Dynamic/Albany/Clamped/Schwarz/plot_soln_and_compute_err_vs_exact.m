
close all;
clear all;

[x0,y0,z0,dispz0] = extract_fields_from_exo('clamped_left.e', 3);
[Nx, Nt] = size(dispz0);
ind0 = find((x0 == 0.0005) & (y0 == 0.0005));

[x1,y1,z1,dispz1] = extract_fields_from_exo('clamped_right.e', 3);
[Nx, Nt] = size(dispz1);
ind1 = find((x1 == 0.0005) & (y1 == 0.0005));

z = [0:0.001:1];
c = sqrt(1e9/1e3);
a = 0.01;
b = 0.5;
s = 0.02;
T = 1e-3;

fig1 = figure(1);
winsize = get(fig1,'Position');
Movie=moviein(Nt,fig1);
set(fig1,'NextPlot','replacechildren')
j = 1;
for i=1:Nt
  time = 10*(i-1)*(1e-6);
  times(i) = time;
  clearvars dispz;
  dispz = 1/2*a*(exp(-(z-c*time-b).^2/2/s^2) + exp(-(z+c*time-b).^2/2/s^2))...
      - 1/2*a*(exp(-(z-c*(T-time)-b).^2/2/s^2) + exp(-(z+c*(T-time)-b).^2/2/s^2));
  plot(z, dispz,'b');
  hold on;
  dz0 = dispz0(ind0, i);
  [zsort0, I] = sort(z0(ind0));
  dispzEx0 = 1/2*a*(exp(-(zsort0-c*time-b).^2/2/s^2) + exp(-(zsort0+c*time-b).^2/2/s^2))...
      - 1/2*a*(exp(-(zsort0-c*(T-time)-b).^2/2/s^2) + exp(-(zsort0+c*(T-time)-b).^2/2/s^2));
  norm_sol0(i) = norm(dispzEx0);
  abs_err0(i) = norm(dispzEx0-dz0(I));
  if (norm_sol0(i)  < 4.0e-3)
    rel_err0(i) = abs_err0(i);
  else
    rel_err0(i) = abs_err0(i)/norm_sol0(i);
  end
  plot(zsort0, dz0(I), '--r');
  hold on;
  dz1 = dispz1(ind1, i);
  [zsort1, I] = sort(z1(ind1));
  dispzEx1 = 1/2*a*(exp(-(zsort1-c*time-b).^2/2/s^2) + exp(-(zsort1+c*time-b).^2/2/s^2))...
      - 1/2*a*(exp(-(zsort1-c*(T-time)-b).^2/2/s^2) + exp(-(zsort1+c*(T-time)-b).^2/2/s^2));
  norm_sol1(i) = norm(dispzEx1);
  abs_err1(i) = norm(dispzEx1-dz1(I));
  if (norm_sol1(i)  < 4.0e-3)
    rel_err1(i) = abs_err1(i);
  else
    rel_err1(i) = abs_err1(i)/norm_sol1(i);
  end
  plot(zsort1, dz1(I), '-.g');
  set(fig1,'NextPlot','replacechildren')
  xlabel('z');
  ylabel('z-disp');
  title(['Time = ', num2str(time)]);
  axis([0 1.0 -0.01, 0.01]);
  legend('Exact', '\Omega_0', '\Omega_1', 'Location','Best');
  pause(0.1)
  Movie(:,j)=getframe(fig1);
  mov(j) = getframe(gcf);
  j = j+1;
end
figure();
plot(times, rel_err0) ;
hold on;
plot(times, rel_err1);
legend('\Omega_0','\Omega_1', 'Location','Best');
xlabel('time');
ylabel('Relative error (displacement)');
figure();
plot(times,norm_sol0);
hold on;
plot(times, norm_sol1);
legend('\Omega_0','\Omega_1', 'Location','Best');
xlabel('time');
ylabel('Norm displacement');
fprintf('Mean relative error in displacement in left domain = %f\n', mean(rel_err0));
fprintf('Max relative error in displacement in left domain = %f\n', max(rel_err0));
fprintf('Mean relative error in displacement in right domain  = %f\n', mean(rel_err1));
fprintf('Max relative error in displacement in right domain = %f\n', max(rel_err1));
fprintf('Mean relative error in displacement = %f\n', (mean(rel_err0)+mean(rel_err1))/2);
fprintf('Max relative error in displacement = %f\n', max(max(rel_err0), max(rel_err1)));

%movie2avi(Movie,'clamped_zdisp_1000_E1e9_zpts.avi','fps',7,'quality',10,'Compression','None');
