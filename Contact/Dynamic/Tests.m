% 
% % Test 1: 2 subdomains symmetric problem {Lr=Ll=0.25, v0_r=100, v0_r=-100, dt_r=dt_l=dt=10-7, EXP-EXP}
% close all; 
% clear all;
% delete U*txt
% delete V*txt
% delete A*txt
% delete KE*.txt
% delete PE*.txt
% fprintf('Starting test 1...\n')
% [nodal_fields, times] = Driver_run_tests('input/2sd_S_L-v100_L0.25_t10-7_EXP-R-v100_L0.25_t10-7_EXP.yaml');
% save('nodal_fields', 'times'); 
% fprintf('...finished test 1\n'); 
% 
% %% Test 2: 2 subdomains symmetric problem {Lr=Ll=0.25, v0_r=100, v0_r=-100, dt_r=dt_l=dt=10-7, IMP-IMP}
% close all; 
% clear all;
% delete U*txt
% delete V*txt
% delete A*txt
% delete KE*.txt
% delete PE*.txt
% fprintf('Starting test 2...\n')
% [nodal_fields, times] = Driver_run_tests('input/2sd_S_L-v100_L0.25_t10-7_IMP-R-v100_L0.25_t10-7_IMP.yaml');
% save('nodal_fields', 'times'); 
% fprintf('...finished test 2\n'); 
% 
% % Test 3: 2 subdomains non-symmetric problem {Lr=Ll=0.25, v0_r=0, v0_r=-200, dt_r=dt_l=dt=10-7, EXP-EXP}
% close all; 
% clear all;
% delete U*txt
% delete V*txt
% delete A*txt
% delete KE*.txt
% delete PE*.txt
% fprintf('Starting test 3...\n')
% [nodal_fields, times] = Driver_run_tests('input/2sd_NS_L-v0_L0.25_t10-7_EXP-R-v200_L0.25_t10-7_EXP.yaml');
% save('nodal_fields', 'times'); 
% fprintf('...finished test 3\n'); 
% 
% %% Test 4: 2 subdomains non-symmetric problem {Lr=Ll=0.25, v0_r=0, v0_r=-200, dt_r=dt_l=dt=10-7, IMP-IMP}
% close all; 
% clear all;
% delete U*txt
% delete V*txt
% delete A*txt
% delete KE*.txt
% delete PE*.txt
% fprintf('Starting test 4...\n')
% [nodal_fields, times] = Driver_run_tests('input/2sd_NS_L-v0_L0.25_t10-7_IMP-R-v200_L0.25_t10-7_IMP.yaml');
% save('nodal_fields', 'times'); 
% fprintf('...finished test 4\n'); 
% % 
% % Test 5: 2 subdomains non-symmetric problem {Lr=0.25, Ll=0.5, v0_r=100, v0_r=-100, dt_r=dt_l=dt=10-7, IMP-IMP}
% close all; 
% clear all;
% delete U*txt
% delete V*txt
% delete A*txt
% delete KE*.txt
% delete PE*.txt
% fprintf('Starting test 5...\n')
% [nodal_fields, times] = Driver_run_tests('input/2sd_NS_L-v100_L0.25_t10-7_IMP-R-v100_L0.5_t10-7_IMP.yaml');
% save('nodal_fields', 'times'); 
% fprintf('...finished test 5\n'); 

% Test 6: 2 subdomains non-symmetric problem {Lr=0.25, Ll=0.5, v0_r=100, v0_r=-100, dt_r=dt_l=dt=10-7, EXP-EXP}
close all; 
clear all;
delete U*txt
delete V*txt
delete A*txt
delete KE*.txt
delete PE*.txt
fprintf('Starting test 6...\n')
[nodal_fields, times] = Driver_run_tests('input/2sd_NS_L-v100_L0.25_t10-7_EXP-R-v100_L0.5_t10-7_EXP.yaml');
save('nodal_fields', 'times'); 
fprintf('...finished test 6\n'); 
% 
% % Test 7: 2 subdomains non-symmetric problem {Lr=0.25, Ll=0.5, v0_r=0, v0_r=-200, dt_r=dt_l=dt=10-7, EXP-EXP}
close all; 
clear all;
delete U*txt
delete V*txt
delete A*txt
delete KE*.txt
delete PE*.txt
fprintf('Starting test 7...\n')
[nodal_fields, times] = Driver_run_tests('input/2sd_NS_L-v0_L0.25_t10-7_EXP-R-v200_L0.5_t10-7_EXP.yaml');
save('nodal_fields', 'times'); 
fprintf('...finished test 7\n'); 

% % Test 8: 2 subdomains non-symmetric problem {Lr=0.25, Ll=0.5, v0_r=0, v0_r=-200, dt_r=dt_l=dt=10-7, IMP-IMP}
% close all; 
% clear all;
% delete U*txt
% delete V*txt
% delete A*txt
% delete KE*.txt
% delete PE*.txt
% fprintf('Starting test 8...\n')
% [nodal_fields, times] = Driver_run_tests('input/2sd_NS_L-v0_L0.25_t10-7_IMP-R-v200_L0.5_t10-7_IMP.yaml');
% save('nodal_fields', 'times'); 
% fprintf('...finished test 8\n'); 
% % %
% close all;
% 
