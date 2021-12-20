dt_values = [5.0e-7, 2.5e-7, 1.25e-7, 6.25e-8];
nx_values = [25, 50, 100, 200];
num_runs = length(nx_values);

%settings_filename = 'contact_settings_carpenter.yaml';
settings_filename = 'contact_settings.yaml';

settings = ReadYaml(settings_filename);

base_output_filename = 'output.mat';

base_output_dir = 'convergence_analysis/';

mkdir(base_output_dir)

for irun = 1:num_runs
  % modify the settings
  settings.dt = dt_values(irun);
  settings.nx{1} = nx_values(irun);
  settings.nx{2} = nx_values(irun);

  % run the code with the new parameters 
  schwarz_output = SchwarzMethod(settings);

  run_dir = "run_" + string(irun);
  
  current_output_dir = base_output_dir + run_dir;
  
  % make the directory
  mkdir(current_output_dir);
  
  % save the data
  save(current_output_dir + "/" + base_output_filename, 'settings', 'schwarz_output'); 


end
