import numpy
import matplotlib
import matplotlib.pyplot as plt
import scipy.io as sio

fig_size = (12, 9)
fig_dpi = 200
font_size = 10

matplotlib.rcParams.update({'font.size': font_size})

method_to_analyze = "schwarz_output"

#number_of_cells = [20, 40, 80, 160, 320]
#time_step_sizes = [1.0e-6, 5.0e-7, 2.5e-7, 1.25e-7, 6.25e-8]

time_step_sizes = [5.0e-7, 2.5e-7, 1.25e-7, 6.25e-8]
number_of_cells = [25, 50, 100, 200]
num_runs = len(number_of_cells)

num_runs = len(number_of_cells)

convergence_fig = plt.figure(figsize=fig_size, dpi=fig_dpi)
convergence_axes = convergence_fig.add_subplot(111)


base_output_dir = 'convergence_analysis/'

for irun in range(num_runs):
  output_filename = base_output_dir + "run_" + str(irun + 1) + "/output.mat"

  mat_contents = sio.loadmat(output_filename)

  tV = mat_contents[method_to_analyze]["tV"][0,0][0].flatten()
  kinetic_energy = mat_contents[method_to_analyze]["total_energy"][0,0][0,0].flatten()

  print(tV.shape)
  print(kinetic_energy)

  plt_label = "$N_x = " + str(number_of_cells[irun]) + ", \Delta t = " + str(time_step_sizes[irun]) + "$"

  convergence_axes.plot(tV, kinetic_energy, label=plt_label)

convergence_axes.set_title("Convergence study of total energy over time for various time step and grid sizes")
convergence_axes.set_xlabel("time, t, sec")
convergence_axes.set_ylabel("total energy, E, J")

convergence_axes.legend()

plt.show()
   
