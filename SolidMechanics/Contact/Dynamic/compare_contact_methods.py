import numpy
import matplotlib
import matplotlib.pyplot as plt
import scipy.io as sio

# crude algorithm to determine impact and release times
# not super accurate for coarse time step problems
def compute_impact_and_release_times(tV, contact):
  contact_indices = numpy.argwhere(contact)

  first_contact_idx = contact_indices[0]
  last_contact_idx = contact_indices[-1]

  t_impact = tV[first_contact_idx][0]
  t_release = tV[last_contact_idx][0]

  return (t_impact, t_release)



# where we place the generated figures
output_dir = "output_low_speed/"

# name of the matlab output file to read with the sim data
output_filename = output_dir + "output.mat"


fig_size = (12, 9)
fig_dpi = 200
font_size = 12

matplotlib.rcParams.update({'font.size': font_size})

mat_contents = sio.loadmat(output_filename)

# name the legend entries
legend_entries = {}
legend_entries["analytical_output"] = "Analytical"
legend_entries["explicit_lagrange_output"] = "Explicit LM"
legend_entries["explicit_penalty_output"] = "Explicit Penalty"
legend_entries["implicit_penalty_output"] = "Implicit Penalty"
legend_entries["implicit_implicit_schwarz_output"] = "Implicit Implicit Schwarz"
legend_entries["explicit_implicit_schwarz_output"] =  "Explicit Implicit Schwarz"

# get the methods to plot
methods_to_plot = list(legend_entries.keys())


# add a dictionary for hard coded ylim for specific plots that need it
y_limits = {}

y_limits_high_speed = {}
y_limits_low_speed = {}

y_limits_high_speed["contact_force"] = [-50.0, 250.0]

#y_limits = y_limits_high_speed
y_limits = y_limits_low_speed



y_labels = {}

y_labels["contact"] = "contact"
y_labels["contact_force"] = "contact force, $\lambda_n$, (N)"
# XV_1 is the velocity of the rightmost node in domain 1
y_labels["XV_1"] = "bar 1 right position, $x_{1,r}$ (m)"
# VV_1 is the velocity of the rightmost node in domain 2
y_labels["VV_1"] = "bar 1 right velocity, $v_{1,r}$ (m)"
y_labels["kinetic_energy"] = "bar 1 kinetic energy, $KE$, (J)"
y_labels["potential_energy"] = "bar 1 potential energy,  $PE$, (J)"
y_labels["total_energy"] = "bar 1 total energy,  $E$, (J)"
y_labels["mass_averaged_velocity"] = "bar 1 mass averaged velocity, $\\bar{u}$, (m/s)"
y_labels["total_energy_error"] = "bar 1 total energy error"
y_labels["schwarz_iterations"] = "number of Schwarz iterations"

fields_to_plot = list(y_labels.keys())


figs = {}
axes = {}

for field in fields_to_plot:
  figs[field] = plt.figure(figsize=fig_size, dpi=fig_dpi)
  axes[field] = figs[field].add_subplot(111)

  #axes[field].ticklabel_format(axis="x", style="sci")


for idx in range(len(methods_to_plot)):

  method = methods_to_plot[idx]
  plot_label = legend_entries[method]

  tV = mat_contents[method]["tV"][0,0][0]
  tV_1 = tV

  contact = mat_contents[method]["contact"][0,0][:,0]
  contact_force = mat_contents[method]["contact_force"][0,0][:,0]

  #(t_impact, t_release) = compute_impact_and_release_times(tV, contact)

  #print(method + " |  t_release: " + str(t_release) + " | t_impact: " + str(t_impact))
  #print(tV.shape)
  #print(tV)
  #axes_contact.plot(tV, contact)
  #axes_contact_force.plot(tV, contact_force)

  analytical_total_energy = mat_contents["analytical_output"]["total_energy"][0,0][0]

  if method == "analytical_output":
    XV_1 = mat_contents[method]["XV_1"][0,0].flatten()
    VV_1 = mat_contents[method]["VV_1"][0,0].flatten()
    kinetic_energy = mat_contents[method]["kinetic_energy"][0,0].flatten()
    potential_energy = mat_contents[method]["potential_energy"][0,0].flatten()
    total_energy = mat_contents[method]["total_energy"][0,0].flatten()
    mass_averaged_velocity = mat_contents[method]["mass_averaged_velocity"][0,0].flatten()

    nt = len(total_energy)
    total_energy_error = numpy.zeros((nt, 1))
  elif "schwarz" in method:

    schwarz_iterations = mat_contents[method]["num_schwarz_iterations"][0,0].flatten()
    axes["schwarz_iterations"].plot(tV, schwarz_iterations, label=plot_label)

    XV_1 = mat_contents[method]["XV"][0,0][0,0][-1, :].flatten()
    VV_1 = mat_contents[method]["VV"][0,0][0,0][-1, :].flatten()

    kinetic_energy = mat_contents[method]["kinetic_energy"][0,0][0,0].flatten()
    potential_energy = mat_contents[method]["potential_energy"][0,0][0,0].flatten()
    total_energy = mat_contents[method]["total_energy"][0,0][0,0].flatten()
    mass_averaged_velocity = mat_contents[method]["mass_averaged_velocity"][0,0][0,0].flatten()

    total_energy_error = (total_energy - analytical_total_energy) / analytical_total_energy

  else:
        
    kinetic_energy = mat_contents[method]["kinetic_energy"][0,0][0,0].flatten()
    potential_energy = mat_contents[method]["potential_energy"][0,0][0,0].flatten()
    total_energy = mat_contents[method]["total_energy"][0,0][0,0].flatten()

    total_energy_error = (total_energy - analytical_total_energy) / analytical_total_energy

    mass_averaged_velocity = mat_contents[method]["mass_averaged_velocity"][0,0][0,0].flatten()
    XV_1 = mat_contents[method]["XV"][0,0][0,0][-1, :].flatten()
    VV_1 = mat_contents[method]["VV"][0,0][0,0][-1, :].flatten()

  print(method) 

  axes["contact"].plot(tV, contact)
  axes["contact_force"].plot(tV, contact_force, label=plot_label)
  axes["XV_1"].plot(tV, XV_1, label=plot_label)
  axes["VV_1"].plot(tV, VV_1, label=plot_label)
  axes["kinetic_energy"].plot(tV_1, kinetic_energy, label=plot_label)
  axes["potential_energy"].plot(tV, potential_energy, label=plot_label)
  axes["total_energy"].plot(tV, total_energy, label=plot_label)
  axes["total_energy_error"].plot(tV, total_energy_error, label=plot_label)
  axes["mass_averaged_velocity"].plot(tV, mass_averaged_velocity, label=plot_label)


for idx in range(len(fields_to_plot)):
  
  field = fields_to_plot[idx]
  axes[field].legend()
  axes[field].set_xlabel("time, t, s")
  axes[field].set_ylabel(y_labels[field])
  if field in y_limits:
    axes[field].set_ylim(y_limits[field])

  figs[field].savefig(output_dir + field)
  plt.close(figs[field])

