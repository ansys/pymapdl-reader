"""
.. _ref_pontoon:

Shell Static Analysis
~~~~~~~~~~~~~~~~~~~~~

Visualize a shell static analysis
"""
# download the pontoon example
from ansys.mapdl.reader import examples
pontoon = examples.download_pontoon()

###############################################################################
# Print the pontoon result
print(pontoon)

###############################################################################
# Plot the nodal displacement
pontoon.plot_nodal_solution(0, show_displacement=True, displacement_factor=100000)


###############################################################################
# print the available result types
pontoon.available_results


###############################################################################
# Plot the shell elements
pontoon.plot()

###############################################################################
# Plot the elastic strain and show exaggerated displacement
pontoon.plot_nodal_elastic_strain(0, 'eqv', show_displacement=True,
                                  displacement_factor=100000,
                                  overlay_wireframe=True,
                                  lighting=False,
                                  add_text=False,
                                  show_edges=True)
# Note: lighting is disabled here as it's too dark
