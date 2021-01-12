"""
.. _ref_load_shaft_result:

Shaft Modal Analysis
~~~~~~~~~~~~~~~~~~~~

Visualize a shaft modal analysis

"""
# sphinx_gallery_thumbnail_number = 6

from ansys.mapdl.reader import examples

# Download an example shaft modal analysis result file
shaft = examples.download_shaft_modal()

###############################################################################
# Mesh is stored within the result object
print(shaft.mesh)

###############################################################################
# ...and contains a VTK unstructured grid
print(shaft.mesh._grid)

###############################################################################
# Plot the shaft
cpos = shaft.plot()

# list shaft node components
###############################################################################
print(shaft.element_components.keys())

###############################################################################
# Plot a node component
#
# This camera angle was saved interactively from ``shaft.plot``
cpos = [(-115.35773008378118, 285.36602704380107, -393.9029392590675),
        (126.12852038381345, 0.2179228023931401, 5.236408799851887),
        (0.37246222812978824, 0.8468424028124546, 0.37964435122285495)]
shaft.plot(element_components=['SHAFT_MESH'], cpos=cpos)
# get cpos from cpos = shaft.plot()


###############################################################################
# Plot a node component as a wireframe
shaft.plot(element_components=['SHAFT_MESH'], cpos=cpos, style='wireframe',
           lighting=False)


###############################################################################
# Plot the shaft with edges and with a blue color
shaft.plot(show_edges=True, color='cyan')

###############################################################################
# Plot the shaft without lighting but with edges and with a blue color
shaft.plot(lighting=False, show_edges=True, color='cyan')


###############################################################################
# plot a mode shape without contours using the "bwr" color map
shaft.plot_nodal_solution(9, element_components=['SHAFT_MESH'],
                          show_displacement=True, cmap='bwr',
                          displacement_factor=0.3, stitle=None,
                          overlay_wireframe=True, cpos=cpos)

###############################################################################
# plot a mode shape with contours and the default colormap
shaft.plot_nodal_solution(1, element_components=['SHAFT_MESH'],
                          n_colors=10, show_displacement=True,
                          displacement_factor=1, stitle=None,
                          overlay_wireframe=True, cpos=cpos)

###############################################################################
# Animate a mode of a component the shaft
#
# Set ``loop==True`` to plot continuously.
# Disable ``movie_filename`` and increase ``nangles`` for a smoother plot
shaft.animate_nodal_solution(5, element_components='SHAFT_MESH',
                             comp='norm', displacement_factor=1,
                             show_edges=True, cpos=cpos,
                             loop=False, movie_filename='demo.gif',
                             nangles=30)
