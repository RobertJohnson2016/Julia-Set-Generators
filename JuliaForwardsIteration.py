import matplotlib.pyplot as plt
import numpy as np
import cmath

"""Renders the Juila set for value c with a given colour map.

The image will by default be plotted from -iteration_bound to 
+iteration_bound on both axes. This will contain the entire Julia set.
"""
def render_julia_set(c,
                     colour_map,
                     iterations,
                     real_bounds=[],
                     imag_bounds=[],
                     resolution=0.01):
    iteration_bound = max(2, np.sqrt(2 * abs(c)))
    if not real_bounds:
        real_bounds = [-iteration_bound, iteration_bound]
    if not imag_bounds:
        imag_bounds = [-iteration_bound, iteration_bound]

    image_width = int(((real_bounds[1] - real_bounds[0]) / resolution)
                      + 1)
    image_height = int(((imag_bounds[1] - imag_bounds[0]) / resolution)
                       + 1)
    image_matrix = np.zeros([image_height, image_width])

    for m in range(image_height):
        for n in range(image_width):
            real_val = real_bounds[0] + n*resolution
            imag_val = imag_bounds[0] + m*resolution
            image_matrix[m,n] = colour_map(c,
                                           real_val + 1j*imag_val,
                                           iterations)
        print(str(round(100 * m / image_matrix.shape[0], 1)) + "%")
        #Print the progress percentage.
    plt.imshow(image_matrix,
               "bone_r",
               extent=[real_bounds[0],
                       real_bounds[1],
                       imag_bounds[0],
                       imag_bounds[1]])
    plt.savefig("image_name.pdf", dpi=1000, bbox_inches="tight")
    # Save image
    plt.show() # View the image in a window

"""A colour map for render_julia_set which colours the Julia set of
z^2 + c by the Level Set Method.
"""

def escape_time(c, z, max_iterations):
    iteration_bound = max(4, np.sqrt(2*abs(c)))
    for i in range(max_iterations):
        z = z * z + c
        if abs(z) >= iteration_bound:
            return(i)
    return(i)

"""A colour map for render_julia_set which colours the Julia set of 
z^2 + c by the Distance Estimation Method.
 
Points in the Julia set are given a colour value of -2, not 0, to give
contrast on the boundary of the set.
"""
def distance_estimation(c, z, max_iterations):
    iteration_bound = max(4, np.sqrt(2*abs(c)))
    z_old = z
    dz_old = 1

    for i in range(max_iterations):
        z_new = z_old*z_old + c
        dz_new = 2*z_old*dz_old
        z_old = z_new
        dz_old = dz_new
        if abs(z_new) >= iteration_bound:
            return (abs(z_new) * np.log(abs(z_new)) / abs(dz_new))
    return (-2) #If bound isn't exceeded then z_0 is likely in the set

#Example of usage:
render_julia_set(0.3, escape_time, 50, resolution=0.005)