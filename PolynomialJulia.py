import matplotlib.pyplot as plt
import numpy as np
import cmath

"""Renders, with a given colour map, the Juila set for the polynomial
with the given coefficients.

The image will by default be plotted from -iteration_bound to
+iteration_bound on both axes. This will contain the entire Julia set.
"""
def render_julia_set(coeffs,
                     colour_map,
                     iterations,
                     real_bounds=[],
                     imag_bounds=[],
                     resolution=0.01):
    iteration_bound = max(1,
                          (4/abs(coeffs[0]))**(1/(len(coeffs)-2)),
                          (2/abs(coeffs[0])) *
                            sum([abs(z) for z in coeffs]))
    if not real_bounds:
        real_bounds = [-iteration_bound, iteration_bound]
    if not imag_bounds:
        imag_bounds = [-iteration_bound, iteration_bound]

    image_width = int(((real_bounds[1] - real_bounds[0])
                       / resolution) + 1)
    image_height = int(((imag_bounds[1] - imag_bounds[0])
                        / resolution) + 1)
    image_matrix = np.zeros([image_height, image_width])

    for m in range(image_height):
        for n in range(image_width):
            real_val = real_bounds[0] + n*resolution
            imag_val = imag_bounds[0] + m*resolution
            image_matrix[m,n] = colour_map(coeffs,
                                           real_val + 1j*imag_val,
                                           iterations,
                                           iteration_bound)
        print(str(round(100 * m / image_matrix.shape[0], 1)) + "%")
    plt.imshow(image_matrix,
               "bone_r",
               extent=[real_bounds[0],
                       real_bounds[1],
                       imag_bounds[0],
                       imag_bounds[1]])
    plt.savefig("image_name.pdf", dpi=1000, bbox_inches="tight")
    # Save image
    plt.show() # View the image in a window

"""A colour map for render_julia_set which colours the Julia set by the
Level Set Method.
"""
def escape_time(coeffs, z, max_iterations, iteration_bound):
    for i in range(max_iterations):
        z = sum([coeffs[i] * z**(len(coeffs)-1 - i )
                 for i in range(len(coeffs))])
        if abs(z) >= iteration_bound:
            return(i)
    return(i)

#Example of usage
coeffs = [1,0,0.25]
render_julia_set(coeffs, escape_time, 20, resolution=0.005)