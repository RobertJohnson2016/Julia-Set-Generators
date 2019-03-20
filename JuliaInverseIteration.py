import matplotlib.pyplot as plt
import numpy as np
import cmath

"""Rounds a complex number to the nearest multiple of increment in both
the real and imaginary parts.
"""
def round_complex(num, increment):
    real_rounded = round(num.real/increment) * increment
    imag_rounded = round(num.imag/increment) * increment
    return(real_rounded + 1j*imag_rounded)

"""Generates an image whith points plotted at given coordinates.

complex_coordinates: list of complex numbers to plot
resolution: the width of each pixel in the complex plane.

Note: if the coordinates have already been rounded to some increment it
is recommended you use this increment for the resolution.
e.g. The output of the Modified Inverse Iteration Method is rounded to
some maximum level of detail.
"""
def image_plot(complex_coords, resolution = 0.01):
    complex_coords = [round_complex(z, resolution)
                      for z in complex_coords]

    real_coords = [z.real for z in complex_coords]
    imag_coords = [z.imag for z in complex_coords]

    real_max = max(real_coords)
    real_min = min(real_coords)
    imag_max = max(imag_coords)
    imag_min = min(imag_coords)

    image_width =  int(((real_max - real_min) / resolution) + 1)
    image_height = int(((imag_max-imag_min)/resolution) + 1)
    image_matrix = np.zeros([image_height, image_width])

    point_indexes = [[int((z.imag-imag_min)/resolution),
                      int((z.real-real_min)/resolution)]
                     for z in complex_coords]

    for index in point_indexes:
        image_matrix[index[0], index[1]] = 1

    plt.imshow(image_matrix,
               "bone_r",
               extent=[real_min, real_max, imag_min, imag_max])
    plt.savefig("image_name.pdf", dpi=500, bbox_inches="tight")
    # Save image
    plt.show() # View the image in a window

""" Generates and plots the approximation of the Julia set of z^2 + c
 with 'iteration_depth' iterations of the Inverse Iteration Method.
"""
def inverse_iteration(c, iteration_depth, resolution):
    julia_set = set()
    previous_preimage = {complex(1)}

    for i in range(iteration_depth):
        preimage = set()
        for z in previous_preimage:
            preimage.add(cmath.sqrt(z-c))
            preimage.add(-cmath.sqrt(z-c))
        julia_set = julia_set.union(preimage)
        previous_preimage = preimage

    coordinates = list(julia_set)
    # Convert to list to guarantee stability of iteration order
    image_plot(coordinates, resolution)

""" Generates and plots the approximation of the Julia set of z^2 + c 
with 'iteration_depth' iterations of the Modified Inverse Iteration
Method.
"""
def modified_inverse_iteration(c, iteration_depth, resolution):
    julia_set = set()
    previous_preimage = {complex(1)}

    for i in range(iteration_depth):
        preimage = set()
        for z in previous_preimage:
            preimage.add(round_complex(cmath.sqrt(z-c), resolution))
            preimage.add(round_complex(-cmath.sqrt(z-c), resolution))
        previous_preimage = preimage - julia_set
        julia_set = julia_set.union(preimage)

    coordinates = list(julia_set)
    # Convert to list to guarantee stability of iteration order
    image_plot(coordinates, resolution)

#Example of usage:
modified_inverse_iteration(0.3+0.55j, 500, 0.001)