
import math
import matplotlib.pyplot as plt
import numpy as np
import os

from matplotlib import cm
from matplotlib.ticker import LinearLocator

from src.utility import get_start, get_terms

def get_values(xg, yg, terms, limits):
    """
    Calculate z=f(x,y) for x,y value, vector of values, or grid of values
    @param xg: x-axis coordinate (as value, vector, or grid)
    @param yg: y-axis coordinate (as value, vector, or grid)
    @param terms: list of list of terms=[amplitude, x_freq, y_freq, x_ctr, y_ctr]
    @param limits: x_min, x_max, y_min, y_max
    @return zg: value z=f(x,y) according to sum of implemented equation given terms

    """

    zg = np.array(0.0*xg) # Same size grid (force array type for single scalar)
    zg += terms[0][0] # Constant term A0
    for terms_list in terms[1:]:
       zg += add_term(xg, yg, terms_list )

    # Add barrier terms at the limits
    zg += add_limits(xg, yg, limits)

    return zg


def add_term(xg, yg, terms_list):
    """
    Calculate component for z=f(x,y) given compenent terms
    @param xg: x-axis coordinate (as value, vector, or grid)
    @param yg: y-axis coordinate (as value, vector, or grid)
    @param terms_list: list of terms=[amplitude, x_freq, y_freq, x_ctr, y_ctr]
    @return values z=f(x,y) given terms
    """
    # Extract components from list
    amp, wx, wy, xc, yc = terms_list

    # Calculate offset from centroid of this wave
    xv = xg - xc
    yv = yg - yc

    # Calculate wave from this term using numpy
    component_values = amp*np.cos(wx*xv)*np.sin(wy*yv)

    return component_values

def add_limits(xg, yg, limits):
    """
    A value due to boundary barrier term
    https://en.wikipedia.org/wiki/Barrier_function
    http://www.seas.ucla.edu/~vandenbe/ee236a/lectures/barrier.pdf
    @param xg: x-axis coordinate (as value, vector, or grid)
    @param yg: y-axis coordinate (as value, vector, or grid)
    @param terms: list of list of terms=[amplitude, x_freq, y_freq, x_ctr, y_ctr]
    @param limits: x_min, x_max, y_min, y_max
    @return zl : component due to boundaries
    """

    zl = np.array(0.0*xg)
    positions = [xmin_posn, xmax_posn, ymin_posn, ymax_posn]
    for posn in limits:
        half_vector = np.array((0 - posn[0], 0.0 - posn[1]))  # vector pointing toward origin
        half_vector /= np.linalg.norm(half_vector)
        xv = np.array(xg - posn[0])
        yv = np.array(yg - posn[1])

        # Add 0.001 so that barrier is 0.001 behind the specified position
        dist = np.array(0.001 + xv*half_vector[0] + yv*half_vector[1])

        # Use Numpy logical indices
        zd = np.array(dist > 0.0)

        zl[zd] += -0.01*np.log(dist[zd])
        zl[~zd] = np.inf
    return zl

def  enforce_bounds(current, limits, diff_vector):
    """
    Check current position is in bounds
    If out of bounds, pull back into bounds and clear the momentum vector
    @param current: current position
    @param limits: x_min, x_max, y_min, y_max
    @return tuple of position and momentum vector
    """

    zl = np.array(0.0*xg)
    positions = [xmin_posn, xmax_posn, ymin_posn, ymax_posn]
    for posn in limits:
        half_vector = np.array((0 - posn[0], 0.0 - posn[1]))  # vector pointing toward origin
        half_vector /= np.linalg.norm(half_vector)
        xv = current[0] - posn[0]
        yv = current[1] - posn[1]

        dist =0.001 + xv*half_vector[0] + yv*half_vector[1]

        if dist < 0:
            curr_posn = np.array(current)
            current += (-dist + 1e-9)*half_vector
            diff_vector = 0.0*diff_vector
            print(f"Enforcing limit at {posn} : {curr_posn} ==> {current}")

    return current, diff_vector


def calc_gradient(position, terms, limits):
    """
    Calculate gradient of z=f(x,y) for x,y value, vector of values, or grid of values
    @param position: numpy array of (x,y) position
    @param terms: list of list of terms=[amplitude, x_freq, y_freq, x_ctr, y_ctr]
    @param limits: x_min, x_max, y_min, y_max
    @return numpy array of 2 values (df/dx, df/dy) values as components of gradient

    """

    # Calculate the gradient due to boundary limits
    gradient = calc_gradient_limits(position, limits)

    # Note: Gradient of 0th term is (0,0)
    for terms_list in terms[1:]:
        grad_i = calc_gradient_term(position, terms_list )
        gradient += grad_i

    #print(" gradient=", gradient)
    return gradient

def calc_gradient_limits(position, limits):
    """
    Calculate gradient of z=f(x,y) for x,y value, vector of values, or grid of values
    @param position: numpy array of (x,y) position
    @param limits: list of limit positions =[xmin_posn, xmax_posn, ymin_posn, ymax_posn]
    @return tuple of (df/dx, df/dy) values as components of gradient
    """
    grad_limit = 0.0*position # wrong, but correct shape
    positions = [xmin_posn, xmax_posn, ymin_posn, ymax_posn]
    for posn in limits:
        half_vector = np.array((0 - posn[0], 0.0 - posn[1]))  # vector pointing toward origin
        half_vector /= np.linalg.norm(half_vector)
        xv = position[0] - posn[0]
        yv = position[1] - posn[1]

        dist = 0.001 + xv*half_vector[0] + yv*half_vector[1]

        # clamp how big we allow this gradient component to be
        if dist > 1e-3:
            grad_limit[0] += -0.01*half_vector[0]/dist
            grad_limit[1] += -0.01*half_vector[1]/dist
        else:
            grad_limit[0] += -0.01*half_vector[0]/1e-3
            grad_limit[1] += -0.01*half_vector[1]/1e-3

    #print(" grad_limit=", grad_limit)
    return grad_limit

def calc_gradient_term(position, terms_list):
    """
    Calculate gradient of z=f(x,y) for x,y value, vector of values, or grid of values
    @param position: numpy array of (x,y) position
    @param terms_list: list of terms=[amplitude, x_freq, y_freq, x_ctr, y_ctr]
    @return tuple of (df/dx, df/dy) values as components of gradient
    """

    # Extract terms and set up x,y relative to center point
    amp, wx, wy, xc, yc = terms_list
    xv = position[0] - xc
    yv = position[1] - yc

    grad_i = 0.0*position # wrong, but sets the correct array shape

    #@TODO - this is your part!
    #        Should be similar to add_term above, but will need 2 components for gradient

    component_values_DfDx = amp* (-np.sin(wx*xv))*np.sin(wy*yv)
    component_values_DfDy = amp*np.cos(wx*xv)*np.cos(wy*yv)

    grad_i[0] = component_values_DfDx  # fix this for df/dx
    grad_i[1] = component_values_DfDy  # fix this for df/dy

    # No need to change this return if you do calculation above correctly
    return grad_i


if __name__ == '__main__':


    # Parameters you should change during homework
    alpha =  0.05 # step size (learning rate) (1.e-8 < alpha < 0.99)
    mom   = 0.8 # Momentum term (used to blend current gradient with prior)
                 #  0 <= mom < 1.0

    # Change to choose different starting files
    # Or assign a different name to generate a random values
    start_file_name = 'start_point.txt'
    # start_file_name = 'start_point_two.txt'

    # Choose a file with terms for the equations
    # Note: the first (0th) term is a constant offset with frequency=0
    # terms_file_name = 'two_terms.txt'
    terms_file_name = 'three_terms.txt'
    #terms_file_name = 'four_terms.txt'

    # You can set to more than 4 if creating a new random file
    # If less terms in file, then only those if file are used
    # If more terms in file than max, then only first max_terms are used
    max_terms_in_equation = 4 # Keep at least 2 and no more than 6
    # (can be used to select subset of existing file)

    #-----------------------------------------------------
    # Terms you can change, but don't need to for this HW
    epsilon = 0.0001 # Convergence factor for norm (length) of gradient
                     # Controls how precise we want to be, but may take longer to converge
                     # also, for larger alphas, might just bounce around extrema

    max_step_count = 5000 # Might not converge if too low, but controls time
                          # and plotted points
                          # If not converging within this count for two term,
                          # then likely your gradient calculation is wrong

    minimize = True  # False # Maximize or minimize?


    num_points_per_axis = 100 # Change details on the graphs
                              # More points show details better, but
                              # make plots less responsive

    # x_min_value, x_max_value, y_min_value, y_max_value
    limits = (-10, 10, -10, 10)

    num_contours = 20

    ###############################################################
    ## No need to change below here, but you should read carefully
    ##

    # Re-assign as limit positions of half-space constraints
    x_min_value, x_max_value, y_min_value, y_max_value = limits
    xmin_posn = (x_min_value, 0.0)
    xmax_posn = (x_max_value, 0.0)
    ymin_posn = (0.0, y_min_value)
    ymax_posn = (0.0, y_max_value)
    limits = [xmin_posn, xmax_posn, ymin_posn, ymax_posn]

    # Get terms from file, or generate random for new file name
    terms = get_terms(terms_file_name, max_terms_in_equation)

    # Get start point from file, or generate random for new file name
    start = get_start(start_file_name)


    # Print terms of equations (without the barrier terms)
    print(' '.join([f'{val:8s}' for val in ['Amp', 'freq_x', 'freq_y', 'x_ctr', 'y_ctr']]))
    print(5*9*' ')
    for term in terms:
        print(' '.join([f'{val:8.4f}' for val in term]))
    print(5*9*' ')

    print(f"Starting point = ({start[0]:.3f}, {start[1]:.3f})")

    # Points along each axis
    xv = np.linspace(x_min_value, x_max_value, num_points_per_axis)
    yv = np.linspace(y_min_value, y_max_value, num_points_per_axis)

    # Make a 2D grid of points (num_points_per_axis squared values)
    # 2D array of values, so corresponding row, column from each xg, yg represents
    # a Cartesian pair (x, y)
    xg, yg = np.meshgrid(xv, yv, indexing='xy')

    # z=f(x, y) for all pairs in grid
    zg = get_values(xg, yg, terms, limits)


    # Initial values
    # Convert point values to individual numpy arrays to match the grid usage
    # so the code can be used for both grid and single points
    start_value = get_values(np.array((start[0],)), np.array((start[1],)), terms, limits)
    print(f"Starting value f({start[0]:.3f}, {start[1]:.3f})={start_value[0]:.6f}")

    # Initialize the optimization variables
    print("Setup for optimization ... ")
    current = np.array(start)  # Create copy of starting pose for updating
    value = np.array(start_value)
    gradient = calc_gradient(current, terms, limits)
    grad_norm = np.linalg.norm(gradient)
    diff_vector = 0.0*current  # zero initial vector for momentum calc


    # Initial surface plot without
    # Plot the surface.
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(xg, yg, zg,
                        alpha=0.5,
                        rstride=1, cstride=1,
                        cmap=cm.get_cmap('viridis'),
                        linewidth=0.01,
                        edgecolor='black',
                        antialiased=True,
                        rasterized=False)

    ax.set_zlim(1.01*np.min(np.min(zg)), 1.01*np.max(np.max(zg)))
    ax.zaxis.set_major_formatter('{x:.02f}')
    plt.xlabel('x')
    plt.ylabel('y')

    ax.plot(start[0], start[1], start_value[0], 'g.')
    # Update and show, but will not allow interaction yet
    plt.show(block=False)


    # Store path as list of x, list of y, list of z=f(x,y) values
    path = [ [start[0]], [start[1]], [start_value[0]]]

    # Main optimization loop
    step_count = 0
    print("Begin optimization loop ... ")
    while step_count < max_step_count and grad_norm > epsilon:
        step_count += 1


        if minimize:
            # Gradient descent to minimize (e.g. cost)
            #   vs. ascent to maximize (utility)
            # So just negate the gradient direction if minimizing
            gradient = -gradient

        # Calculate step with momentum
        # https://towardsdatascience.com/gradient-descent-with-momentum-59420f626c8f
        # Andrew Ng - https://www.youtube.com/watch?v=k8fTYJPd3_I
        diff_vector = mom*diff_vector + (1.0 - mom)*alpha*gradient

        #if step_count % 100 == 0:
        #    print(f" step_count = {step_count} |grad| = {grad_norm:.6f} ")
        print(f"{step_count:3d}: posn=({current[0]:.6f}, {current[1]:.6f}) dp=({diff_vector[0]:.6f}, {diff_vector[1]:.6f}) |grad| = {grad_norm:.6f} " )

        # Update the current position based on diff_vector
        current += diff_vector

        # Limit any steps that take us out of designated bounds
        # and reset the diff_vector in that case
        current, diff_vector = enforce_bounds(current, limits, diff_vector)

        # Update current value and gradient for the next step calculation
        # Using np.array for single x,y to allow logical indexing inside function
        value = get_values(np.array((current[0],)), np.array((current[1],)), terms, limits)
        gradient = calc_gradient(current, terms, limits)
        grad_norm = np.linalg.norm(gradient)

        # Store our path of ((x,y), z=f(x,y)) for later plotting
        path[0].append(current[0]) # x
        path[1].append(current[1]) # y
        path[2].append(value[0])   # z=f(x,y)

        # Plot progress as invidual points
        #  -- Forget for now to speed results --
        #if minimize:
        #    ax.scatter(current[0], current[1], value, 'r.')
        #else:
        #    ax.scatter(current[0], current[1], value, 'g.')
        #
        # Update but will not allow interaction yet
        #plt.show(block=False)

        # Check for runaway and break
        #   (should be protected by barrier and bounds check now)
        if np.fabs(current[0]) > 15.0:
            print(f"Going out of bounds, quit at {step_count}! {current}")
            break

        if np.fabs(current[1]) > 15.0:
            print(f"Going out of bounds, quit at {step_count}! {current}")
            break


    if minimize:
        print("Minimizing ...")
        print(f"    grid min = {np.min(zg):.6f}")  # Show true min of grid with given resolution
    else:
        print("Maximizing ...")
        print(f"    grid max = {np.max(zg):.6f}")  # Show true max of grid with given resolution

    print(f"alpha={alpha:.6f}, momentum={mom:.6f}, convergence={epsilon:.6f}")
    print(f" Steps={step_count}")
    print(f" Start f({start[0]:.4f}, {start[1]:.4f}) = {start_value[0]:.6f}")
    print(f" Final f({current[0]:.4f}, {current[1]:.4f}) = {value[0]:.6f} |grad|={grad_norm:.8f}")

    # Plot path with connected dashed line segments
    if minimize:
        ax.scatter(path[0], path[1], path[2], 'b.')
    else:
        ax.scatter(path[0], path[1], path[2], 'r.')

    # Show as contour plot
    fig, ax = plt.subplots()
    contours = ax.contour(xg, yg, zg, num_contours)
    ax.clabel(contours, inline=True, fontsize=10)
    if minimize:
        ax.plot(path[0], path[1], 'b:.')
    else:
        ax.plot(path[0], path[1], 'r:.')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_aspect('equal')

    base_terms_file = terms_file_name.split(".")[0]
    base_start_file = start_file_name.split(".")[0]

    title = f"{base_terms_file} alpha={alpha:.3f}, momentum={mom:.3f} steps={step_count} (justin.sanders.20)"
    ax.set_title(title)
    plot_file = f"{base_terms_file}_{base_start_file}_a{int(1000*alpha):04d}_m{int(1000*mom):04d}.png"
    # added relative path here, my machine was not saving
    fig.savefig(os.path.join("../plots", plot_file))

    # Interactive plots
    plt.show()
