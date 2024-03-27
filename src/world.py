import itertools
import math
from typing import Tuple
import numpy as np
import scipy.sparse
from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix


# TODO step and solver methods should not be in this class
class World2D:
    def __init__(self, size: Tuple[float, float], shape: Tuple[int, int], viscosity: float, diffusion: float,
                 dissipation_rate: float):
        # World parameters
        self.ndim = 2
        self.size = size  # L
        self.shape = shape  # N
        self.vshape = (self.ndim,) + shape  # Velocity fields are vector fields

        # Field
        self.velocity = np.zeros(self.vshape)  # U0

        # Substance
        self.density = np.zeros(shape)  # S0
        self.viscosity = viscosity  # visc
        self.diffusion = diffusion  # kS
        self.dissipation_rate = dissipation_rate  # aS

        self.force = np.zeros(self.vshape)  # TODO for future use
        self.force_source = (0, 0)  # TODO for future use

        self.voxel_size = [0] * len(size)  # D
        for i in range(len(size)):
            self.voxel_size[i] = size[i] / shape[i]

    def apply_force(self, position: Tuple[int, int], force: Tuple[float, float]):
        pass

    def step(self, dt: float):
        velocity_diffusion_poisson_array = self.create_2d_poisson_array(self.shape,
                dt * self.viscosity / (self.voxel_size[0] * self.voxel_size[0]),
                dt * self.viscosity / (self.voxel_size[1] * self.voxel_size[1]))

        field_diffusion_poisson_array = self.create_2d_poisson_array(self.shape,
                dt * self.diffusion / (self.voxel_size[0] * self.voxel_size[0]),
                dt * self.diffusion / (self.voxel_size[1] * self.voxel_size[1]))

        field_diffusion_poisson_array += np.identity(self.shape[0] * self.shape[1])  # TODO Check for sign errors
        velocity_diffusion_poisson_array += np.identity(self.shape[0] * self.shape[1])  # TODO Check for sign errors

        self.projection_poisson_array = self.create_2d_poisson_array(self.shape,
                                                                     -1 / (self.voxel_size[0] * self.voxel_size[0]),
                                                                     -1 / (self.voxel_size[1] * self.voxel_size[
                                                                         1]))  # TODO Check for sign errors

        velocity_diffusion_poisson_array = csc_matrix(velocity_diffusion_poisson_array)
        field_diffusion_poisson_array = csc_matrix(field_diffusion_poisson_array)
        self.projection_poisson_array = csc_matrix(self.projection_poisson_array)

        # TODO implement gravity and buoyancy

        # velocity steps (Vstep)
        for i in range(self.ndim):
            self.velocity[i] = self.add_force(self.velocity[i], self.force[i], dt)
        for i in range(self.ndim):
            self.velocity[i] = self.transport(self.velocity[i], self.velocity, dt, i+1)
        for i in range(self.ndim):
            self.velocity[i] = self.diffuse(self.velocity[i], velocity_diffusion_poisson_array)
        self.velocity = self.project_velocity(self.velocity, dt)

        # scalar field steps (Sstep)
        # self.density = self.add_force(self.density, self.force, dt)
        self.density = self.transport(self.density, self.velocity, dt, 0)
        self.density = self.diffuse(self.density, field_diffusion_poisson_array)
        self.density = self.dissipate_field(self.density, self.dissipation_rate, dt)

    def add_force(self, field, force, dt: float):
        return field + (force * dt)  # Elementwise addition and multiplication

    # The velocity and density stays constant after advection.
    # So trace back the particle in time and take its velocity and density
    # Note that this method is invoked for only one component of the velocity field at a time.
    # Appendix A
    def transport(self, field, velocity_total, dt: float, bound_type):
        # TODO Check for sign errors, seems to not work on the edges
        # FIXME Edge handling: for advecting velocities:
        #  The edges should "push back" against a velocity towards it. Note that the direction of the velocity needs
        #  to be a parameter for that. See also http://www.multires.caltech.edu/teaching/demos/java/FluidSolver.java
        back = np.zeros(field.shape)
        for i, j in itertools.product(range(field.shape[0]), range(field.shape[1])):
            x = (i, j)      # TODO Paper suggests adding +0.5, but this makes out results weird
            x_prev = self.trace_particle(x, velocity_total, dt)
            # x_prev = (x_prev[0], x_prev[1])
            back[i, j] = self.interpolate_field(x_prev, field, bound_type)
        return back

    # Appendix B, p123. right
    def project_velocity(self, velocity0, dt: float):
        div = self.divergence_velocity(velocity0, self.voxel_size)
        sol = self.solve_sparse_system(self.projection_poisson_array, div)
        sol = sol.reshape(self.shape)
        return velocity0 - np.gradient(sol, self.voxel_size[0], self.voxel_size[1])  # TODO implement bounds

    # Appendix B
    def diffuse(self, density0, poisson_array):
        return self.solve_sparse_system(poisson_array, density0)

    def dissipate_field(self, density0, dissipation_rate: float, dt: float):
        dissipation_matrix = np.full(self.shape, 1 + dt * dissipation_rate)
        density1 = density0 / dissipation_matrix  # Elementwise division
        return density1.reshape(density0.shape)  # Check whether correct reshape

    # Calculates the position of a particle in the field at the last timestep. Used for advection(transport)
    # p. 125 bottom right
    def trace_particle(self, pos: (int, int), velocity, dt, steps=1) -> (float, float):
        h = -dt / steps  # Note that we use -dt to traverse the fild in the opposite direction
        curr_pos = pos
        # TODO Cleaner code
        for i in range(steps):
            k1 = self.interpolate_velocity(curr_pos, velocity)
            k1 = (h * k1[0], h * k1[1])
            k2_pos = (curr_pos[0] + k1[0] / 2, curr_pos[1] + k1[1] / 2)
            k2 = self.interpolate_velocity(k2_pos, velocity)
            k2 = (h * k2[0], h * k2[1])
            curr_pos = (curr_pos[0] + k2[0], curr_pos[1] + k2[1])
        return curr_pos

    def interpolate_field(self, curr_pos, field, bound_type):
        # Variable names are labled as if top left is the least coordinate with x incrasing to the right and y increasing down
        # TODO Current boundaries are set to the last field inside
        # with this we can get coords outside the field and boundaries
        int_coords = (int(math.floor(curr_pos[0])), int(math.floor(curr_pos[1])))

        top_left = self.access_field_with_bound(field, int_coords, bound_type)
        top_right = self.access_field_with_bound(field, (int_coords[0] + 1, int_coords[1]), bound_type)
        bottom_left = self.access_field_with_bound(field, (int_coords[0], int_coords[1] + 1), bound_type)
        bottom_right = self.access_field_with_bound(field, (int_coords[0] + 1, int_coords[1] + 1), bound_type)

        ratio_down = curr_pos[0] - int_coords[0]
        ratio_right = curr_pos[1] - int_coords[1]

        back = (1 - ratio_down) * ((1 - ratio_right) * top_left + ratio_right * top_right) + \
               ratio_down * ((1 - ratio_right) * bottom_left + ratio_right * bottom_right)

        return back

    # bound_type 0 for density, 1 for velocity field horizontal, 2 for velocity field vertical
    def access_field_with_bound(self, field, int_coords, bound_type):
        coords_in_field = (np.clip(int_coords[0], 0, field.shape[0] - 1), np.clip(int_coords[1], 0, field.shape[1] - 1))
        # we are in a boundary corner
        if ((int_coords[0] < 0 and int_coords[1] < 0) or
                (int_coords[0] >= field.shape[0] and int_coords[1] < 0) or
                (int_coords[0] < 0 and int_coords[1] >= field.shape[1]) or
                (int_coords[0] >= field.shape[0] and int_coords[1] >= field.shape[1])):
            if bound_type != 0:
                back = 0
            else:
                back = field[coords_in_field]
        # we are on the left or right boundary
        elif int_coords[0] < 0 or int_coords[0] >= field.shape[0]:
            if bound_type == 1:
                back = -field[coords_in_field]
            else:
                back = field[coords_in_field]
        # we are on the upper or lower boundary
        elif int_coords[1] < 0 or int_coords[1] >= field.shape[1]:
            if bound_type == 2:
                back = -field[coords_in_field]
            else:
                back = field[coords_in_field]
        # we are inside the field
        else:
            back = field[int_coords]
        return back

    def interpolate_velocity(self, curr_pos, velocity):
        return self.interpolate_field(curr_pos, velocity[0], 1), self.interpolate_field(curr_pos, velocity[1], 2)

    def divergence_velocity(self, velocity, voxel_size):
        # TODO Test following line: divergence = np.sum(np.gradient(velocity, dx, dy, dz), axis=0)
        # TODO Implement bounds
        back = np.zeros((velocity.shape[1], velocity.shape[2]))
        for i, j in itertools.product(range(velocity.shape[1]), range(velocity.shape[2])):
            back[i, j] += self.access_field_with_bound(velocity[0], (i + 1, j), 1) / voxel_size[0]
            back[i, j] -= self.access_field_with_bound(velocity[0], (i - 1, j), 1) / voxel_size[0]

            back[i, j] += self.access_field_with_bound(velocity[1], (i, j + 1), 2) / voxel_size[1]
            back[i, j] -= self.access_field_with_bound(velocity[1], (i, j - 1), 2) / voxel_size[1]
        return back

    def solve_sparse_system(self, lhs, rhs):
        flat_rhs = rhs.flatten('C')
        sol = spsolve(lhs, flat_rhs)    # TODO SparseEfficiencyWarning: spsolve requires A be CSC or CSR matrix format
        return sol.reshape(rhs.shape)

    # Creates the LHS of the 2D Poisson system of linear equations
    # This stays constant over the runtime of the program and therefor needs to be invoked only once
    # Comparable to the top equation in Appendix B, BUT THE SIGNS OF A AND B ARE FLIPPED!!
    # TODO Implement other values for out of bounds values. Current implementation lets the bounds be walls (I think...)
    def create_2d_poisson_array(self, shape, x_factor, y_factor):
        back = np.zeros((shape[0] * shape[1], shape[0] * shape[1]))
        cr = 0  # Current row
        s1 = shape[1]
        for i in range(shape[0]):
            for j in range(shape[1]):
                # -u[i+1, j] + 2u[i, j] - u[i-1, j]
                if j - 1 >= 0:
                    back[cr, j - 1 + i * s1] += -x_factor
                else:
                    back[cr, j + i * s1] += -x_factor
                back[cr, j + i * s1] += 2 * x_factor
                if j + 1 < shape[1]:
                    back[cr, j + 1 + i * s1] += -x_factor
                else:
                    back[cr, j + i * s1] += -x_factor
                # -u[i, j+1] + 2u[i, j] - u[i, j-1]
                if i - 1 >= 0:
                    back[cr, j + (i - 1) * s1] += -y_factor
                else:
                    back[cr, j + i * s1] += -y_factor
                back[cr, j + i * s1] += 2 * y_factor
                if i + 1 < shape[0]:
                    back[cr, j + (i + 1) * s1] += -y_factor
                else:
                    back[cr, j + i * s1] += -y_factor
                cr += 1
        return back
