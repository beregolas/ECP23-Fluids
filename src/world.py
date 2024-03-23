import itertools
from typing import Tuple
import numpy as np
from scipy.sparse.linalg import spsolve


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
        self.diffusion_poisson_array = self.create_2d_poisson_array(self.shape, dt * self.diffusion / (
                    self.voxel_size[0] * self.voxel_size[0]), dt * self.diffusion / (
                                                                                self.voxel_size[1] * self.voxel_size[
                                                                            1]))
        self.diffusion_poisson_array += np.identity(self.shape[0] * self.shape[1])  # TODO Check for sign errors

        self.divergence_poisson_array = self.create_2d_poisson_array(self.shape,
                                                                     1 / (self.voxel_size[0] * self.voxel_size[0]),
                                                                     1 / (self.voxel_size[1] * self.voxel_size[
                                                                         1]))  # TODO Check for sign errors

        # velocity steps (Vstep)
        for i in range(self.ndim):
            self.velocity[i] = self.add_force_velocity(self.velocity[i], self.force[i], dt)
        for i in range(self.ndim):
            self.velocity[i] = self.transport_velocity(self.velocity[i], self.velocity, dt)
        for i in range(self.ndim):
            self.velocity[i] = self.diffuse_velocity(self.velocity[i], self.viscosity, dt)
        self.project_velocity(self.velocity, dt)

        # scalar field steps (Sstep)
        # self.density = self.add_force_field(self.density, self.force, dt)
        self.density = self.transport_field(self.density, self.velocity, dt)
        self.density = self.diffuse_field(self.density, self.diffusion, dt)
        self.density = self.dissipate_field(self.density, self.dissipation_rate, dt)

    def add_force_velocity(self, velocity0, force, dt: float):
        return velocity0 + (force * dt)  # Elementwise addition and multiplication

    # The velocity stays constant after advection. So trace back the particle in time and take its velocity
    # Note that this method is invoked for only one component of the velocity field at a time
    # Appendix A
    def transport_velocity(self, velocity0, velocity_total, dt: float):
        velocity1 = np.zeros(velocity0.shape)
        for i, j in itertools.product(range(velocity0.shape[0]), range(velocity0.shape[1])):
            x = (i, j)
            x_prev = self.trace_particle(x, velocity_total, dt)
            velocity1[i, j] = self.interpolate_field(x_prev, velocity0)
        return velocity1

    # Appendix B, p123. right
    def diffuse_velocity(self, velocity0, visc, dt: float):
        return self.solve_sparse_system(self.diffusion_poisson_array, velocity0)

    # Appendix B, p123. right
    def project_velocity(self, velocity0, dt: float):
        div = self.divergence_velocity(velocity0, self.voxel_size)
        sol = self.solve_sparse_system(self.divergence_poisson_array, div)
        sol = sol.reshape(self.shape)
        return velocity0 - np.gradient(sol, self.voxel_size[0], self.voxel_size[1])  # TODO validate

    def add_force_field(self, density0, source, dt: float):
        # source seems to be the external forces
        return density0 + source * dt

    # p 125 transport function
    def transport_field(self, density0, velocity, dt: float):
        density1 = np.zeros(density0.shape)
        for i, j in itertools.product(range(density0.shape[0]), range(density0.shape[1])):
            # TODO original code works on coordinated instead of indices. Check whether there is a difference
            x = (i, j)
            x_prev = self.trace_particle(x, velocity, dt)
            density1[i, j] = self.interpolate_field(x_prev, density0)
        return density1

    def diffuse_field(self, density0, diffusion: float, dt: float):
        return self.solve_sparse_system(self.diffusion_poisson_array, density0)

    def dissipate_field(self, density0, dissipation_rate: float, dt: float):
        dissipation_matrix = np.full(self.shape, 1 + dt * dissipation_rate)
        density1 = density0 / dissipation_matrix  # Elementwise division
        return density1.reshape(density0.shape)  # Check whether correct reshape

    # Calculates the position of a particle in the field at the last timestep. Used for advection(transport)
    # p. 125 bottom right
    def trace_particle(self, pos: (int, int), velocity, dt, steps=1) -> (float, float):
        h = -dt / steps  # Note that we use -dt to traverse the fild in the opposite direction
        curr_pos = pos
        # FIXME: Correct implementation of tuple addition/multiplication (pairwise)
        # FIXME: Find out how to handle float indices (maybe interpolate between neighboring cells?, Maybe round to next int?)
        # FIXME: velocity is a tuple of arrays
        for i in range(steps):
            k1 = self.interpolate_velocity(curr_pos, velocity)
            k1 = (h * k1[0], h * k1[1])
            k2_pos = (curr_pos[0] + k1[0] / 2, curr_pos[1] + k1[1] / 2)
            k2 = self.interpolate_velocity(k2_pos, velocity)
            k2 = (h * k2[0], h * k2[1])
            curr_pos = (curr_pos[0] + k2[0], curr_pos[1] + k2[1])
        return curr_pos

    def interpolate_field(self, curr_pos, field):
        # Variable names are labled as if top left is the least coordinate with x incrasing to the right and y increasing down
        # TODO Current boundaries are set to the last field inside
        top_left = (int(curr_pos[0]), int(curr_pos[1]))
        top_right = (min(top_left[0] + 1, field.shape[0] - 1), top_left[1])
        bottom_left = (top_left[0], min(top_left[1] + 1, field.shape[1] - 1))
        bottom_right = (min(top_left[0] + 1, field.shape[0] - 1), min(top_left[1] + 1, field.shape[1] - 1))

        ratio_down = curr_pos[0] - top_left[0]
        ratio_right = curr_pos[1] - top_left[1]

        back = (1 - ratio_down) * ((1 - ratio_right) * field[top_left] + ratio_right * field[top_right]) + \
               ratio_down * ((1 - ratio_right) * field[bottom_left] + ratio_right * field[bottom_right])

        return back

    def interpolate_velocity(self, curr_pos, velocity):
        return self.interpolate_field(curr_pos, velocity[0]), self.interpolate_field(curr_pos, velocity[1])

    def divergence_velocity(self, velocity, voxel_size):
        # TODO Test following line: divergence = np.sum(np.gradient(velocity, dx, dy, dz), axis=0)
        back = np.zeros((velocity.shape[1], velocity.shape[2]))
        for i, j in itertools.product(range(velocity.shape[1]), range(velocity.shape[2])):
            if i + 1 < velocity.shape[1]:
                back[i, j] += velocity[0, i + 1, j] / voxel_size[0]
            if i > 0:
                back[i, j] -= velocity[0, i - 1, j] / voxel_size[0]

            if j + 1 < velocity.shape[2]:
                back[i, j] += velocity[1, i, j + 1] / voxel_size[1]
            if j > 0:
                back[i, j] -= velocity[1, i, j - 1] / voxel_size[1]
        return back

    # Poisson 2d cartesian differential equation:
    # δ^2*u/δ^2x + δ^2*u/δ^2x = f(x, y)
    #
    # Finite difference equation TODO: validate
    # δ^2*u/δ^2x = (u[i+1, j] - 2u[i, j] + u[i-1, j]) / dx^2
    # δ^2*u/δ^2y = (u[i, j+1] - 2u[i, j] + u[i, j-1]) / dy^2
    def cartesian_poisson_2d(self, u, pos: (int, int)):
        pass

    def solve_sparse_system(self, lhs, rhs):
        flat_rhs = rhs.flatten('C')
        sol = spsolve(lhs, flat_rhs)
        return sol.reshape(rhs.shape)

    # Creates the LHS of the 2D Poisson system of linear equations
    # This stays constant over the runtime of the program and therefor needs to be invoked only once
    # Comparable to the top equation in Appendix B, BUT THE SIGNS OF A AND B ARE FLIPPED!!
    # TODO Implement other values for out of bounds values. Current implementation lets the bounds be walls (I think...)
    # TODO Find out why the original equations multiply with k
    def create_2d_poisson_array(self, shape, x_factor, y_factor):
        back = np.zeros((shape[0] * shape[1], shape[0] * shape[1]))
        cr = 0  # Current row
        s1 = shape[1]
        for i in range(shape[0]):
            for j in range(shape[1]):
                # -u[i+1, j] + 2u[i, j] - u[i-1, j]
                if j - 1 >= 0:
                    back[cr, j - 1 + i * s1] += -x_factor
                back[cr, j + i * s1] += 2 * x_factor
                if j + 1 < shape[1]:
                    back[cr, j + 1 + i * s1] += -x_factor

                # -u[i, j+1] + 2u[i, j] - u[i, j-1]
                if i - 1 >= 0:
                    back[cr, j + (i - 1) * s1] += -y_factor
                back[cr, j + i * s1] += 2 * y_factor
                if i + 1 < shape[0]:
                    back[cr, j + (i + 1) * s1] += -y_factor
                cr += 1
        return back
