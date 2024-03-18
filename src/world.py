import itertools
from typing import Tuple
import numpy as np


class World2D:
    def __init__(self, size: Tuple[float, float], shape: Tuple[int, int], viscosity: float, diffusion: float,
                 dissipation_rate: float):
        # World parameters
        self.size = size                    # L
        self.shape = shape                  # N

        # Field
        self.velocity0 = np.zeros(shape)    # U0
        self.velocity1 = np.zeros(shape)    # U1

        # Substance
        self.density0 = np.zeros(shape)     # S0
        self.density1 = np.zeros(shape)     # S1
        self.viscosity = viscosity          # visc
        self.diffusion = diffusion          # kS
        self.dissipation_rate = dissipation_rate    # aS

        self.force = 0  # TODO for future use
        self.force_source = (0, 0)  #TODO for future use

        self.voxel_size = [0] * len(size)   # D
        for i in range(len(size)):
            self.voxel_size[i] = size[i] / shape[i]

    def apply_force(self, position: Tuple[int, int], force: Tuple[float, float]):
        pass

    def step(self, dt: float):
        # switch fields, older field is in [...]1
        self.velocity0, self.velocity1 = self.velocity1, self.velocity0
        self.density0, self.density1 = self.density1, self.density0
        # velocity steps (Vstep)
        for i in range(2):
            self.add_force_velocity(self.velocity0[i], self.force, dt)
        for i in range(2):
            self.transport_velocity(self.velocity0[i], self.velocity1[i], self.velocity0, dt)
        for i in range(2):
            self.diffuse_velocity(self.velocity0[i], self.velocity1[i], self.viscosity, dt)
        self.project_velocity(self.velocity0, self.velocity1, dt)

        # scalar field steps (Sstep)
        self.density0 = self.add_force_field(self.density0, self.force_source, dt)
        self.transport_field(self.density0, self.density1, self.velocity1, self.dt)
        self.diffuse_field(self.density0, self.density1, self.diffusion, self.dt)
        self.dissipate_field(self.density0, self.density1, self.dissipation_rate, self.dt)

    def add_force_velocity(self, velocity0, force, dt: float) -> None:
        pass

    #Appendix A
    def transport_velocity(self, velocity0, velocity1, velocity_total, dt: float) -> None:
        pass

    def diffuse_velocity(self, velocity0, velocity1, visc, dt: float) -> None:
        pass

    def project_velocity(self, velocity0, velocity1, dt: float) -> None:
        pass

    def add_force_field(self, density0, source, dt: float) -> None:
        return density0 + source * dt

    # p 125 transport function
    def transport_field(self, density0, density1, velocity, dt: float) -> None:
        for i, j in itertools.product(range(self.shape[0]), range(self.shape[1])):
            x = (i+0.5, j+0.5) #FIXME: Original formula is multiplied by "D". What is D?
            x_prev = self.trace_particle(x, velocity, dt)
            density1[i, j] = density0[x_prev]   # FIXME Original formula interpolates linearly? How, what and why?

    def diffuse_field(self, density0, density1, diffusion: float, dt: float) -> None:
        pass

    def dissipate_field(self, density0, density1, dissipation_rate: float, dt: float) -> None:
        pass

    # Calculates the position of a particle in the field at the last timestep. Used for advection(transport)
    # p. 125 bottom right
    def trace_particle(self, pos: (int, int), velocity, dt, steps=10) -> (int, int):
        h = -dt / steps     # Note that we use -dt to traverse the fild in the opposite direction
        curr_pos = pos
        # FIXME: Correct implementation of tuple addition/multiplication (pairwise)
        # FIXME: Find out how to handle float indices (maybe interpolate between neighboring cells?, Maybe round to next int?)
        # FIXME: velocity is a tuple of arrays
        for i in range(steps):
            k1 = h * velocity[curr_pos]
            k2 = h * velocity[curr_pos + k1/2]
            curr_pos = curr_pos + k2

    # Poisson 2d cartesian differential equation:
    # δ^2*u/δ^2x + δ^2*u/δ^2x = f(x, y)
    #
    # Finite difference equation TODO: validate
    # δ^2*u/δ^2x = (u[i+1, j] - 2u[i, j] + u[i-1, j]) / dx^2
    # δ^2*u/δ^2y = (u[i, j+1] - 2u[i, j] + u[i, j-1]) / dy^2
    def cartesian_poisson_2d(self, u, pos:(int, int)):
        pass

