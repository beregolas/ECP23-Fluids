from typing import Tuple
import numpy as np


class World2D:
    def __init__(self, size: Tuple[float, float], shape: Tuple[int, int], viscosity: float, diffusion: float,
                 dissipation_rate: float):
        # World parameters
        self.size = size
        self.shape = shape

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

    def apply_force(self, position: Tuple[int, int], force: Tuple[float, float]):
        pass

    def step(self, dt: float):
        # switch fields, older field is in [...]1
        self.velocity0, self.velocity1 = self.velocity1, self.velocity0
        self.density0, self.density1 = self.density1, self.density0
        # velocity steps
        for i in range(2):
            self.add_force_velocity(self.velocity0[i], self.force, dt)
        for i in range(2):
            self.transport_velocity(self.velocity0[i], self.velocity1[i], self.velocity0, dt)
        for i in range(2):
            self.diffuse_velocity(self.velocity0[i], self.velocity1[i], self.viscosity, dt)

        # scalar field steps
        self.add_force_field(self.density0, self.force_source, dt)
        self.transport_field(self.density1, self.density0, )

        # density steps

        pass

    def add_force_velocity(self, velocity0, force, dt: float):
        pass

    def transport_velocity(self, velocity0, velocity1, velocity_total, dt: float):
        pass

    def diffuse_velocity(self, velocity0, velocity1, visc, dt: float):
        pass

    def project_velocity(self, velocity0, velocity1, dt: float):
        pass

    def add_force_field(self, density0, source, dt: float):
        pass

    def transport_field(self, density0, density1, velocity, dt: float):
        pass

    def diffuse_field(self, density0, density1, diffusion_const: float, dt: float):
        pass

    def dissipate_field(self, density0, density1, dissipation_rate: float, dt: float):
        pass
