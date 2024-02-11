from typing import Tuple
import numpy as np


class World2D:
    def __init__(self, size: Tuple(float, float), shape: Tuple(int, int)):
        self.size = size
        self.data = np.zeros(shape)

    pass
