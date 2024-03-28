import numpy as np
import imageio
import matplotlib.pyplot as plt
import time
import matplotlib.animation as animation

from world import World2D

if __name__ == "__main__":

    #TODO: get input image
    
    shape = (15, 15)
    world = World2D(shape, shape, 0.0005, 0.0, 0.01)
    # world.density = np.random.randn(shape[0], shape[1])
    # world.density[5, 5] = 1
    # world.density[12, 8] = 1
    world.density[6:7, 2] = 1
    world.velocity[0] = np.random.randn(shape[0], shape[1])
    world.velocity[1] = np.random.randn(shape[0], shape[1])
    world.velocity[0] = np.zeros(shape)
    world.velocity[1] = np.zeros(shape)
    world.velocity[0, 6:8, 1:5] = 10
    world.velocity[1, 10:12, 2:4] = -5

    x = np.linspace(0, shape[0], shape[0])
    y = np.linspace(0, shape[1], shape[1])
    X, Y = np.meshgrid(x, y)

    fig, (ax1, ax2) = plt.subplots(1, 2)
    quiv = ax1.quiver(X, Y, world.velocity[0], world.velocity[1])
    cmesh = ax2.pcolormesh(X, Y, world.density)

    def animate_velocity(num):
        world.step(.1)

        quiv.set_UVC(world.velocity[0], world.velocity[1])
        return quiv

    def animate_density(num):
        cmesh.set_array(world.density)
        return cmesh

    anim1 = animation.FuncAnimation(fig, animate_velocity, interval=1000, blit=False)
    anim2 = animation.FuncAnimation(fig, animate_density, interval=1000, blit=False)
    plt.show()

    # TODO: output image with imageio
    # maybe one image per step for an animation(gif)

    pass
