import numpy as np
import imageio
import matplotlib.pyplot as plt
import time
import matplotlib.animation as animation
import matplotlib.colors as colors

from world import World2D

if __name__ == "__main__":
    # =========|CONSTANT PARAMETERS |====================
    shape = (20, 20)        # Size of the animation
    fps = 1                 # Frames per second. Each frame is on computation step. Increasing this increases the load
    animation_speed = 1     # Ratio of animation tim to real time

    viscosity = 0.0005
    diffusion = 0.01
    dissipation_rate = 0.01
    # ===================================================

    world = World2D(shape, shape, viscosity, diffusion, dissipation_rate)

    # =========|INITIAL CONDITIONS |=====================
    world.density[6:7, 2] = 1
    world.velocity[0, 6:8, 1:5] = 10
    world.velocity[1, 10:12, 2:4] = -5
    # ===================================================

    x = np.linspace(0, shape[0], shape[0])
    y = np.linspace(0, shape[1], shape[1])
    X, Y = np.meshgrid(x, y)

    plt.style.use("dark_background")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), dpi=200)
    ax1.set_box_aspect(1)
    ax1.set_title('Velocity')
    ax2.set_box_aspect(1)
    ax2.set_title('Density')
    quiv = ax1.quiver(X, Y, world.velocity[0], world.velocity[1], color="white")

    shadeopts = {'cmap': 'inferno', 'shading': 'gouraud'}

    def create_contourf(ax):
        return ax.contourf(X, Y, world.density, norm=colors.SymLogNorm(linthresh=0.001, linscale=1), levels=150,
                           **shadeopts)

    def animate_velocity(num):
        world.step(fps * animation_speed)

        quiv.set_UVC(world.velocity[0], world.velocity[1])
        return quiv


    def animate_density(num):
        ax2.clear()
        return create_contourf(ax2)


    anim1 = animation.FuncAnimation(fig, animate_velocity, interval=100 / fps, blit=False)
    anim2 = animation.FuncAnimation(fig, animate_density, interval=100 / fps, blit=False)
    plt.show()
