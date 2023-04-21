import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob

img_paths = sorted(glob.glob('plot/*.png'))
images = [plt.imread(path) for path in img_paths]

def animate(i):
    return plt.imshow(images[i])

plt.axis('off')

fps = 5  # Desired frame rate
interval = 1000 / fps  # Compute the interval

anim = animation.FuncAnimation(plt.gcf(), animate, frames=len(images), interval=interval)

anim.save('animation.gif', writer='imagemagick', fps=fps)

