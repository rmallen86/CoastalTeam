"""Make gif from images."""
try:
    import imageio
    from pygifsicle import optimize
except Exception:
    raise ImportError('Required libraries could not be found: \n'
                      'imageio and pygifsicle')

# first make the gif - set range based on number of output images
images = []
for i in range(0, 30):
    fname = str(i) + '.png'
    images.append(imageio.imread(fname))

imageio.mimsave('demo.gif', images, duration=0.25)

# then optimize it
optimize('demo.gif')
