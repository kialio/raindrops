"""
Fermi Rain Animation

Author: Jeremy S. Perkins (Fermi mods)

Original Author: Nicolas P. Rougier

"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from astropy.io import fits
import matplotlib as mpl
from astropy.time import Time

# Open all of the needed files and setup some things
hdu = fits.open('prob_ft1_00.fits')
mask = hdu[1].data['3FGL J0509.4+0541'] > 0.85
RA = (hdu[1].data['RA'][mask]-np.average(hdu[1].data['RA'][mask]) + 1)/2.
DEC = (hdu[1].data['DEC'][mask]-np.average(hdu[1].data['DEC'][mask]) + 1)/2.
TIMEMET = hdu[1].data['TIME'][mask]
TIME = TIMEMET - np.min(TIMEMET)
ENERGY = hdu[1].data['ENERGY'][mask]

# Calculate the dates
MJDREFI = hdu[0].header['MJDREFI']
MJDREFF = hdu[0].header['MJDREFF']
TIMEMJD = hdu[1].data['TIME'] / 86400. + MJDREFI + MJDREFF

# This is the width (in seconds) of each frame
WIDTH = 60. * 10.
FRAMES = np.int(TIME[-1] / WIDTH)


# This is the setup for the colormap of the raindrops
cm_array = np.genfromtxt('colormap2.txt', delimiter=',')
cm = mpl.colors.ListedColormap(cm_array/255)
cm_norm = mpl.colors.Normalize(vmin=np.log10(ENERGY.min()),
                               vmax=np.log10(ENERGY.max()))
cmap = mpl.cm.ScalarMappable(norm=cm_norm, cmap=cm)

# Create new Figure and an Axes which fills it.
my_dpi = 72
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi),
                 dpi=my_dpi, facecolor='black')
ax = fig.add_axes([0, 0, 1, 1], frameon=False)
ax.set_xlim(0, 1), ax.set_xticks([])
ax.set_ylim(0, 1), ax.set_yticks([])

# Create rain data
n_drops = len(TIME)
rain_drops = np.zeros(n_drops, dtype=[('position', float, 2),
                                      ('size',     float, 1),
                                      ('growth',   float, 1),
                                      ('color',    float, 4)])

rain_drops['position'] = zip(RA, DEC)
rain_drops['growth'] = ENERGY/10.
rain_drops['color'] = cmap.to_rgba(np.log10(ENERGY))
rain_drops['color'][:, 3] = 0


# Construct the scatter which we will update during animation
# as the raindrops develop.
scat = ax.scatter(rain_drops['position'][:, 0], rain_drops['position'][:, 1],
                  s=rain_drops['size'], lw=1.0, edgecolors=rain_drops['color'],
                  facecolors='none')


# Construct the histogram for the counts map
H, xedges, yedges = np.histogram2d(rain_drops['position'][:, 0],
                                   rain_drops['position'][:, 1], bins=15)
im = plt.imshow(H, interpolation='none', origin='low',
                extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                cmap='bone', vmin=5, vmax=H.max()*2.0)

# Uncomment this out if you want the energy text printed
# E_txt = [ax.text(rain_drops['position'][i, 0],
#                 rain_drops['position'][i,1], Es,
#                 color=rain_drops['color'][0])
#         for i,Es in enumerate(E_str)]


# This is the inital time step
TIMEMJD = (TIMEMET[0])/86400. + MJDREFI + MJDREFF
t = Time(TIMEMJD, format='mjd')
txt_time = ax.text(0.8, 0.9, t.datetime.strftime("%d %b %Y %H:%M"))
txt_time.set_color((1, 1, 1, 1))

# This is the line showing the scale
txt_scale = ax.text(0.13, 0.07, r'0.2$^\circ$')
txt_scale.set_color((1, 1, 1, 1))
line, = ax.plot([0.1, 0.2], [0.1, 0.1], 'k-')
line.set_color((1, 1, 1, 1))


def update(frame_number):
    # Setup the masks for the events we want
    mask = TIME < (frame_number + 1) * WIDTH
    new = (TIME > (frame_number) * WIDTH) & (TIME < (frame_number + 1) * WIDTH)

    # Figure out the time step of the current frame and update the time string
    TIMEMJD = (TIMEMET[0] + (frame_number + 1) * WIDTH)/86400.
    TIMEMJD += MJDREFI + MJDREFF
    t = Time(TIMEMJD, format='mjd')
    txt_time.set_text(t.datetime.strftime("%d %b %Y %H:%M"))

    # Make all colors more transparent as time progresses.
    rain_drops['color'][:, 3][mask] -= 1.0/50
    rain_drops['color'][:, 3][mask] = np.clip(rain_drops['color'][:, 3][mask],
                                              0, 1)

    # Fade out the energy texts
    # [E_txt[i].set_color(color) for i, color
    #  in enumerate(rain_drops['color'][:15])]

    # Fade the scale line out if you want
    # line_color = np.clip(line.get_color()[3] - 1.0/100,0,1)
    # line.set_color((1,1,1,line_color))
    # txt_scale.set_color((1,1,1,line_color))

    # Make all circles bigger and get new circles
    rain_drops['size'][mask] += rain_drops['growth'][mask]
    rain_drops['size'][new] = 5
    rain_drops['color'][new] = rain_drops['color'][new] + (0, 0, 0, 1)
    rain_drops['color'][:, 3][new] = np.clip(rain_drops['color'][:, 3][new],
                                             0, 1)

    # Update the histogram with a new frame
    hist = np.histogram2d(rain_drops['position'][mask][:, 0],
                          rain_drops['position'][mask][:, 1])
    im.set_data(hist[0])

    # Update the scatter collection, with the new colors, sizes and positions.
    scat.set_edgecolors(rain_drops['color'])
    scat.set_sizes(rain_drops['size'])
    scat.set_offsets(rain_drops['position'])

    
# Construct the animation, using the update function as the animation
# director.
animation = FuncAnimation(fig, update, frames=FRAMES,
                          interval=50, blit=False, repeat=False)
# animation.save('3C279_scale.gif', writer='imagemagick',
#                fps=30,savefig_kwargs={'facecolor':'black'})
animation.save('TXS_0506_full.mp4', fps=30,
               extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p'],
               savefig_kwargs={'facecolor': 'black'})
plt.show()
