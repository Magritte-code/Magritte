import matplotlib.pyplot as plt                     # Plotting
import magritte.core     as magritte                # Core functionality
import numpy             as np                      # Data structures
import os                                           # Creating directories

from matplotlib.gridspec  import GridSpec           # Plot layout
from astropy              import constants, units   # Unit conversions
from scipy.interpolate    import griddata           # Grid interpolation
from palettable.cubehelix import cubehelix2_16      # Nice colormap
from tqdm                 import tqdm               # Progress bars
from ipywidgets           import interact           # Interactive plots


def image(
        model,
        image_nr =  -1,
        zoom     = 1.3,
        npix_x   = 300,
        npix_y   = 300,
        x_unit   = units.au,
        v_unit   = units.km/units.s
    ):
    """
    Plot image.
    """
    
    # Check if there are images
    if (len(model.images) < 1):
        print('No images in model.')
        return
    
    # Get path of image directory
    im_dir = os.path.dirname(os.path.abspath(model.parameters.model_name())) + '/images/'
    
    # If no image directory exists yet
    if not os.path.exists(im_dir):
        # Create image directory
        os.makedirs(im_dir)
        print('Created image directory:', im_dir)
    
    # Extract data of last image
    imx = np.array(model.images[image_nr].ImX)
    imy = np.array(model.images[image_nr].ImY)
    imI = np.array(model.images[image_nr].I)
    imv = np.array(model.radiation.frequencies.nu)[0]

    # Extract the number of frequency bins
    nfreqs = model.parameters.nfreqs()
    
    # Set image boundaries
    x_min, x_max = np.min(imx)/zoom, np.max(imx)/zoom
    y_min, y_max = np.min(imy)/zoom, np.max(imy)/zoom

    # Create image grid values
    xs = np.linspace(x_min, x_max, npix_x)
    ys = np.linspace(y_min, y_max, npix_y)
    
    # Extract the spectral / velocity data
    freqs = np.array(model.radiation.frequencies.nu)[0]
    f_ij  = np.mean(freqs)
    velos = (freqs - f_ij) / f_ij * constants.c.to(v_unit).value   # [km/s]

    # Interpolate the scattered data to an image (regular grid)
    Is = np.zeros((nfreqs))
    zs = np.zeros((nfreqs, npix_x, npix_y))
    for f in range(nfreqs):
        # Nearest neighbor interpolate scattered image data
        zs[f] = griddata(
            (imx, imy),
            imI[:,f],
            (xs[None,:], ys[:,None]),
            method='nearest'
        )
        Is[f] = np.sum(zs[f])
    Is = Is / np.max(Is)

    # Get the logarithm of the data (matplotlib has a hard time handling logarithmic data.)
    log_zs     = np.log(zs)
    log_zs_min = np.min(log_zs)
    log_zs_max = np.max(log_zs)
    
    figs = []
    gs   = GridSpec(1,2, wspace=.1, width_ratios=[2, 1])

    for f in tqdm(range(nfreqs)):
        fig = plt.figure(dpi=300)
        ax1 = fig.add_subplot(gs[0])
        ax1.contourf(
            xs / (1.0 * x_unit).si.value,
            ys / (1.0 * x_unit).si.value,
            log_zs[f],
            cmap=cubehelix2_16.mpl_colormap,
            vmin=log_zs_min,
            vmax=log_zs_max,
            levels=250
        )
        ax1.set_aspect('equal')
        ax1.set_xlabel(f'image x [{x_unit}]', labelpad = 10)
        ax1.set_ylabel(f'image y [{x_unit}]', labelpad = 10)
    
        ax2 = fig.add_subplot(gs[1])
        ax2.plot(velos, Is/np.max(Is))
        ax2.yaxis.set_label_position("right")
        ax2.yaxis.tick_right()
        ax2.axvline(velos[f], c='red')
        ax2.set_ylabel(f'Relative intensity',  labelpad=15)
        ax2.set_xlabel(f'velocity [{v_unit}]', labelpad=10)
        asp = 2*np.diff(ax2.get_xlim())[0] / np.diff(ax2.get_ylim())[0]
        ax2.set_aspect(asp)

        plt.savefig(f"{im_dir}/image_{f:0>3d}.png", bbox_inches='tight')
    
        figs.append(fig)
    
        plt.close()
        
    return figs