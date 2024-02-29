import matplotlib.pyplot    as plt                  # Mpl plotting
import matplotlib                                   # Mpl
import plotly.graph_objects as go                   # Plotly plotting
import magritte.core        as magritte             # Core functionality
import numpy                as np                   # Data structures
import os                                           # Creating directories
import warnings                                     # Hide warnings
warnings.filterwarnings('ignore')                   # especially for yt

from matplotlib.gridspec  import GridSpec           # Plot layout
from plotly.subplots      import make_subplots      # Plotly subplots
from astropy              import constants, units   # Unit conversions
from scipy.interpolate    import griddata           # Grid interpolation
from palettable.cubehelix import cubehelix2_16      # Nice colormap
from tqdm                 import tqdm               # Progress bars
from ipywidgets           import interact           # Interactive plots
from ipywidgets.embed     import embed_minimal_html # Store interactive plots
from magritte.core        import ImageType, ImagePointPosition  # Image type, point position
from math                 import floor, ceil        # Math helper functions
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# Port matplotlib colormaps to plotly
cubehelix2_16_rgb = []
for i in range(0, 255):
    cubehelix2_16_rgb.append(
        matplotlib.colors.colorConverter.to_rgb(
            cubehelix2_16.mpl_colormap(
                matplotlib.colors.Normalize(vmin=0, vmax=255)(i)
            )
        )
    )

def matplotlib_to_plotly(cmap, pl_entries):
    """
    Converts matplotlib cmap to plotly colorscale.
    """
    h = 1.0/(pl_entries-1)
    pl_colorscale = []

    for k in range(pl_entries):
        C = list(map(np.uint8, np.array(cmap(k*h)[:3])*255))
        pl_colorscale.append([k*h, 'rgb'+str((C[0], C[1], C[2]))])
    return pl_colorscale

# Plotly compatible colorscale
cubehelix2_16_plotly = matplotlib_to_plotly(cubehelix2_16.mpl_colormap, 255)

# Plotly standard config
modeBarButtonsToRemove = ['select2d', 'lasso2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'hoverClosestCartesian', 'hoverCompareCartesian', 'zoom3d', 'pan3d', 'resetCameraDefault3d', 'resetCameraLastSave3d', 'hoverClosest3d', 'orbitRotation', 'tableRotation', 'zoomInGeo', 'zoomOutGeo', 'resetGeo', 'hoverClosestGeo', 'sendDataToCloud', 'hoverClosestGl2d', 'hoverClosestPie', 'toggleHover', 'resetViews', 'toggleSpikelines', 'resetViewMapbox']


def image_mpl(
        model,
        image_nr   =  -1,
        zoom       = 1.3,
        npix_x     = 300,
        npix_y     = 300,
        x_unit     = units.au,
        v_unit     = units.km/units.s,
        method     = 'nearest'
    ):
    """
    Create plots of the channel maps of a synthetic observation (image) with matplotlib.

    Parameters
    ----------
    model : object
        Magritte model object.
    image_nr : int
        Number of the synthetic observation to plot. (Use -1 to indicate the last one.)
    zoom : float
        Factor with which to zoom in on the middel of the image.
    npix_x : int
        Number of pixels in the image in the horizontal (x) direction.
    npix_y : int
        Number of pixels in the image in the vertical (y) direction.
    x_unit : astropy.units object
        Unit of length for the horixontal (x) axis.
    y_unit : astropy.units object
        Unit of length for the vertical (y) axis.
    method : str
        Method to interpolate the scattered intensity data onto a regular image grid.

    Returns
    -------
    None
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

    # Workaround for model images
    if (model.images[image_nr].imagePointPosition == ImagePointPosition.AllModelPoints):
    # if (False):
        # Filter imaging data originating from boundary points
        bdy_indices = np.array(model.geometry.boundary.boundary2point)
        imx = np.delete(imx, bdy_indices)
        imy = np.delete(imy, bdy_indices)
        imI = np.delete(imI, bdy_indices, axis=0)

    # Extract the number of frequency bins
    nfreqs = model.images[image_nr].nfreqs

    # Set image boundaries
    deltax = (np.max(imx) - np.min(imx))/zoom
    midx = (np.max(imx) + np.min(imx))/2.0
    deltay = (np.max(imy) - np.min(imy))/zoom
    midy = (np.max(imy) + np.min(imy))/2.0

    x_min, x_max = midx - deltax/2.0, midx + deltax/2.0
    y_min, y_max = midy - deltay/2.0, midy + deltay/2.0

    # Create image grid values
    xs = np.linspace(x_min, x_max, npix_x)
    ys = np.linspace(y_min, y_max, npix_y)

    # Extract the spectral / velocity data
    freqs = np.array(model.images[image_nr].freqs)
    f_ij  = np.mean(freqs)
    # Velocity axis definition in images is flipped compared to the direction the observer is looking
    velos = (f_ij - freqs) / f_ij * constants.c.to(v_unit).value

    # Interpolate the scattered data to an image (regular grid)
    Is = np.zeros((nfreqs))
    zs = np.zeros((nfreqs, npix_x, npix_y))
    for f in range(nfreqs):
        # Nearest neighbor interpolate scattered image data
        zs[f] = griddata(
            (imx, imy),
            imI[:,f],
            (xs[None,:], ys[:,None]),
            method=method,
            fill_value = 0.0 #for non-nearest neighbor interpolation, otherwise the ceil/floor functions will complain
        )
        Is[f] = np.sum(zs[f])
    Is = Is / np.max(Is)

    # Put zero/negative values to the smallest positive value
    zs[zs<=0.0] = np.min(zs[zs>0.0])
    # Put nan values to smallest positive value
    zs[np.isnan(zs)] = np.min(zs[zs>0.0])

    # Get the logarithm of the data (matplotlib has a hard time handling logarithmic data.)
    log_zs     = np.log10(zs)
    log_zs_min = np.min(log_zs)
    log_zs_max = np.max(log_zs)

    lzmin = ceil (log_zs_min)
    lzmax = floor(log_zs_max)

    lz_25 = ceil (log_zs_min + 0.25*(log_zs_max - log_zs_min))
    lz_50 = ceil (log_zs_min + 0.50*(log_zs_max - log_zs_min))
    lz_75 = floor(log_zs_min + 0.75*(log_zs_max - log_zs_min))

    ticks  = [lzmin, lz_25, lz_50, lz_75, lzmax]
    levels = np.linspace(log_zs_min, log_zs_max, 250)

    figs = []
    gs   = GridSpec(1,2, wspace=.1, width_ratios=[2, 1])

    for f in tqdm(reversed(range(nfreqs)), total=nfreqs):
        fig = plt.figure(dpi=300)
        ax1 = fig.add_subplot(gs[0])
        ax  = ax1.contourf(
            xs / (1.0 * x_unit).si.value,
            ys / (1.0 * x_unit).si.value,
            log_zs[f],
            cmap=cubehelix2_16.mpl_colormap,
            levels=levels
        )
        ax0 = inset_axes(
                  ax1,
                  width="100%",
                  height="5%",
                  loc='lower left',
                  bbox_to_anchor=(0, 1.025, 1, 1),
                  bbox_transform=ax1.transAxes,
                  borderpad=0
        )

        cbar = fig.colorbar(ax, cax=ax0, orientation="horizontal")
        ax0.xaxis.set_ticks_position('top')
        ax0.xaxis.set_label_position('top')
        ax0.xaxis.set_ticks         (ticks)
        ax0.xaxis.set_ticklabels    ([f'$10^{{{t}}}$' for t in ticks])

        ax1.set_aspect('equal')
        ax1.set_xlabel(f'image x [{x_unit}]', labelpad = 10)
        ax1.set_ylabel(f'image y [{x_unit}]', labelpad = 10)

        ax2 = fig.add_subplot(gs[1])
        ax2.plot(velos, Is/np.max(Is))
        ax2.yaxis.set_label_position("right")
        ax2.yaxis.tick_right()
        ax2.axvline(velos[f], c='red')
        ax2.set_xlabel(f'velocity [{v_unit}]', labelpad=10)
        asp = 2*np.diff(ax2.get_xlim())[0] / np.diff(ax2.get_ylim())[0]
        ax2.set_aspect(asp)

        if   (model.images[image_nr].imageType == ImageType.Intensity):
            ax0.set_xlabel('Intensity [W m$^{-2}$ sr$^{-1}$ Hz$^{-1}$]', labelpad=11)
            ax2.set_ylabel('Relative intensity',                         labelpad=15)
        elif (model.images[image_nr].imageType == ImageType.OpticalDepth):
            ax0.set_xlabel('Optical depth [.]',      labelpad=11)
            ax2.set_ylabel('Relative optical depth', labelpad=15)

        plt.savefig(f"{im_dir}/image_{f:0>3d}.png", bbox_inches='tight')

        figs.append(fig)

        plt.close()

    # Create a widget for plots
    widget = interact(lambda v: figs[v], v=(0, len(figs)-1))

    return widget


def image_plotly(
        model,
        image_nr   =  -1,
        zoom       = 1.3,
        npix_x     = 300,
        npix_y     = 300,
        x_unit     = units.au,
        v_unit     = units.km/units.s,
        method     = 'nearest',
        width      = 620,   # Yields approx square channel map
        height     = 540    # Yields approx square channel map
    ):
    """
    Create plots of the channel maps of a synthetic observation (image) with plotly.

    Parameters
    ----------
    model : object
        Magritte model object.
    image_nr : int
        Number of the synthetic observation to plot. (Use -1 to indicate the last one.)
    zoom : float
        Factor with which to zoom in on the middel of the image.
    npix_x : int
        Number of pixels in the image in the horizontal (x) direction.
    npix_y : int
        Number of pixels in the image in the vertical (y) direction.
    x_unit : astropy.units object
        Unit of length for the horixontal (x) axis.
    y_unit : astropy.units object
        Unit of length for the vertical (y) axis.
    method : str
        Method to interpolate the scattered intensity data onto a regular image grid.
    width : float
        Width of the resulting figure.
    height : float
        Height of the resulting figure.

    Returns
    -------
    None
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

    # Workaround for model images
    if (model.images[image_nr].imagePointPosition == ImagePointPosition.AllModelPoints):
        # Filter imaging data originating from boundary points
        bdy_indices = np.array(model.geometry.boundary.boundary2point)
        imx = np.delete(imx, bdy_indices)
        imy = np.delete(imy, bdy_indices)
        imI = np.delete(imI, bdy_indices, axis=0)

    # Extract the number of frequency bins
    nfreqs = model.images[image_nr].nfreqs

    # Set image boundaries
    deltax = (np.max(imx) - np.min(imx))/zoom
    midx = (np.max(imx) + np.min(imx))/2.0
    deltay = (np.max(imy) - np.min(imy))/zoom
    midy = (np.max(imy) + np.min(imy))/2.0

    x_min, x_max = midx - deltax/2.0, midx + deltax/2.0
    y_min, y_max = midy - deltay/2.0, midy + deltay/2.0

    # Create image grid values
    xs = np.linspace(x_min, x_max, npix_x)
    ys = np.linspace(y_min, y_max, npix_y)

    # Extract the spectral / velocity data
    freqs = np.array(model.images[image_nr].freqs)
    f_ij  = np.mean(freqs)
    # Velocity axis definition in images is flipped compared to the direction the observer is looking
    velos = (f_ij - freqs) / f_ij * constants.c.to(v_unit).value

    # Interpolate the scattered data to an image (regular grid)
    Is = np.zeros((nfreqs))
    zs = np.zeros((nfreqs, npix_x, npix_y))
    for f in range(nfreqs):
        # Nearest neighbor interpolate scattered image data
        zs[f] = griddata(
            (imx, imy),
            imI[:,f],
            (xs[None,:], ys[:,None]),
            method=method,
            fill_value = 0.0 #for non-nearest neighbor interpolation, otherwise the ceil/floor functions will complain
        )
        Is[f] = np.sum(zs[f])
    Is = Is / np.max(Is)

    # Put zero-values to the smallest non-zero value
    zs[zs<=0.0] = np.min(zs[zs>0.0])
    # Put nan values to smallest positive value
    zs[np.isnan(zs)] = np.min(zs[zs>0.0])

    # Get the logarithm of the data (matplotlib has a hard time handling logarithmic data.)
    log_zs     = np.log10(zs)
    log_zs_min = np.min(log_zs)
    log_zs_max = np.max(log_zs)

    if   (model.images[image_nr].imageType == ImageType.Intensity):
        # Create plotly plot
        fig = make_subplots(
            rows               = 1,
            cols               = 2,
            column_widths      = [0.7, 0.3],
            horizontal_spacing = 0.05,
            subplot_titles     = ['Intensity', ''],
        )
    elif (model.images[image_nr].imageType == ImageType.OpticalDepth):
        # Create plotly plot
        fig = make_subplots(
            rows               = 1,
            cols               = 2,
            column_widths      = [0.7, 0.3],
            horizontal_spacing = 0.05,
            subplot_titles     = ['Optical depth', ''],
        )


    fig.add_vrect(
        row        = 1,
        col        = 1,
        x0         = -1.0e+99,
        x1         = +1.0e+99,
        line_width = 0,
        fillcolor  = "black"
    )

    # Convert to given units
    xs = xs / (1.0 * x_unit).si.value
    ys = ys / (1.0 * x_unit).si.value

    # Convert to given units
    x_max = np.max(xs)
    x_min = np.min(xs)
    y_max = np.max(ys)
    y_min = np.min(ys)

    # Build up plot
    for f in reversed(range(nfreqs)):
        fig.add_trace(
            go.Heatmap(
                x          = xs       .astype(float),
                y          = ys       .astype(float),
                z          = log_zs[f].astype(float),
                visible    = False,
                hoverinfo  = 'none',
                zmin       = log_zs_min.astype(float),
                zmax       = log_zs_max.astype(float),
                colorscale = cubehelix2_16_plotly,
                showscale  = False
            ),
            row = 1,
            col = 1
        )

        fig.add_trace(
            go.Scatter(
                x          = (velos)        .astype(float),
                y          = (Is/np.max(Is)).astype(float),
                visible    = False,
                hoverinfo  = 'none',
                line_color = '#1f77b4',
                showlegend = False
            ),
            row = 1,
            col = 2
        )

        fig.add_trace(
            go.Scatter(
                x          = np.array([velos[f], velos[f]], dtype=float),
                y          = np.array([-1.0e+10, +1.0e+10], dtype=float),
                visible    = False,
                hoverinfo  = 'none',
                line_color = 'red',
                showlegend = False
            ),
            row = 1,
            col = 2
        )

    # Boxes around plots (more liek mpl)
    fig.update_xaxes(
        showline=True, linewidth=2, linecolor='rgba(1,1,1,1)', mirror=True, ticks='outside'
    )
    fig.update_yaxes(
        showline=True, linewidth=2, linecolor='rgba(1,1,1,1)', mirror=True, ticks='outside'
    )

    v_min = float(min(velos)) * 1.05
    v_max = float(max(velos)) * 1.05

    # Black background for channel map
    fig.add_vrect(
        row        = 1,
        col        = 1,
        x0         = 1000.0*x_min,  # large enough so you don't see edges
        x1         = 1000.0*x_max,  # large enough so you don't see edges
        line_width = 0,
        fillcolor  = "black",
        layer="below"
    )

    # Plot axes
    fig.update_xaxes(
        row        = 1,
        col        = 1,
        title_text = f'image x [{x_unit}]',
        range    = [x_min, x_max],
        showgrid = False,
        zeroline = False
    )
    fig.update_yaxes(
        row         = 1,
        col         = 1,
        title_text  = f'image y [{x_unit}]',
        scaleanchor = "x",
        scaleratio  = 1,
        showgrid    = False,
        zeroline    = False
    )
    fig.update_xaxes(
        row        = 1,
        col        = 2,
        title_text = f'velocity [{v_unit}]',
        range      = [v_min, v_max]
    )

    if   (model.images[image_nr].imageType == ImageType.Intensity):
        fig.update_yaxes(
            row        = 1,
            col        = 2,
            title_text = "Relative intensity",
            side       = 'right',
            range      = [-0.05, +1.05]
        )
    elif (model.images[image_nr].imageType == ImageType.OpticalDepth):
        fig.update_yaxes(
            row        = 1,
            col        = 2,
            title_text = "Relative opacity",
            side       = 'right',
            range      = [-0.05, +1.05]
        )

    # Subplot titles are annotations
    fig.update_annotations(
        font_size = 16,
        borderpad = 7
    )

    fig.update_layout(
        width        = width,
        height       = height,
        plot_bgcolor = 'rgba(0,0,0,0)',
        dragmode     = 'pan',
        font         = dict(family="Calibri", size=14, color='black')
    )

    # Make 3 middle traces visible
    fig.data[3*nfreqs//2  ].visible = True
    fig.data[3*nfreqs//2+1].visible = True
    fig.data[3*nfreqs//2+2].visible = True

    # Create and add slider
    steps = []
    for f in range(nfreqs):
        step = dict(
            method = "restyle",
            args = [
                {"visible": [False] * len(fig.data)},
                {"title": "Channel map: " + str(f)}
            ],
            label = ''
        )
        # Toggle f'th trace to "visible"
        step["args"][0]["visible"][3*f:3*f+3] = [True, True, True]
        steps.append(step)

    sliders = [
        dict(
            active     = nfreqs//2,
            pad        = {"t": 75},
            steps      = steps,
            tickcolor  = 'white',
            transition = {'duration': 0}
        )
    ]

    fig.update_layout(
        sliders=sliders
    )

    # Config for modebar buttons
    config = {
        "modeBarButtonsToRemove": modeBarButtonsToRemove,
        "scrollZoom": True
    }

    # Save figure as html file
    fig.write_html(f"{im_dir}/image.html", config=config)

    return fig.show(config=config)


def image_channel(
        model,
        channels: list[int],
        figure_size: tuple[int, int], # (width, height), product should be equal to len(channels)
        image_nr   =  -1,
        zoom       = 1.3,
        npix_x     = 200,
        npix_y     = 200,
        x_unit     = units.au,
        v_unit     = units.km/units.s,
        x_ticks    = None,
        y_ticks    = None,
        method     = 'nearest'
):
    """Creates a single image of channel maps plotted next to each other, with the line profile below.

    Args:
        model (_type_): The Magritte model object.
        channels (list[int]): The channels to plot
        figure_size (tuple[int, int]): Number of columns, rows to use for channel maps. product should be equal to len(channels)
        image_nr (int, optional): Image number to plot. Defaults to -1.
        zoom (float, optional): Zoom level. Defaults to 1.3.
        npix_x (int, optional): Number of pixels in x direction. Defaults to 300.
        npix_y (int, optional): Number of pixels in y direction. Defaults to 300.
        x_unit (_type_, optional): Units for position coordinate. Defaults to units.au.
        v_unit (_type_, optional): Units for velocity coordinate. Defaults to units.km/units.s.
        x_ticks (float, optional): Ticks for x-axis in x_unit units. Default to 'None', in which case the automatic matplotlib ticks are used.
        y_ticks (float, optional): Ticks for y-axis in x_unit units. Default to 'None', in which case the automatic matplotlib ticks are used.
        method (str, optional): Interpolation method. Defaults to 'nearest'.
    """
    if (len(model.images) < 1):
        print('No images in model.')
        return
    if (figure_size[0]*figure_size[1] != len(channels)):
        print('Figure size does not match number of channels.')
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

    # Workaround for model images
    if (model.images[image_nr].imagePointPosition == ImagePointPosition.AllModelPoints):
        # Filter imaging data originating from boundary points
        bdy_indices = np.array(model.geometry.boundary.boundary2point)
        imx = np.delete(imx, bdy_indices)
        imy = np.delete(imy, bdy_indices)
        imI = np.delete(imI, bdy_indices, axis=0)

    # Extract the number of frequency bins
    nfreqs = model.images[image_nr].nfreqs

    # Set image boundaries
    deltax = (np.max(imx) - np.min(imx))/zoom
    midx = (np.max(imx) + np.min(imx))/2.0
    deltay = (np.max(imy) - np.min(imy))/zoom
    midy = (np.max(imy) + np.min(imy))/2.0

    x_min, x_max = midx - deltax/2.0, midx + deltax/2.0
    y_min, y_max = midy - deltay/2.0, midy + deltay/2.0

    # Create image grid values
    xs = np.linspace(x_min, x_max, npix_x)
    ys = np.linspace(y_min, y_max, npix_y)

    # Extract the spectral / velocity data
    freqs = np.array(model.images[image_nr].freqs)
    f_ij  = np.mean(freqs)
    velos = (f_ij - freqs) / f_ij * constants.c.to(v_unit).value

    # Interpolate the scattered data to an image (regular grid)
    Is = np.zeros((nfreqs))
    zs = np.zeros((nfreqs, npix_x, npix_y))
    for f in range(nfreqs):
        # Nearest neighbor interpolate scattered image data
        zs[f] = griddata(
            (imx, imy),
            imI[:,f],
            (xs[None,:], ys[:,None]),
            method=method,
            fill_value = 0.0 #for non-nearest neighbor interpolation, otherwise the ceil/floor functions will complain
        )
        Is[f] = np.sum(zs[f])
    Is = Is / np.max(Is)

    # Put zero/negative values to the smallest positive value
    zs[zs<=0.0] = np.min(zs[zs>0.0])
    # Put nan values to smallest positive value
    zs[np.isnan(zs)] = np.min(zs[zs>0.0])

    # Get the logarithm of the data (matplotlib has a hard time handling logarithmic data.)
    log_zs     = np.log10(zs)
    log_zs_min = np.min(log_zs)
    log_zs_max = np.max(log_zs)

    lzmin = ceil (log_zs_min)
    lzmax = floor(log_zs_max)

    lz_25 = ceil (log_zs_min + 0.25*(log_zs_max - log_zs_min))
    lz_50 = ceil (log_zs_min + 0.50*(log_zs_max - log_zs_min))
    lz_75 = floor(log_zs_min + 0.75*(log_zs_max - log_zs_min))

    ticks  = [lzmin, lz_25, lz_50, lz_75, lzmax]
    levels = np.linspace(log_zs_min, log_zs_max, 250)

    colorbarwidth = 0.05
    lineprofileheight = 1.0

    fig = plt.figure(figsize = (1.5*(figure_size[1]+colorbarwidth),1.5*(figure_size[0]+lineprofileheight)), dpi=300)
    colorbarsplit = fig.subfigures(1,2, width_ratios=[1, colorbarwidth], wspace=0)
    colorbarsplitrelevant = colorbarsplit[1].subfigures(2,1, height_ratios=[figure_size[1],lineprofileheight], hspace = 0)
    colorbarax = colorbarsplitrelevant[0].subplots(1,1)
    subfigs = colorbarsplit[0].subfigures(2, 1, height_ratios=[figure_size[0], 1], hspace = 0.0)
    channelmapsfig = subfigs[0].subplots(figure_size[0], figure_size[1], sharex=True, sharey=True)
    subfigs[0].subplots_adjust(wspace = 0.05, hspace = 0.05)
    if type(channelmapsfig)==np.ndarray:
        channelmapsfig = channelmapsfig.flatten()
    lineprofilefig = subfigs[1].subplots(1, 1)

    for f in range(len(channels)):
        im = channelmapsfig[f].contourf(
            xs / (1.0 * x_unit).si.value,
            ys / (1.0 * x_unit).si.value,
            log_zs[channels[f]],
            cmap=cubehelix2_16.mpl_colormap,
            levels=levels
        )
        lineprofilefig.axvline(velos[channels[f]], linestyle=':', color='r')
        channelmapsfig[f].text(0.07,0.07,f'{velos[channels[f]]:.2f} '+str(v_unit), transform = channelmapsfig[f].transAxes, color='white')
        if x_ticks is not None:
            channelmapsfig[f].xaxis.set_ticks(x_ticks)
        if y_ticks is not None:
            channelmapsfig[f].yaxis.set_ticks(y_ticks)

    #now also plot the line profile
    lineprofilefig.plot(velos, Is/max(Is))
    #and some indicators of the line profile
    lineprofilefig.axvline(min(velos[channels]), color='r')
    lineprofilefig.axvline(max(velos[channels]), color='r')
    subfigs[1].supxlabel(f'velocity [{str(v_unit)}]', y = -0.20)

    #also add the figure axes
    subfigs[0].supxlabel(f'X [{str(x_unit)}]', y = 0.03)
    subfigs[0].supylabel(f'Y [{str(x_unit)}]')
    subfigs[0].suptitle("Intensity [W m$^{-2}$ sr$^{-1}$ Hz$^{-1}$]", y = 0.95)
    subfigs[1].suptitle("Relative intensity", y = 1.01)
    # Fix the colorbar ticks, as we plot the logarithm of the data
    cbar = fig.colorbar(im, cax=colorbarax)
    cbarticks = np.arange(lzmin, lzmax+1)
    cbar.ax.set_yticks(cbarticks)
    cbar.ax.set_yticklabels([f'{10**float(str(t)):.0e}' for t in cbarticks])

    return fig


def plot_velocity_1D(model, xscale='log', yscale='linear'):
    """
    Plot the velocity profile of the model (in 1D, radially).

    Parameters
    ----------
    model : object
        Magritte model object.
    xscale : str
        Scale of the xaxis ("linear", "log", "symlog", "logit", ...)
    yscale : str
        Scale of the yaxis ("linear", "log", "symlog", "logit", ...)

    Returns
    -------
    None
    """
    rs = np.linalg.norm(model.geometry.points.position, axis=1)
    vs = np.linalg.norm(model.geometry.points.velocity, axis=1)

    fig = plt.figure(dpi=300)
    plt.plot  (rs, vs*constants.c.si.value)
    plt.xlabel('radius [m]',     labelpad=10)
    plt.ylabel('velocity [m/s]', labelpad=10)
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.show  ()


def plot_temperature_1D(model, xscale='log', yscale='linear'):
    """
    Plot the temperature profile of the model (in 1D, radially).

    Parameters
    ----------
    model : object
        Magritte model object.
    xscale : str
        Scale of the xaxis ("linear", "log", "symlog", "logit", ...)
    yscale : str
        Scale of the yaxis ("linear", "log", "symlog", "logit", ...)

    Returns
    -------
    None
    """
    rs   = np.linalg.norm(model.geometry.points.position, axis=1)
    temp = np.array      (model.thermodynamics.temperature.gas)

    fig = plt.figure(dpi=300)
    plt.plot  (rs, temp)
    plt.xlabel('radius [m]',      labelpad=10)
    plt.ylabel('temperature [K]', labelpad=10)
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.show  ()


def plot_turbulence_1D(model, xscale='log', yscale='linear'):
    """
    Plot the (micro) turbulence profile of the model (in 1D, radially).

    Parameters
    ----------
    model : object
        Magritte model object.
    xscale : str
        Scale of the xaxis ("linear", "log", "symlog", "logit", ...)
    yscale : str
        Scale of the yaxis ("linear", "log", "symlog", "logit", ...)

    Returns
    -------
    None
    """
    rs     = np.linalg.norm(model.geometry.points.position, axis=1)
    vturb2 = np.array      (model.thermodynamics.turbulence.vturb2)

    fig = plt.figure(dpi=300)
    plt.plot  (rs, np.sqrt(vturb2)*constants.c.si.value)
    plt.xlabel('radius [m]',       labelpad=10)
    plt.ylabel('turbulence [m/s]', labelpad=10)
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.show  ()


def plot_number_densities_1D(model, xscale='log', yscale='log'):
    """
    Plot the number densities of all species in the model (in 1D, radially).

    Parameters
    ----------
    model : object
        Magritte model object.
    xscale : str
        Scale of the xaxis ("linear", "log", "symlog", "logit", ...)
    yscale : str
        Scale of the yaxis ("linear", "log", "symlog", "logit", ...)

    Returns
    -------
    None
    """
    rs   = np.linalg.norm(model.geometry.points.position, axis=1)
    abns = np.array      (model.chemistry.species.abundance)
    syms = np.array      (model.chemistry.species.symbol)

    for s in range(1, model.parameters.nspecs()-2):
        fig = plt.figure(dpi=300)
        plt.plot  (rs, abns[:,s])
        plt.xlabel('radius [m]',                               labelpad=10)
        plt.ylabel(f'{syms[s]} number density [m$^{{{-3}}}$]', labelpad=10)
        plt.xscale(xscale)
        plt.yscale(yscale)
        plt.show  ()


def plot_populations_1D(model, lev_max=7, xscale='log', yscale='log'):
    """
    Plot the relative populations in the model (in 1D, radially).

    Parameters
    ----------
    model : object
        Magritte model object.
    lev_max : int
        Number of levels to plot.
    xscale : str
        Scale of the xaxis ("linear", "log", "symlog", "logit", ...)
    yscale : str
        Scale of the yaxis ("linear", "log", "symlog", "logit", ...)

    Returns
    -------
    None
    """
    rs      = np.linalg.norm(model.geometry.points.position, axis=1)
    npoints = model.parameters.npoints()

    for lspec in model.lines.lineProducingSpecies:
        nlev     = lspec.linedata.nlev
        pops     = np.array(lspec.population).reshape((npoints,nlev))
        pops_tot = np.array(lspec.population_tot)

        plt.figure(dpi=300)

        for i in range(min([lev_max, nlev])):
            plt.plot  (rs, pops[:,i]/pops_tot, label=f'i={i}')
            plt.ylabel('fractional level populations [.]', labelpad=10)
            plt.xlabel('radius [m]',                       labelpad=10)
            plt.xscale(xscale)
            plt.yscale(yscale)
            plt.legend()
        plt.show()
