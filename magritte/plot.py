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
        image_nr =  -1,
        zoom     = 1.3,
        npix_x   = 300,
        npix_y   = 300,
        x_unit   = units.au,
        v_unit   = units.km/units.s
    ):
    """
    Plot channel maps of synthetic observation (image) with matplotlib.
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
    velos = (freqs - f_ij) / f_ij * constants.c.to(v_unit).value

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
        
    return interact(lambda v: figs[v], v=(0, len(figs)-1))


def image_plotly(
        model,
        image_nr =  -1,
        zoom     = 1.3,
        npix_x   = 300,
        npix_y   = 300,
        x_unit   = units.au,
        v_unit   = units.km/units.s,
        width    = 620,   # Yields approx square channel map
        height   = 540,   # Yields approx square channel map
    ):
    """
    Plot channel maps of synthetic observation (image) with plotly.
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
    velos = (freqs - f_ij) / f_ij * constants.c.to(v_unit).value

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
    
    # Create plotly plot
    fig = make_subplots(
        rows               = 1, 
        cols               = 2,
        column_widths      = [0.7, 0.3],
        horizontal_spacing = 0.05,
        subplot_titles     = ['Channel map', 'Spectrum'],
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
    x_min = np.min(xs)
    
    # Build up plot
    for f in range(nfreqs):
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
    fig.update_yaxes(
        row        = 1,
        col        = 2,
        title_text = "relative intensity",
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