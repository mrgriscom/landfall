


def display(poly):
    from matplotlib import pyplot
    from matplotlib.patches import Circle
    from shapely.geometry import Polygon
    from descartes.patch import PolygonPatch

    def plot_coords(ax, ob):
        x, y = ob.xy
        ax.plot(x, y, 'o', color='#999999', zorder=1)
    
    fig = pyplot.figure(1, figsize=[5, 5], dpi=90)

    ax = fig.add_subplot(111)

    for poly in explode_poly(poly):
        for inter in poly.interiors:
            plot_coords(ax, inter)
        plot_coords(ax, poly.exterior)

        patch = PolygonPatch(poly, facecolor='#888888', edgecolor='#222222')
        ax.add_patch(patch)

    #ax.set_aspect(1)
    pyplot.show()
