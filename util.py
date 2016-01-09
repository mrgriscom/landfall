
def display_polygon(poly):
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

def render_index(decomp, depth):
    import process_data as pd

    full, partial = decomp

    dim = 2**depth
    BITMAP = [[0 for i in xrange(dim)] for j in xrange(dim)]
    def paint(ix, val):
        z, x, y = 0, 0, 0
        for k in ix:
            k = int(k)
            z += 1
            x = 2*x + k % 2
            y = 2*y + k // 2
        c = {0: 0, 1: 128, 2: 255}[val]
        celldim = 2**(depth - z)
        for i in xrange(celldim):
            for j in xrange(celldim):
                BITMAP[y*celldim + j][x*celldim + i] = c
    for ix in full:
        paint(ix, pd.OVERLAP_FULL)
    for ix in partial.keys():
        paint(ix, pd.OVERLAP_PARTIAL)

    import os
    import tempfile
    raw = tempfile.mktemp('.grey')
    img = tempfile.mktemp('.png')

    with open(raw, 'w') as f:
        f.write(''.join(''.join(chr(col) for col in row) for row in BITMAP))
    os.popen('convert -size %dx%d -depth 8 gray:%s %s' % (dim, dim, raw, img))
    os.popen('gnome-open %s' % img)

