


def dump_ix(root, path=None):
    import pickle
    if not path:
        path = tempfile.mktemp()
    with open(path, 'w') as f:
        pickle.dump(root, f)
    return path

def load_ix(path):
    import pickle
    with open(path) as f:
        return pickle.load(f)

def proc_index(node, handler, depth=0, tile=(0, 0)):
    if not node:
        return
    if node.children:
        for child, subtile in zip(node.children, [(2*tile[0] + i, 2*tile[1] + j) for j in xrange(2) for i in xrange(2)]):
            proc_index(child, handler, depth + 1, subtile)
    else:
        handler(node, depth, tile)


def render(ix, depth):
    dim = 2**depth
    BITMAP = [[0 for i in xrange(dim)] for j in xrange(dim)]
    def paint(val, z, (x, y)):
        c = {None: 0, 0: 128, 1: 255}[val.land]
        celldim = 2**(depth - z)
        for i in xrange(celldim):
            for j in xrange(celldim):
                BITMAP[y*celldim + j][x*celldim + i] = c
    proc_index(ix, paint)

    import os
    import tempfile
    raw = tempfile.mktemp('.grey')
    img = tempfile.mktemp('.png')

    with open(raw, 'w') as f:
        f.write(''.join(''.join(chr(col) for col in row) for row in BITMAP))
    os.popen('convert -size %dx%d -depth 8 gray:%s %s' % (dim, dim, raw, img))
    os.popen('gnome-open %s' % img)

if __name__ == "__main__":
    max_depth = int(sys.argv[1])

    coast = load(COAST_DATA)
    #admin = load(ADMIN_DATA, True)
    print 'loaded'
    ix = build_index(coast, GLOBAL, max_depth)
    print 'index built'
    print dump_ix(ix)
