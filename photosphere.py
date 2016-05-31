from optparse import OptionParser
import os.path
import tempfile
import subprocess
from PIL import Image
import math

def write_proj_params(input_path, span, baseline):
    im = Image.open(input_path)
    width, height = im.size

    radius = width / math.radians(span)
    px_above = height * (1. - baseline)
    px_below = height * baseline
    top_angle = math.degrees(math.atan(px_above / radius))
    bottom_angle = math.degrees(math.atan(px_below / radius))
    proj_height = int(round((top_angle + bottom_angle) * width / span))
    top_offset = int(round((90. - top_angle) * width / span))
    print 'vfov: %s' % (top_angle + bottom_angle)

    proj_path = tempfile.mktemp('.pto')
    with open(proj_path, 'w') as f:
        # TODO support offset baseline
        f.write('p f2 w%d h%d v%.3f\n' % (width, proj_height, span))
        f.write('m g2.2 i0\n')
        f.write('i w%d h%d f1 v%.3f r0 p0 y0 n"%s"\n' % (width, height, span, input_path))
    return proj_path, top_offset

def reproject(input_pano, span, baseline, bgcolor):
    with_bg = tempfile.mktemp('.png')
    subprocess.check_call(['convert', input_pano, '-background', bgcolor, '-flatten', with_bg])

    proj_path, crop_top = write_proj_params(with_bg, span, baseline)
    proj_out = tempfile.mktemp()
    subprocess.check_call(['nona', '-o', proj_out, proj_path])
    proj_out += '0000.tif'
    assert os.path.exists(proj_out)

    return proj_out, crop_top

def make_photosphere(final_out, proj_out, crop_top, span, left):
    subprocess.check_call(['convert', '-quality', '98', proj_out, final_out])

    im = Image.open(final_out)
    width, height = im.size
    full_width = int(round(width / span * 360.))
    full_height = int(round(full_width * .5))
    crop_left = int(round(full_width * .5 * (1. - span / 360.)))

    xmp_vars = {
        'UsePanoramaViewer': True,
        'ProjectionType': 'equirectangular',
        'CroppedAreaImageWidthPixels': width,
        'CroppedAreaImageHeightPixels': height,
        'FullPanoWidthPixels': full_width,
        'FullPanoHeightPixels': full_height,
        'CroppedAreaLeftPixels': crop_left,
        'CroppedAreaTopPixels': crop_top,
        'PoseHeadingDegrees': left + .5 * span,
    }
    xmp_path = tempfile.mktemp()
    with open(xmp_path, 'w') as f:
        f.write('reg GPano http://ns.google.com/photos/1.0/panorama/\n')
        for k, v in sorted(xmp_vars.items()):
            f.write('set Xmp.GPano.%s XmpText %s\n' % (k, v))
    subprocess.check_call(['exiv2', '-m', xmp_path, final_out])

def geo_tag(lat, lon, path):
    # offset slightly from pole so bearings still work as expected
    POLE_MARGIN = 1e-5
    POLE_THRESH_LAT = 90. - POLE_MARGIN
    if lat > POLE_THRESH_LAT:
        lat, lon = POLE_THRESH_LAT, 180.
    elif lat < -POLE_THRESH_LAT:
        lat, lon = -POLE_THRESH_LAT, 0.
        
    ns = 'N' if lat >= 0 else 'S'
    ew = 'E' if lon >= 0 else 'W'
    args = {
        'GPSLatitudeRef': ns,
        'GPSLatitude': abs(lat),
        'GPSLongitudeRef': ew,
        'GPSLongitude': abs(lon),
    }
    subprocess.check_call(['exiftool', '-overwrite_original'] + ['-%s=%s' % (k, v) for k, v in args.iteritems()] + [path])

def load_params(path):
    with open(path) as f:
        lines = filter(None, f.read().split('\n'))
    return dict(ln.split(': ', 1) for ln in lines)

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option("-o", "--output", dest="output",
                      help="output path of photosphere")
    parser.add_option("-p", "--params", dest="params",
                      help="path of params file")
    parser.add_option("-b", "--baseline", dest="baseline", type="float", default=.5,
                      help="location of 'equator' of image, from 0 (bottom) to 1 (top)")
    parser.add_option("--background", dest="bgcolor", default="#ffffff", help="background color")

    (options, args) = parser.parse_args()
    input_pano = args[0]
    assert options.baseline == .5, 'can\'t support different baselines yet'

    params_path = options.params or (os.path.splitext(input_pano)[0] + '.params')
    assert os.path.exists(params_path)

    output_path = options.output or (os.path.splitext(input_pano)[0] + '.photosphere.jpg')

    params = load_params(params_path)
    lat, lon = map(float, params['origin'].split(','))
    left = float(params['lonleft'])
    span = float(params['lonspan'])
    assert span <= 360., 'cannot make photosphere with lonspan > 360'

    proj_out, crop_top = reproject(input_pano, span, options.baseline, options.bgcolor)
    make_photosphere(output_path, proj_out, crop_top, span, left)
    geo_tag(lat, lon, output_path)

    print output_path
