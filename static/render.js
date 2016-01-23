EARTH_MEAN_RAD = 6371009.;
EARTH_FARTHEST = Math.PI * EARTH_MEAN_RAD;
EARTH_CIRCUMF = 2 * Math.PI * EARTH_MEAN_RAD;

function init() {
    // DATA
    // PARAMS
    // ADMIN_INFO

    var canvas = render();
    setupMap(canvas);
}

function render() {
    var width = PARAMS.dim[0];
    var height = PARAMS.dim[1];

    var c = mk_canvas(DATA.postings.length, height)
    var ctx = c.context;
    //$('body').append(c.canvas);

    var logmin = Math.log(DATA.min_dist); // TODO logmin for sat
    var logmax = Math.log(EARTH_CIRCUMF);



    var disty = function(dist) {
        var logdist = Math.log(dist);
        var k = Math.max((logdist - logmin) / (logmax - logmin), 0);
        var y = height * k;
        return {k: k, y: y};
    }

    // land
    for (var i = 0; i < DATA.postings.length; i++) {
        var dist = DATA.postings[i][0];
        var admin = admin_postings[i];
        if (dist >= 0) {
            var _d = disty(dist);
            // TODO randomize coloring
            var C = PARAMS.colors[PARAMS.hues[colors[admin]]][Math.floor(_d.k * (PARAMS.color_steps - 1))];
            ctx.fillStyle = C;
            ctx.fillRect(i, _d.y, 1, height - _d.y);
        }
    }

    var c2 = mk_canvas(DATA.postings.length, height)
    var ctx2 = c2.context;
    //$('body').append(c2.canvas);

    // creases
    for (var i = 0; i < DATA.postings.length; i++) {
        var dist0 = DATA.postings[i][0];        
        var dist1 = DATA.postings[(i + 1) % DATA.postings.length][0];
        var diff = Math.abs(dist0 - dist1);
        var closer = Math.min(dist0, dist1);
        var farther = Math.max(dist0, dist1);
        var quantum = closer * Math.PI / 180. * DATA.res;
        if (diff >= 10 * quantum) {
            var logdiff = Math.abs(Math.log(dist0) - Math.log(dist1));
            var k = (Math.log(farther) - logmin) / (logmax - logmin);
            var y = height * k;
            
            var pixelw = DATA.postings.length / width;
            var weight = 6 * logdiff / (logmax - logmin);
            var weight2 = Math.min(diff / 5e5, .3);
            weight = Math.max(weight, weight2);
            var _w = pixelw * Math.min(weight, 1.) / 3;
            var x0 = i+1 - (dist0 > dist1 ? _w : 0);
            ctx2.fillStyle = '#000';
            ctx2.fillRect(x0, y, _w, height);
        }

    }

    var downsize = function(src) {
        var cdest = mk_canvas(width, height);
        pica.WW = false;
        pica.resizeCanvas(src.canvas, cdest.canvas, {
            quality: 3,
            alpha: true,
            unsharpAmount: 0,
            unsharpThreshold: 0,
            transferable: false //true
        }, function (err) {});
        return cdest;
    }

    var fin = mk_canvas(width, height);
    $('body').append(fin.canvas);

    var x0 = 0;
    var drawRules = function(first_pass) {
        if (PARAMS.dist_unit == 'none') {
            return;
        }

        var antipode = EARTH_FARTHEST;
        var units = {
            km: 1000.,
            mi: 1609.344,
            deg: antipode / 90.,
        };
        var stops = {
            km: [.03, .1, .3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000, antipode / units.km],
            mi: [150/5280., 500/5280., 1500/5280., 1, 3, 10, 30, 100, 300, 1000, 3000, antipode / units.mi],
            deg: [1/3600., 1/1200., 1/360., 1/120., 1/60., 1/20., 1/6., .5, 1, 3, 10, 30, 90],
        }
        $.each(stops[PARAMS.dist_unit], function(i, e) {
            var dist = e * units[PARAMS.dist_unit];
            var is_antipode = Math.abs(dist - antipode) < 1e-6;
            if (PARAMS.dist_unit == 'km') {
                if (is_antipode) {
                    var label = 'antipode'
                } else if (e < 1.) {
                    var label = (1000. * e) + ' m';
                } else {
                    var label = e + ' km';
                }
            } else if (PARAMS.dist_unit == 'mi') {
                if (is_antipode) {
                    var label = 'antipode'
                } else if (e < 1.) {
                    var label = Math.round(5280. * e) + ' ft';
                } else {
                    var label = e + ' mi';
                }
            } else if (PARAMS.dist_unit == 'deg') {
                if (is_antipode) {
                    var label = 'nadir'
                } else if (e < 1/60.) {
                    var label = Math.round(e * 3600.) + '\u2033 below horizon'
                } else if (e < 1.) {
                    var label = Math.round(e * 60.) + '\u2032 below horizon'
                } else {
                    var label = e + '\xb0 down'
                }
            }
            var y = disty(dist).y;

            fin.context.fillStyle = 'rgba(0, 0, 0, ' + (first_pass ? .2 : .04) + ')';
            fin.context.fillRect(0, y, width, 1);
            fin.context.fillStyle = 'rgba(0, 0, 0, .6)';
            fin.context.font = '8pt sans-serif';
            fin.context.fillText(label, x0 + 2, y - 2);
        });
    }
    drawRules(true);

    var compose = function(canv) {
        fin.context.drawImage(canv.canvas, 0, 0, width, height);
    };

    compose(downsize(c));
    var edges = downsize(c2);
    for (var i = 0; i < 2; i++) {
        compose(edges);
    }

    drawRules(false);

    fin.context.fillStyle = 'rgba(0, 0, 0, .6)';
    fin.context.textAlign = 'center';
    var bearing_tick = 15;
    var bearing_label_min = bearing_tick * (Math.floor(data.range[0] / bearing_tick) - 1);
    var bearing_label_max = bearing_tick * (Math.ceil(data.range[1] / bearing_tick) + 1);
    var is_polar = Math.abs(data.origin[0]) > 90. - 1e-6;
    for (var bearing = bearing_label_min; bearing <= bearing_label_max; bearing += bearing_tick) {
        var nbear = fixmod(bearing, 360);
        var x = (bearing - data.range[0]) * width / (data.range[1] - data.range[0]);
        var major = (bearing % 45 == 0);

        if (!is_polar) {
            if (major) {
                var label = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'][nbear / 45];
            } else {
                var label = nbear + '\xb0';
            }
        } else {
            var lon = (data.origin[0] > 0 ? fixmod(-nbear, 360) : nbear);
            if (lon == 0 || lon == 180) {
                var dir = '';
            } else if (lon < 180) {
                var dir = 'E';
            } else {
                var dir = 'W';
                lon = 360. - lon;
            }
            var label = lon + '\xb0' + dir;
        }

        if (major) {
            fin.context.font = 'bold 18pt sans-serif';
        } else {
            fin.context.font = 'bold 12pt sans-serif';
        }
        fin.context.fillText(label, x, height - 5);
    }

    fin.context.fillStyle = 'rgba(0, 0, 0, .6)';
    fin.context.textAlign = 'left';
    fin.context.font = '10pt sans-serif';
    fin.context.rotate(-.5*Math.PI);
    fin.context.fillText('data \xa9 OpenStreetMap contributors', -height + 35, width - 5);
    fin.context.rotate(.5*Math.PI);

    return fin.canvas;
}

function mk_canvas(w, h) {
    var $c = $('<canvas />');
    $c.attr('width', w);
    $c.attr('height', h);
    var c = $c[0];
    var ctx = c.getContext('2d');
    return {canvas: c, context: ctx};
}

function setupMap(data, canv) {
    var $canv = $(canv);

    COMPANION = window.open('/companion.html', 'companion', 'width=600,height=600,location=no,menubar=no,toolbar=no,status=no,personalbar=no');
    COMPANION.onbeforeunload = function() { COMPANION = null; };

    window.onbeforeunload = function() {
        if (window.COMPANION) {
            COMPANION.close();
        }
    };

    var MAX = Array($canv.width());
    var MIN = Array($canv.width());
    for (var i = 0; i < data.postings.length; i++) {
        var dist = data.postings[i][0];
        var fpx = (i + .5) / data.postings.length * $canv.width();
        var px = Math.floor(fpx);
        if (MAX[px] === undefined || (dist > MAX[px][1])) {
            MAX[px] = [fpx, dist];
        }
        if (MIN[px] === undefined || (dist < MIN[px][1])) {
            MIN[px] = [fpx, dist];
        }
    }

    MMTIMER = null;
    $canv.mousemove(function(e) {
        if (MMTIMER != null) {
            clearTimeout(MMTIMER);
        }
        MMTIMER = setTimeout(function() {
            var x = e.pageX - $canv.offset().left;
            //var y = e.pageY - $canv.offset().top;

            var xn = MIN[x][0] / $canv.width();
            var dist = MIN[x][1];
            var bearing = data.range[0] + data.lonspan * xn;

            var target = line_plotter(data.origin, bearing)(dist);
            console.log(target);
            COMPANION.postMessage({
                pos: target,
            }, '*');
        }, 100);
    });

    $(document).keydown(function(e) {
        if (e.which == 187) {
            COMPANION.postMessage({action: 'zoomin'}, '*');
        } else if (e.which == 189) {
            COMPANION.postMessage({action: 'zoomout'}, '*');
        }
    });
}


function init_companion() {
    var map = new L.Map('map', {
        attributionControl: false,
        //fadeAnimation: false,
        //zoomAnimation: false,
    });
    map.addControl(new L.Control.Scale({
        maxWidth: 125,
	    position: 'bottomright',                       
	}));
    map.setZoom(7);

    var layers = {
        'gmap': L.tileLayer('https://mts{s}.google.com/vt/lyrs=m&x={x}&y={y}&z={z}', {subdomains: '0123'}),
        'gsat': L.tileLayer('https://mts{s}.google.com/vt/lyrs=s&x={x}&y={y}&z={z}', {subdomains: '0123'}),
        'osm': L.tileLayer('http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png'),
        'mc': L.tileLayer('http://{s}.tile.cloudmade.com/1a1b06b230af4efdbb989ea99e9841af/999/256/{z}/{x}/{y}.png'),
    };
    L.control.layers(layers).addTo(map);
    map.addLayer(layers['gmap']);

    window.addEventListener("message", function(e) {
        if (e.data.action) {
            map[{zoomin: 'zoomIn', zoomout: 'zoomOut'}[e.data.action]]();
        } else {
            companion_update(map, e.data.pos);
        }
    }, false);
}

function companion_update(map, pos) {
    map.setView(pos);
}







function fixmod(a, b) {
    return ((a % b) + b) % b;
}
