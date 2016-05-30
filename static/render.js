EARTH_MEAN_RAD = 6371009.;
EARTH_ANTIPODE = Math.PI * EARTH_MEAN_RAD;
EARTH_CIRCUMF = 2 * Math.PI * EARTH_MEAN_RAD;

UNITS = {
    km: 1000.,
    mi: 1609.344,
    nmi: 1852.,
    deg: EARTH_ANTIPODE / 90.,
};

function init() {
    // DATA
    // PARAMS

    var canvas = render();
    $('#map').click(function() {
        setupMap(canvas);
    });
    $('#png').click(function() {
        var timestamp = Math.floor(new Date().getTime() / 1000.);
        var filename = DATA.tag + '_' + timestamp;
        save_canvas(canvas, filename + '.png');
        write_params(filename + '.params');
    });
}

function render() {
    var width = PARAMS.dim[0];
    var height = PARAMS.dim[1];

    var raw_width = Math.round(width * PARAMS.res / DATA.res);
    var raw_offset = Math.round(fixmod(PARAMS.bear0 - DATA.range[0], 360.) / DATA.res);

    var c = mk_canvas(raw_width, height);
    var ctx = c.context;
    //$('body').append(c.canvas);

    var logmin = Math.log(PARAMS.min_dist || DATA.min_dist);
    var logmax = Math.log(PARAMS.max_dist || EARTH_CIRCUMF);
    var colorlogmin = (PARAMS.colordists[0] == 0 ? logmin : Math.log(PARAMS.colordists[0]));
    var colorlogmax = (PARAMS.colordists[1] == 0 ? logmax : Math.log(PARAMS.colordists[1]));

    var disty = function(dist, no_clip) {
        var logdist = Math.log(dist);
        var k = Math.max((logdist - logmin) / (logmax - logmin), no_clip ? -1e9 : 0);
        return height * k;
    }
    var distramp = function(dist) {
        var logdist = Math.log(dist);
        var k = (logdist - colorlogmin) / (colorlogmax - colorlogmin);
        return Math.min(Math.max(k, 0.), 1.);
    }

    console.log(DATA.colors);

    // land
    for (var i = 0; i < raw_width; i++) {
        var i_posting = fixmod(i + raw_offset, DATA.postings.length);
        var dist = DATA.postings[i_posting][0];
        var admin = DATA.admin_postings[i_posting];
        if (dist >= 0) {
            var y = disty(dist);
            var ramp = PARAMS.colors[DATA.colors[admin]];
            var C = ramp[Math.floor(distramp(dist) * (ramp.length - 1))];
            ctx.fillStyle = C;
            ctx.fillRect(i, y, 1, height - y);
        }
    }

    var c2 = mk_canvas(raw_width, height)
    var ctx2 = c2.context;
    //$('body').append(c2.canvas);

    // creases
    for (var i = 0; i < raw_width - 1; i++) {
        var i_posting = fixmod(i + raw_offset, DATA.postings.length);
        var dist0 = DATA.postings[i_posting][0];
        var dist1 = DATA.postings[(i_posting + 1) % DATA.postings.length][0];
        var closer = Math.min(dist0, dist1);
        var farther = Math.max(dist0, dist1);
        var quantum = closer * Math.PI / 180. * DATA.res;

        var edge_weight_px = function(farther) {
            if (farther < closer) {
                return 0;
            }
            var diff = Math.abs(closer - farther);
            if (diff < 10 * quantum) {
                return 0;
            }
            var logdiff = Math.abs(Math.log(closer) - Math.log(farther));
            var pixelw = raw_width / width;
            var weight = 6 * logdiff / (logmax - logmin);
            var weight2 = Math.min(diff / 5e5, .3);
            weight = Math.max(weight, weight2);
            return pixelw * Math.min(weight, 1.) / 3;
        };

        var _w = edge_weight_px(farther);
        if (_w == 0) {
            continue;
        }
        for (var t = 1; t < _w; t++) {
            var k = (dist0 > dist1 ? i - t : i+1 + t);
            if (k < 0 || k >= raw_width) {
                break;
            }
            var k_posting = fixmod(k + raw_offset, DATA.postings.length);
            var dist = DATA.postings[k_posting][0];
            var _w2 = edge_weight_px(dist);
            if (_w2 < t) {
                _w = t;
                break;
            }
            _w = Math.min(_w, _w2);
        }

        var x0 = i+1 - (dist0 > dist1 ? _w : 0);
        ctx2.fillStyle = '#000';
        ctx2.fillRect(x0, disty(farther), _w, height);
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

    var ylabel_x0 = 0;
    var drawRules = function(first_pass) {
        if (PARAMS.dist_unit == null) {
            return;
        }

        var FT_NM = .3048 / 1852.;
        var stops = {
            km: [.01, .03, .1, .3, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000, 'antip'],
            mi: [50/5280., 150/5280., 500/5280., 1500/5280., 1, 3, 10, 30, 100, 300, 1000, 3000, 'antip'],
            nmi: [60*FT_NM, 180*FT_NM, 600*FT_NM, 1800*FT_NM, 1, 3, 10, 30, 100, 300, 1000, 3000, 'antip'],
            deg: [1/3600., 1/1200., 1/360., 1/120., 1/60., 1/20., 1/6., .5, 1, 3, 10, 30, 'antip'],
        }
        $.each(PARAMS.yticks || stops[PARAMS.dist_unit], function(i, e) {
            var dist = (e == 'antip' ? EARTH_ANTIPODE : e * UNITS[PARAMS.dist_unit]);
            var is_antipode = Math.abs(dist - EARTH_ANTIPODE) < 1e-6;
            if (PARAMS.dist_unit == 'km') {
                if (is_antipode) {
                    var antip_units = dist / UNITS[PARAMS.dist_unit];
                    var label = 'antipode (' + Math.round(antip_units) + ' km)'
                } else if (e < 1.) {
                    var label = (1000. * e) + ' m';
                } else {
                    var label = e + ' km';
                }
            } else if (PARAMS.dist_unit == 'mi') {
                if (is_antipode) {
                    var antip_units = dist / UNITS[PARAMS.dist_unit];
                    var label = 'antipode (' + Math.round(antip_units) + ' mi)'
                } else if (e < 1.) {
                    var label = Math.round(5280. * e) + ' ft';
                } else {
                    var label = e + ' mi';
                }
            } else if (PARAMS.dist_unit == 'nmi') {
                if (is_antipode) {
                    var antip_units = dist / UNITS[PARAMS.dist_unit];
                    var label = 'antipode (' + Math.round(antip_units) + ' nmi)'
                } else if (e < 1.) {
                    var label = Math.round(e / FT_NM / 3) + ' yd';
                } else {
                    var label = e + ' nmi';
                }
            } else if (PARAMS.dist_unit == 'deg') {
                var _bh = ' below horizon';
                if (is_antipode) {
                    var label = 'nadir'
                } else if (e < 1/60.) {
                    var label = Math.round(e * 3600.) + '\u2033' + _bh;
                } else if (e < 1.) {
                    var label = Math.round(e * 60.) + '\u2032' + _bh;
                } else {
                    var label = e + '\xb0' + _bh;
                }
            }
            var y = disty(dist, true);
            if (y < -1e-6) {
                return;
            }

            var ylabel_fontsize = 8;
            var ylabel_padding = 2;
            var label_below = (y < ylabel_fontsize + 2 * ylabel_padding);

            fin.context.fillStyle = 'rgba(0, 0, 0, ' + (first_pass ? .2 : .04) + ')';
            fin.context.fillRect(0, y, width, 1);
            fin.context.fillStyle = 'rgba(0, 0, 0, .6)';
            fin.context.font = ylabel_fontsize + 'pt sans-serif';
            fin.context.textBaseline = (label_below ? 'top' : 'alphabetic');
            for (var i = 0; i < PARAMS.ylabelreps; i++) {
                var x = width * i / PARAMS.ylabelreps + ylabel_x0;
                fin.context.fillText(label, x + ylabel_padding, y - (label_below ? -1 : 1) * ylabel_padding);
            }
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

    var bearing_preferred_tick_spacing = 150; // px; TODO make param
    var bearing_preferred_ticks = [90, 30, 15, 5, 1];
    var bearing_tick = _.min(bearing_preferred_ticks, function(e) {
        return Math.abs(Math.log(e / PARAMS.res / bearing_preferred_tick_spacing));
    });
    var bearing_label_min = bearing_tick * (Math.floor(PARAMS.bear0 / bearing_tick) - 1);
    var bearing_label_max = bearing_tick * (Math.ceil((PARAMS.bear0 + PARAMS.bearspan) / bearing_tick) + 1);
    var is_polar = Math.abs(DATA.origin[0]) > 90. - 1e-6;
    for (var bearing = bearing_label_min; bearing <= bearing_label_max; bearing += bearing_tick) {
        var nbear = fixmod(bearing, 360);
        var x = (bearing - PARAMS.bear0) * width / PARAMS.bearspan;
        var major = (bearing % 45 == 0);

        if (!is_polar) {
            if (major) {
                var label = ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'][nbear / 45];
            } else {
                var label = nbear + '\xb0';
            }
        } else {
            var lon = (DATA.origin[0] > 0 ? fixmod(-nbear, 360) : nbear);
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
        fin.context.fillStyle = 'rgba(0, 0, 0, .6)';
        fin.context.textAlign = 'center';
        fin.context.fillText(label, x, height - 5);
    }

    fin.context.fillStyle = 'rgba(0, 0, 0, .6)';
    fin.context.font = '10pt sans-serif';
    fin.context.rotate(-.5*Math.PI);

    fin.context.textAlign = 'left';
    fin.context.fillText(PARAMS.primary_annotation, -height + 35, width - 5);

    fin.context.textAlign = 'right';
    fin.context.fillText(PARAMS.secondary_annotation, -5, width - 5);

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

function setupMap(canv) {
    if (window.COMPANION) {
        return;
    }

    var $canv = $(canv);

    COMPANION = window.open('/companion.html', 'companion', 'width=600,height=600,location=no,menubar=no,toolbar=no,status=no,personalbar=no');
    COMPANION.onbeforeunload = function() { COMPANION = null; };

    window.onbeforeunload = function() {
        if (window.COMPANION) {
            COMPANION.close();
        }
    };

    var width = $canv.width();
    var raw_width = Math.round(width * PARAMS.res / DATA.res);
    var raw_offset = Math.round(fixmod(PARAMS.bear0 - DATA.range[0], 360.) / DATA.res);

    var MAX = Array($canv.width());
    var MIN = Array($canv.width());
    for (var i = 0; i < raw_width; i++) {
        var i_posting = fixmod(i + raw_offset, DATA.postings.length);
        var dist = DATA.postings[i_posting][0];
        var admin = DATA.admin_postings[i_posting];
        var fpx = (i + .5) / raw_width * width;
        var px = Math.floor(fpx);
        if (MAX[px] === undefined || (dist > MAX[px][1])) {
            MAX[px] = [fpx, dist, admin];
        }
        if (MIN[px] === undefined || (dist < MIN[px][1])) {
            MIN[px] = [fpx, dist, admin];
        }
    }

    $canv.mousemove(function(e) {
        if (window.FREEZE || window.COMPANION == null) {
            return;
        }

        var x = e.pageX - $canv.offset().left;
        //var y = e.pageY - $canv.offset().top;

        var bearing = PARAMS.bear0 + PARAMS.res * MIN[x][0];
        var dist = MIN[x][1];
        var admin = MIN[x][2];
        var target = line_plotter(DATA.origin, bearing)(dist);
        var unit = {label: PARAMS.dist_unit || 'km'};
        unit.size = UNITS[unit.label];
        COMPANION.postMessage({
            origin: DATA.origin,
            pos: target,
            bearing: bearing,
            dist: dist,
            unit: unit,
            admin_code: admin,
            admin_name: DATA.admin_info[admin].name,
            color: DATA.colors[admin] + 1,
        }, '*');
    });

    register_keyboard(COMPANION);
}

function register_keyboard(w) {
    $(document).keydown(function(e) {
        if (e.which == 187) {
            w.postMessage({action: 'zoomin'}, '*');
        } else if (e.which == 189) {
            w.postMessage({action: 'zoomout'}, '*');
        } else if (e.which == 32) {
            FREEZE = true;
        }
    });
    $(document).keyup(function(e) {
        if (e.which == 32) {
            FREEZE = false;
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
    };
    L.control.layers(layers).addTo(map);
    map.addLayer(layers['gmap']);

    MMTIMER = null;
    LAST_MAP_MOVE = 0;
    window.addEventListener("message", function(e) {
        if (e.data.action) {
            map[{zoomin: 'zoomIn', zoomout: 'zoomOut'}[e.data.action]]();
        } else {
            if (window.FREEZE) {
                return;
            }

            $('#pos').text(fmt_ll(e.data.pos[0], 'NS', 2, 3) + ' ' + fmt_ll(e.data.pos[1], 'EW', 3, 3));
            $('#relpos').text(fmt_ang(e.data.bearing, 3, 2) + ' \xd7 ' + (e.data.dist / e.data.unit.size).toFixed(2) + e.data.unit.label);
            $('#admin').text(e.data.admin_code + ' \u2013 ' + e.data.admin_name);
            $('#color').text(e.data.color);

            var update = function() {
                map.setView(e.data.pos);
            };
            if (MMTIMER != null) {
                clearTimeout(MMTIMER);
            }
            var now = performance.now() / 1000;
            if (now - LAST_MAP_MOVE > .4) {
                update();
                LAST_MAP_MOVE = now;
            } else {
                MMTIMER = setTimeout(update, 100);
            }

            setMapOverlay(map, 'cast', new L.geoJson(mk_geojson(e.data, 256 * Math.pow(2, map.getZoom())), {
                style: function(feature) {
                    return {
                        color: 'red',
                        opacity: .8,
                        weight: 2,
                    }
                },
                onEachFeature: function (feature, layer) {
                    layer.options.smoothFactor = 0;
                }
            }));
        }
    }, false);

    register_keyboard(window);
}








function fixmod(a, b) {
    return ((a % b) + b) % b;
}

function npad(n, pad) {
    var s = '' + n;
    while (s.length < pad) {
        s = '0' + s;
    }
    return s;
}

function fmt_ang(k, pad, prec) {
    return npad(fixmod(k, 360.).toFixed(prec), prec + 1 + pad) + '\xb0';
}

function fmt_ll(k, dir, pad, prec) {
    return dir[k >= 0 ? 0 : 1] + fmt_ang(Math.abs(k), pad, prec);
};




function setMapOverlay(map, name, geom) {
    if (map[name]) {
        map.removeLayer(map[name]);
    }
    map.addLayer(geom);
    map[name] = geom;
}

// nicely split up geometry that crosses the IDL
function process_geo(points) {
    var segments = [];
    var segment = [];
    for (var i = 0; i < points.length - 1; i++) {
        var a = points[i];
        var b = points[i + 1];
        segment.push(a);
        if (Math.abs(a[0] - b[0]) > 180.) {
            segment.push([unwraparound(a[0], b[0], 360), b[1]]);
            segments.push(segment);
            segment = [];
            segment.push([unwraparound(b[0], a[0], 360), a[1]]);
        }
    }
    segment.push(b);
    segments.push(segment);
    return segments;
}

function mk_geojson(data, scale_px) {
    var geojson = {
        type: 'FeatureCollection',
        features: []
    };
    var addMulti = function(props, points) {
        // need to split into separate linestring features because
        // leaflet geojson has issues with multi*
        _.each(process_geo(points), function(e) {
            geojson.features.push({
                type: 'Feature',
                properties: props,
                geometry: {
                    type: 'LineString', //'MultiLineString',
                    coordinates: e,
                }
            });
        });
    }
    addMulti({name: 'line'}, lineplot(data.origin, data.bearing, data.dist, scale_px, 'end')); 
    return geojson;
}





function write_params(filename) {
    var data = 'origin: ' + DATA.origin[0] + ',' + DATA.origin[1] + '\nmindist: ' + DATA.min_dist + '\nlonleft: ' + PARAMS.bear0 + '\nlonspan: ' + PARAMS.bearspan + '\nurl: ' + window.location.href + '\nrandom_seed: ' + PARAMS.random_seed + '\n';
    var blob = new Blob([data], {type : 'text/plain'});
    saveAs(blob, filename);
}

function save_canvas(canvas, filename) {
    var dataurl = canvas.toDataURL();
    var rawb64 = dataurl.substring(dataurl.indexOf(',')+1);
    var raw = atob(rawb64);
    var buf = new Uint8Array(raw.length);
    for (var i = 0; i < raw.length; i++) {
        buf[i] = raw.charCodeAt(i);
    }

    var blob = new Blob([buf], {type : 'image/png'});
    saveAs(blob, filename);
}
