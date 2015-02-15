function init() {
    var data = DATA_CPT;

    var canv;
    loadColors(function() {
        canv = render(data, 3000, 800);
        setupMap(data, canv);
    });
}

NUM_COLORS = 6;
STEPS = 1000;

function loadColors(after) {
    WS = new WebSocket('ws://localhost:8888/socket')
    WS.onmessage = function(e) {
        COLORS = JSON.parse(e.data);
        COLORS = _.shuffle(COLORS); // randomize the hue order
        after();
    };
    setTimeout(function() {
        WS.send(JSON.stringify({
            num_colors: NUM_COLORS,
            lum: .6,
            sat: [.03, .38],
            steps: STEPS,
        }));
    }, 500);
}

function render(data, width, height) {
    var c = mk_canvas(data.postings.length, height)
    var ctx = c.context;
    //$('body').append(c.canvas);

    var logmin = Math.log(data.min_dist);
    var logmax = Math.log(2*Math.PI * 6371009.);

    var SATMIN = .1;
    var SATMAX = .8;

    console.log(logmin, logmax);

    var admin_by_px = [];
    for (var i = 0; i < data.postings.length; i++) {
        var _admins = data.postings[i][1];
        var admin = (_admins.length > 0 ? _admins[Math.floor(Math.random() * _admins.length)] : -1);
        admin_by_px.push(admin);
    }

    // determine administrative areas that appear
    var _admins = {};
    for (var i = 0; i < admin_by_px.length; i++) {
        _admins[admin_by_px[i]] = true;
    }
    var admins = _.map(_.keys(_admins), function(k) { return +k; });

    // build adjacency graph
    var adj = {};
    var adjedge = function(a, b) {
        if (!adj[a]) {
            adj[a] = {};
        }
        adj[a][b] = true;
    };
    for (var i = 0; i < admin_by_px.length; i++) {
        var a0 = admin_by_px[i];
        var a1 = admin_by_px[(i + 1) % admin_by_px.length];
        if (a0 == a1) {
            continue;
        }
        adjedge(a0, a1);
        adjedge(a1, a0);
    }
    _.each(adj, function(v, k) {
        adj[k] = _.map(_.keys(v), function(k) { return +k; });
    });

    var colors = assign_colors(NUM_COLORS, admins, adj);
    console.log(colors);

    for (var i = 0; i < data.postings.length; i++) {
        var dist = data.postings[i][0];
        var admin = admin_by_px[i];
        if (dist >= 0) {
            var logdist = Math.log(dist);
            var k = Math.max((logdist - logmin) / (logmax - logmin), 0);
            var y = height * k;

            var C = COLORS[colors[admin]][Math.floor(k * (STEPS - 1))];
            ctx.fillStyle = C;
            ctx.fillRect(i, y, 1, height - y);
        }
    }

    var c2 = mk_canvas(data.postings.length, height)
    var ctx2 = c2.context;
    //$('body').append(c2.canvas);

    for (var i = 0; i < data.postings.length; i++) {
        var dist0 = data.postings[i][0];        
        var dist1 = data.postings[(i + 1) % data.postings.length][0];
        var diff = Math.abs(dist0 - dist1);
        var closer = Math.min(dist0, dist1);
        var farther = Math.max(dist0, dist1);
        var quantum = closer * Math.PI / 180. * data.res;
        if (diff >= 10 * quantum) {
            var logdiff = Math.abs(Math.log(dist0) - Math.log(dist1));
            var k = (Math.log(farther) - logmin) / (logmax - logmin);
            var y = height * k;
            
            var pixelw = data.postings.length / width;
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
    var compose = function(canv) {
        fin.context.drawImage(canv.canvas, 0, 0, width, height);
    };

    compose(downsize(c));
    var edges = downsize(c2);
    for (var i = 0; i < 2; i++) {
        compose(edges);
    }

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

    COMPANION = window.open('companion.html', 'companion', 'width=600,height=600,location=no,menubar=no,toolbar=no,status=no,personalbar=no');
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
            var bearing = data.range[0] * (1-xn) + data.range[1] * xn;

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

function assign_colors(num_colors, nodes, adj) {
    var matrix = [];
    for (var i = 0; i < nodes.length; i++) {
        var row = [];
        for (var j = 0; j < nodes.length; j++) {
            row.push(0);
        }
        matrix.push(row);
    }
    _.each(adj, function(v, k) {
        k = +k;
        _.each(v, function(e, i) {
            matrix[nodes.indexOf(k)][nodes.indexOf(e)] = 1;
        });
    });

    var coloring = colorizing(matrix);
    var colorcount = _.max(coloring);
    for (var i = 0; i < coloring.length; i++) {
        coloring[i] -= 1;
    }

    // we have randomized the hue ordering, so don't need to
    // worry about any biasing effect
    if (colorcount > num_colors) {
        // too many colors
        for (var i = 0; i < coloring.length; i++) {
            coloring[i] = coloring[i] % num_colors;
        }
    } else if (colorcount < num_colors) {
        // too few colors
        var remap = {};
        for (var i = 0; i < num_colors; i++) {
            var c = i % colorcount;
            if (!remap[c]) {
                remap[c] = [];
            }
            remap[c].push(i);
        }
        for (var i = 0; i < coloring.length; i++) {
            var opts = remap[coloring[i]];
            coloring[i] = opts[Math.floor(Math.random() * opts.length)];
        }
    }

    var colormap = {};
    for (var i = 0; i < nodes.length; i++) {
        colormap[nodes[i]] = coloring[i];
    }
    return colormap;
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






function colorizing(a){  
  // this function finds the unprocessed vertex of which degree is maximum
  var MaxDegreeVertex = function(){
    var max = -1;
    var max_i;
    for (var i =0; i<a.length; i++){
      if ((color[i] === 0) && (degree[i]>max)){
        max = degree[i];
        max_i = i;
      }
    }
    return max_i;
  };

  // find the vertex in NN of which degree is maximum
  var MaxDegreeInNN = function(){
    var max = -1;
    var max_i, i;
    for (var k=0; k<NN.length; k++){
      i = NN[k];
      if ((color[i] === 0) && (degree[i]>max)){
        max = degree[i];
        max_i = i;
      }
    }
    return max_i;
  };

  // this function updates NN array
  var UpdateNN = function(ColorNumber){
    NN = [];
    for (var i=0; i<a.length; i++){
      if (color[i] === 0){
        NN.push(i);
      }
    }
    for (var i=0; i<a.length; i++){
      if (color[i] === ColorNumber){
        for (var j=0; j<NN.length; j++){
          while (a[i][NN[j]] === 1){
            NN.splice(j,1)
          }
        }
      }
    }
  };

  // this function will find suitable y from NN
  var FindSuitableY = function(ColorNumber,VerticesInCommon){
    var temp,tmp_y,y=0;
    var scanned = [];
    VerticesInCommon = 0;
    for (var i=0; i<NN.length; i++) {
      tmp_y = NN[i];
      temp = 0;
      for (var f=0; f<a.length; f++){
        scanned[f] = 0;
      }
      for (var x=0; x<a.length; x++){
        if (color[x] === ColorNumber){
          for (var k=0; k<a.length; k++){
            if (color[k] === 0 && scanned[k] === 0){
              if (a[x][k] === 1 && a[tmp_y][k] === 1){
                temp ++;
                scanned[k] = 1;
              }
            }
          }
        }
      }
      if (temp > VerticesInCommon){
        VerticesInCommon = temp;
        y = tmp_y;
      }
    }
    return [y,VerticesInCommon];
  };

  var color = [];
  var degree = [];
  var NN = [];

  for (var i=0; i<a.length; i++){
    color[i] = 0;
    degree[i] = 0;
    for (var j = 0; j<a.length; j++){
      if (a[i][j] === 1){
        degree[i] ++;
      }
    }
  }
  
  var x,y;
  var result;
  var ColorNumber = 0;
  var VerticesInCommon = 0;
  var unprocessed = a.length;
  while (unprocessed>0){
    x = MaxDegreeVertex();
    color[x] = ++ColorNumber;
    unprocessed --;
    UpdateNN(ColorNumber);
    while (NN.length>0){
      result = FindSuitableY(ColorNumber, VerticesInCommon);
      y = result[0];
      VerticesInCommon = result[1];
      if (VerticesInCommon === 0){
        y = MaxDegreeInNN();
      }
      color[y] = ColorNumber;
      unprocessed --;
      UpdateNN(ColorNumber);
    }
  }
  return color
};