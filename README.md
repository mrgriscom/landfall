Preparation
===========

You must first build the coastline index.
Run `process_data.py` to kick off the process.
This will:

- pull the latest coastline data from openstreetmap servers.
- amend the coastline data with some hard-coded corrections.
- build a spatial index of the coastline.
- tag coastline with administrative areas.
- perform validation checks.

**The whole process may take up to several hours and consume several gigabytes of disk.**

Outputs are stored in `data/tmp/`.
Only the `tagged_coastline` file need be kept once the process completes, but if the pipeline is re-run it will re-use any interim outputs still present in the directory.
The pipeline may also be re-run in the following modes:

- `--refresh_coast` -- pulls new coast data from openstreetmap
- `--update_coast_corrections` -- changes to coast corrections (in `data/no_land/`)
- `--update_admin_areas` -- changes to admin areas (in `data/admin_areas` or `data/admin_index.csv`)

These flags may be combined.

Validation
----------

Several sanity checks are performed in the final stage of building the coastline index.
You may see the following errors:

- **"unexpected unclaimed land"** -- this means there is land that is not claimed by any entity.
  This could mean several things:

  - Most rarely, the land is actually unclaimed -- real _terra nullius_ -- in which case it should be noted in `config.expected_unclaimed` to suppress the warning.

  - The land is actually part of a country and the country border should be edited to include it.
    This scenario can take several forms:

    - The maritime border is shoddily drawn and misses an outlying island.

    - The land is disputed and omitted from either country's territory in an effort to be 'diplomatic'.
      In this case both country's borders should be drawn to include the land and then handle the conflict under "disputed areas" below.

    - Adjacent country borders don't exactly line up and leave a small sliver of unclaimed land between them.

  - Most probably, the land is a phantom island that doesn't actually exist.

- **"parent without a subdivision"** -- this means there is land that is part of a parent entity that has subdivisions, but the land is not part of any subdivision.

  - Some countries are not fully decomposed into subdivisions.
    For such countries, specify `'all'` in `config.incompletely_subdivided` to suppress the warning.

  - Almost all the scenarios for "unclaimed land" still apply:

    - Land is 'unclaimed' by any subdivision within the country, in which case the area should be specified in `config.incompletely_subdivided` to suppress the warning.

    - The subdivision border is drawn such that it misses an outlying island.

    - Borders between adjacent subdivisions may leave slivers.

    - The land isn't actually real.

- **"subdivision not covered by parent"** -- this means the land is part of a country subdivision but is *not* included in the territory of the parent country.
  This is always an error-- the parent country boundary must be a superset of its subdivisions.
  Either the parent or subdivision boundary must be modified to make this invariant hold.

- **"unexpected disputed area"** -- this means the land is claimed by more than one entity.

  - If this is a genuine conflict, mark in `config.disputed_areas` or `config.defacto` to suppress the warning and ensure proper rendering (the disputed area can be assigned to either entity during rendering).

  - Adjacent borders may slightly overlap to create slivers of disputed area.
    The borders should be made to match up.

  - This also occurs when a border exactly follows coastline (which is a common occurrence in real life).
    This is one case where the checker could be smarter to figure out which entity actually claims the 'land' side of the coast.
    However, it is left as a warning because boundaries exactly following the coastline is a bad idea.
    The coastline may be updated dynamically while the boundaries database is static, so if at some future point the coastline is made more precise, the boundary will no longer coincide and you'll get artifacts.
    Set the boundary buffered some distance from the coast to fix this warning.

Usage
=====

Launch the server: `python server.py`

Generate a panorama
-------------------

Visit `http://localhost:8000/landfall?origin=<lat>,<lon>&size=<size>`, where you fill in `<lat>`, `<lon>`, and `<size>`.
`<size>` is the number of raw samples to take and should be 4-8x the desired width of your final panorama.

Other supported parameters:

- `range=<start>,<end>` -- only compute between the compass bearings `<start>` and `<end>` (the default is full 360&deg;).
- `mindist=<n>` -- ignore any land within `<n>` meters of the origin point.

This will load the coastline index and compute the closest landfall in all directions.
Loading the index the first time takes several minutes.
Afterward the index stays loaded in memory and subsequent runs are much faster.

The generated data is saved in `output/` and browseable from the server main page `http://localhost:8000/`.
Once complete, the data will also be rendered into a panorama with default parameters.

Render a panorama
-----------------

There are many settings to change the appearance of the rendered panorama.
Add them as url params to the `http://localhost:8000/render/`... url to use them.

### Sizing and labeling

- `width=<pixels>` -- width of rendered image
- `height=<pixels>` -- height of rendered image
- `distunit=<unit>` -- unit to indicate distance: `km`, `mi`, `nmi` (nautical miles), `deg` (declination -- *not* degrees of arc), or `none` (no y-axis labels)
- `yticks=<n>,<n>,<n>` -- tick marks for the distance (y) scale.
  Comma-separated list of numbers in the current axis unit.
  The string `antip` may also be used in lieu of a number to designate the distance to antipode.
- `min_downscale=<n>` -- the generated panorama data is downscaled for rendering to give a smoother appearance.
  However, if too much downscaling is used, tiny islands may become very faint and even lost.
  If more than `2^(n+1)` factor of downscaling is required, the source data is resampled first.
  Basically, the minimum opacity of the tiniest island will be `100% / 2^(n+1)` (12.5% for the default value of 2).

### Cropping and layout

- `mindist=<meters>` -- near limit of rendered distance
- `maxdist=<meters>` -- far limit of rendered distance
- `left=<a>` -- set the compass bearing of the left edge of the image
- `right=<a>` -- set the compass bearing of the right edge of the image
- `trimleft=<n>` -- trim `<n>` degrees of the left edge of the image
- `trimright=<n>` -- trim `<n>` degrees of the right edge of the image
- `bearcenter=<a>` -- center the image on this compass bearing
- `bearmargin=<n>` -- add `<n>` degrees of wraparound to the edges of a 360&deg; panorama (negative values allowed)

### Admin areas

- `forcecolor=<admin>:<n>,<admin>:<n>` -- force admin area with code `<admin>` to have a given color number `<n>` (0-indexed)
- `diffcolor=<admin>:<admin>,<admin>:<admin>:<admin>` -- force all groups of admin areas separated by colons to have different colors
- `nosubdiv=<admin>,<admin>` -- do not render subdivisions of the specified admin areas
- `resolve=<admin>:<admin>,<admin>:<admin>` -- for each pair of admin areas `A:B`, render `A` as if it were part of `B`

### Color palette

- `numcolors=<n>` -- number of distinct colors to use for admin areas.
- `hues=<hue>,<hue>,<hue>` -- actual hues (0--360) to use for coloring admin areas.
  Overrides `numcolors`.
- `lumclose=<f>` -- luminance (0--1) of the nearest rendered distance
- `lumfar=<f>` -- luminance (0--1) of the farthest rendered distance
- `satclose=<f>` -- saturation (0 -- ~1) of the nearest rendered distance
- `satfar=<f>` -- saturation (0 -- ~1) of the farthest rendered distance
- `colorneardist=<meters>` -- use a different reference 'closest distance' for color rendering
- `colorfardist=<meters>` -- use a different reference 'farthest distance' for color rendering
- `randseed=<n>` -- random seed.
  Specify to get deterministic results across renders.
  The seed used in a given render is saved in the export params.

Annotate a panorama
-------------------

Once you're satisfied with the appearance of the panorama, click 'export png' to save the image to disk and also write out the parameters used to render it.

It is likely you'll want to add labels in an external image editor.
'launch map' is useful to figure out which landforms are which.

Alternate output formats
------------------------

A landfall panorama may also be exported as KML.

Use the `photosphere.py` script to turn a panorama image (even one that has been edited with annotations) into a photosphere.
The `.params` file is needed.
Do not crop the rendered image.