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
Only the `tagged_coastline` file need be kept once the process completes, but if the pipeline is re-run it will re-use any interim outputs still present in the directory. The pipeline may also be re-run in the following modes:

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
