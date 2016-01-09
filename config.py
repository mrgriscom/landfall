
# Path where final land-casts are stored
OUTPUT_PATH = 'output'

# Areas to suppress warnings that the land is not claimed by any entity
expected_unclaimed = {
    'SPR': {
        'name': 'Spratly Islands',
        'bounds': 'POLYGON ((117.18018 9.34067,115.82886 7.89603,113.58765 6.58822,111.92871 6.91552,111.22559 8.58102,112.57690 10.94119,114.71924 12.18970,116.89453 11.96410,117.36694 10.84410,117.18018 9.34067))',
    },
    'DOU': {
        'name': 'Doumeira Islands',
        'bounds': 'POLYGON ((43.13361 12.72098,43.14417 12.70758,43.16279 12.71562,43.15207 12.72550,43.13361 12.72098))',
    },
}

# Areas/countries to suppress warnings that the country's land is not fully
# covered by subdivisions.
incompletely_subdivided = {
    'CN': 'all',
    'NO': 'all',
    'FI': 'all',
    'CA': 'POLYGON ((-76.72852 42.87596,-79.84863 47.54687,-88.46191 50.17690,-94.39453 46.61926,-89.03320 40.01079,-81.43066 40.24599,-76.72852 42.87596))', # great lakes
    'US': '''MULTIPOLYGON (
            ((-76.72852 42.87596,-79.84863 47.54687,-88.46191 50.17690,-94.39453 46.61926,-89.03320 40.01079,-81.43066 40.24599,-76.72852 42.87596)),
            ((-75.44 39.79,-75.44 39.89,-75.23 39.89,-75.23 39.79,-75.44 39.79))
          )''', # great lakes + patch of delaware river in PA
}

# Areas to suppress warnings that land is claimed by multiple entities, i.e.,
# a known conflict.
disputed_areas = {
    'CR': {
        'name': 'Crimea',
        'parties': ['RU', 'UA'],
        'bounds': 'POLYGON ((33.59619 46.31658,31.81641 45.62940,33.57422 44.03232,36.72729 44.73893,36.83716 45.72152,33.59619 46.31658))',
    },
    'HAN': {
        'name': 'Hans Island',
        'parties': ['CA-NU', 'GL'],
        'bounds': 'POLYGON ((-66.46042 80.84053,-66.34781 80.82281,-66.45561 80.81152,-66.56342 80.82937,-66.46042 80.84053))',
    },
}

# Further confliced areas where an entire entity is the conflict area.
defacto = {
    'EH': ['MA'],
    'XAB': ['GE'],
    'XNC': ['CY'],
}

