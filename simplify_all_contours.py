from simplify_contours import main

contours_list = (
    ('data/AmericanSamoa/ASamoa_Tutuila_10m_negdepths_UTM2S84_contours.shp', 'outputs/Tutuila_contours_simplified.shp'),
    ('data/ASamoa_Manua_10m_negdepths_UTM2S84_contours.shp', 'outputs/Manua_contours_simplified.shp'),
    ('data/Guam/Guam_10m_negdepths_UTM55N84_contours.shp', 'outputs/Guam_contours_simplified.shp'),
    ('data/HawaiianIslands/BigIsland/big_island_contours.shp', 'outputs/BigIsland_contours_simplified.shp'),
    ('data/HawaiianIslands/Kahoolawe/kahului_1s_20050909_projected_gis_contours.shp', 'outputs/Kahului_contours_simplified.shp'),
    ('data/HawaiianIslands/LanaiMolokai/LanaiMolokai_10m_negdepths_UTM4N84_contours.shp', 'outputs/LanaiMolokai_contours_simplified.shp'),
    ('data/HawaiianIslands/Maui/MauiCandidate0524_contours.shp', 'outputs/MauiCandidate0524_contours_simplified.shp'),
    ('data/HawaiianIslands/Niihau/niihau_10m_negdepths_UTM4N84_contours.shp', 'outputs/Niihau_contours_simplified.shp'),
    ('data/HawaiianIslands/Oahu/Oahu_10m_nonoaa_contours.shp', 'outputs/Oahu_contours_simplified.shp'),
    ('data/SaipanTinian/SaipanTinian_10m_negdepths_UTM55N84_contours.shp', 'outputs/SaipanTinian_contours_simplified.shp')
)

for i in contours_list:
    main(i[0], i[1], -2.0, 2000, 0.125)
