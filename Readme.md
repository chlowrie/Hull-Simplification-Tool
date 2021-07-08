<h1>Hull Simplification Tool</h1>
This script is to create a constrained exterior hull.  It is similar to a local convex hull, where the convexity is determined by a search radius.
<br>
Flags: (*=required)
<ul>
    <li>in_shp*</li>
    <li>out_shp*</li>
    <li>--contour_line : the contour value to simplify. Default ELEV=-2.0</li>
    <li>--search_distance : Distance to search for the next shell simplification point.</li>
    <li>--length_cutoff : The percentage length of the maximum contour.  Default 0.25.  Should be set lower if the DEM contains multiple islands of different sizes.</li>
</ul>