# CoConverse
Accurate conversion of Earth-Fixed Earth-Centered Coordinates to Geodetic Coordinates.

"A closed form algorithm for the exact transformation of Earth-Centered Earth-Fixed (ECEF)
coordinates to geodetic coordinates is presented that is computationally fast, safe and accurate.
Boosting the computational robustness of Jijie Zhu’s algorithm, it reduces the worst case
transformation error by up to 500 million times."

Karl Osen. Accurate Conversion of Earth-Fixed Earth-Centered Coordinates to Geodetic Coordinates.
[Research Report] Norwegian University of Science and Technology. 2017. 
https://hal.archives-ouvertes.fr/hal-01704943v2/document

Porting the algorithms described in Mr. Osen Paper from C language to php.
+ Adding haversine (found on https://stackoverflow.com/questions/365826/calculate-distance-between-2-gps-coordinates)
+ and Vincenty's algorithm (https://gist.github.com/lkacenja/41b5bd876d61b806433f)

... for calculation of baseline distance between to GNSS marks.
