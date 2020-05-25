<?php
/* NOTE Source of algorithm used:
        Karl Osen. Accurate Conversion of Earth-Fixed Earth-Centered Coordinates to Geodetic Coordinates.
        [Research Report] Norwegian University of Science and Technology. 2017. ffhal-01704943v2f
        https://hal.archives-ouvertes.fr/hal-01704943v2/document					*/


/*****************************************************************************\
|* wgs84def.h
\*****************************************************************************/
// Set to 1 to use degrees in WGS84GEO structure
define("DEGREES", 1);
// Design limits
define("MIN_LAT", -90);
define("MAX_LAT", 90);
define("MIN_LON", -360);
define("MAX_LON", 360);

/*****************************************************************************\
|* wgs84data.h
\*****************************************************************************/
define("WGS84_A", +6.37813700000000000000e+0006); 	/* a */
define("WGS84_INVF", +2.98257223563000000000e+0002);	/* 1/f */
define("WGS84_F", +3.35281066474748071998e-0003); 	/* f */
define("WGS84_INVA", +1.56785594288739799723e-0007); 	/* 1/a */
define("WGS84_INVAA", +2.45817225764733181057e-0014);	/* 1/(a^2) */
define("WGS84_B", +6.35675231424517949745e+0006); 	/* b */
define("WGS84_C", +5.21854008423385332406e+0005);	/* c */
define("WGS84_E", +8.18191908426214947083e-0002); 	/* e */
define("WGS84_EE", +6.69437999014131705734e-0003); 	/* e^2 */
define("WGS84_EED2", +3.34718999507065852867e-0003); 	/* (e^2)/2 */
define("WGS84_EEEE", +4.48147234524044602618e-0005); 	/* e^4 */
define("WGS84_EEEED4", +1.12036808631011150655e-0005);  /* (e^4)/4 */
define("WGS84_AADC", +7.79540464078689228919e+0007); 	/* (a^2)/c */
define("WGS84_BBDCC", +1.48379031586596594555e+0002); 	/* (b^2)/(c^2) */
define("WGS84_P1MEE", +9.93305620009858682943e-0001);	/* 1-(e^2) */
define("WGS84_P1MEEDAA", +2.44171631847341700642e-0014);/* (1-(e^2))/(a^2) */
define("WGS84_P1MEEDB", +1.56259921876129741211e-0007); /* (1-(e^2))/b */
define("WGS84_HMIN", +2.25010182030430273673e-0014);	/* (e^12)/4 */
define("WGS84_INVCBRT2", +7.93700525984099737380e-0001);/* 1/(2^(1/3)) */
define("WGS84_INV3", +3.33333333333333333333e-0001);	/* 1/3 */
define("WGS84_INV6", +1.66666666666666666667e-0001);	/* 1/6 */
define("WGS84_D2R", +1.74532925199432957691e-0002);	/* pi/180 */
define("WGS84_R2D", +5.72957795130823208766e+0001);	/* 180/pi */

define("FLPREC", 18);	/* 18 significant digits ~ 80bit */

# const static double definitions
static $invaa		= WGS84_INVAA;		// 1/(a^2)
static $aadc 		= WGS84_AADC;		// (a^2)/c
static $bbdcc 		= WGS84_BBDCC;		// (b^2)/(c^2)
static $l 		= WGS84_EED2;		// (e^2)/2
static $p1mee 		= WGS84_P1MEE;		// 1-(e^2)
static $p1meedaa 	= WGS84_P1MEEDAA;	// (1-(e^2))/(a^2)
static $Hmin 		= WGS84_HMIN;		// (e^12)/4
static $ll4 		= WGS84_EEEE;		// 4*(l^2) = e^4
static $ll 		= WGS84_EEEED4;		// l^2 = (e^4)/4
static $invcbrt2 	= WGS84_INVCBRT2;	// 1/(2^(1/3))
static $inv3 		= WGS84_INV3;		// 1/3
static $inv6 		= WGS84_INV6;		// 1/6

/* NOTE https://www.oc.nps.edu/oc2902w/coord/llhxyz.htm
        Verified wgs2ecef array('lat'=>49.058346944,'lon'=>8.9124975,'alt'=>292.21) => 4137.165 / 648.787 / 4795.034 (X/Y/Z-ECEF)
*/
function wgs2ecef($geo) {
    bcscale(FLPREC);
    global $aadc, $bbdcc, $p1mee;
    $lat = $geo['lat']; $lon = $geo['lon']; $alt = $geo['alt'];
    if (DEGREES) { $lat = deg2rad($geo['lat']); $lon = deg2rad($geo['lon']); }
    $coslat = cos($lat); $sinlat = sin($lat);
    $coslon = cos($lon); $sinlon = sin($lon);
    $N = bcdiv($aadc,bcsqrt(bcadd(bcpow($coslat,2),$bbdcc)));
    $d = bcmul(bcadd($N,$alt),$coslat);
    return array('x'=>bcmul($d,$coslon),'y'=>bcmul($d,$sinlon),'z'=>bcmul(bcadd(bcmul($p1mee,$N),$alt),$sinlat));
}

function ecef2wga($ecef) {
    bcscale(FLPREC);
    global $invaa, $aadc, $bbdcc, $l, $p1mee, $p1meedaa, $Hmin, $ll4, $ll, $invcbrt2, $inv3, $inv6;
    $x = $ecef['x']; $y = $ecef['y']; $z=$ecef['z'];
    $ww = bcadd(bcpow($x,2),bcpow($y,2));
    $m = $ww * $invaa;			// bcmul(sprintf("%.18f",$ww),sprintf("%.18f",$invaa)) ... doesn't work INVAA is sooo small E-14
    $n = bcpow($z,2) * $p1meedaa;	// again ...
    $mpn = bcadd($m,$n);
    $p = $inv6 * ($mpn - $ll4);		
    $G = $m * $n * $ll;
    $H = bcmul(2,bcpow($p,3))+$G;
    if ($H < $Hmin) { return -1; }
    $C = pow($H + $G + 2 * sqrt($H * $G), $inv3) * $invcbrt2;
    $i = -$ll - 0.5 * $mpn;
    $P = bcpow($p,2);
    $beta = bcsub(bcsub(bcmul($inv3,$i),$C),bcdiv($P,$C));
    $k = $ll * ($ll - $mpn);
    // Compute left part of t
    $t1 = $beta * $beta - $k;
    $t2 = sqrt($t1);
    $t3 = $t2 - 0.5 * ($beta + $i);
    $t4 = sqrt($t3);
    // Compute right part of t
    $t5 = 0.5 * ($beta - $i);
    // t5 may accidentally drop just below zero due to numeric turbulence
    // This only occurs at latitudes close to +- 45.3 degrees
    $t5 = abs($t5);
    $t6 = sqrt($t5);
    $t7 = ($m < $n) ? $t6 : -$t6;
    // Add left and right parts
    $t = $t4 + $t7;
    // Use Newton-Raphson's method to compute t correction
    $j = $l * ($m - $n);
    $g = 2 * $j;
    $tt = $t * $t;
    $ttt = $tt * $t;
    $tttt = $tt * $tt;
    $F = $tttt + 2 * $i * $tt + $g * $t + $k;
    $dFdt = 4 * $ttt + 4 * $i * $t + $g;
    $dt = -$F / $dFdt;
    // compute latitude (range -PI/2..PI/2)
    $u = $t + $dt + $l;
    $v = $t + $dt - $l;
    $w = sqrt($ww);
    $zu = $z * $u;
    $wv = $w * $v;
    $lat = atan2($zu, $wv);
    // compute altitude
    $invuv = 1 / ($u * $v);
    $dw = $w - $wv * $invuv;
    $dz = $z - $zu * $p1mee * $invuv;
    $da = sqrt($dw * $dw + $dz * $dz);
    $alt = ($u < 1) ? -$da : $da;
    // compute longitude (range -PI..PI)
    $lon = atan2($y, $x);
    if (DEGREES) { $lat = rad2deg($lat); $lon = rad2deg($lon); }
    return array('lat'=>$lat,'lon'=>$lon,'alt'=>$alt);
}

function haversine($latitudeFrom, $longitudeFrom, $latitudeTo, $longitudeTo, $earthRadius = 6371000) //GreatCircleDistance
{
  // convert from degrees to radians
  $latFrom = deg2rad($latitudeFrom);
  $lonFrom = deg2rad($longitudeFrom);
  $latTo = deg2rad($latitudeTo);
  $lonTo = deg2rad($longitudeTo);

  $latDelta = $latTo - $latFrom;
  $lonDelta = $lonTo - $lonFrom;

  $angle = 2 * asin(sqrt(pow(sin($latDelta / 2), 2) +
    cos($latFrom) * cos($latTo) * pow(sin($lonDelta / 2), 2)));
  return $angle * $earthRadius;
}

// Inputs in Radians so degtorad degrees first.
function vincenty($lat1, $lon1, $lat2, $lon2) {
        $lat1=deg2rad($lat1); $lon1=deg2rad($lon1);
        $lat2=deg2rad($lat2); $lon2=deg2rad($lon2);
        // Equitorial Radius
        $a = 6378137.0;
        // Polar Radius
        $b = 6356752.31424518;
        //Flattening of the ellipsoid
        $f = 0.00335281066;
        // Difference in longitude
        $L = $lon2 - $lon1;  
        $U1 = atan((1 - $f) * tan($lat1));  //U is 'reduced latitude'
        $U2 = atan((1 - $f) * tan($lat2));
        $sinU1 = sin($U1);
        $sinU2 = sin($U2);
        $cosU1 = cos($U1);
        $cosU2 = cos($U2);

        $lambda = $L;
        $lambdaP = 2 * pi();
        $i = 20;

        while (abs($lambda - $lambdaP) > 1e-12 && --$i > 0) {
            $sinLambda = sin($lambda);
            $cosLambda = cos($lambda);
            $sinSigma = sqrt(($cosU2 * $sinLambda) * ($cosU2 * $sinLambda) + ($cosU1 * $sinU2 - $sinU1 * $cosU2 * $cosLambda) * ($cosU1 * $sinU2 - $sinU1 * $cosU2 * $cosLambda));

            if ($sinSigma == 0)
                return 0;  //co-incident points

            $cosSigma = $sinU1 * $sinU2 + $cosU1 * $cosU2 * $cosLambda;
            $sigma = atan2($sinSigma, $cosSigma);
            $sinAlpha = $cosU1 * $cosU2 * $sinLambda / $sinSigma;
            $cosSqAlpha = 1 - $sinAlpha * $sinAlpha;
            $cos2SigmaM = $cosSigma - 2 * $sinU1 * $sinU2 / $cosSqAlpha;
            if (is_nan($cos2SigmaM))
                $cos2SigmaM = 0;  //equatorial line: cosSqAlpha=0 (6)
            $c = $f / 16 * $cosSqAlpha * (4 + $f * (4 - 3 * $cosSqAlpha));
            $lambdaP = $lambda;
            $lambda = $L + (1 - $c) * $f * $sinAlpha * ($sigma + $c * $sinSigma * ($cos2SigmaM + $c * $cosSigma * (-1 + 2 * $cos2SigmaM * $cos2SigmaM)));
        }

        if ($i == 0)
            return false;  //formula failed to converge

        $uSq = $cosSqAlpha * ($a * $a - $b * $b) / ($b * $b);
        $A = 1 + $uSq / 16384 * (4096 + $uSq * (-768 + $uSq * (320 - 175 * $uSq)));
        $B = $uSq / 1024 * (256 + $uSq * (-128 + $uSq * (74 - 47 * $uSq)));
        $deltaSigma = $B * $sinSigma * ($cos2SigmaM + $B / 4 * ($cosSigma * (-1 + 2 * $cos2SigmaM * $cos2SigmaM) - $B / 6 * $cos2SigmaM * (-3 + 4 * $sinSigma * $sinSigma) * (-3 + 4 * $cos2SigmaM * $cos2SigmaM)));
        $d = $b * $A * ($sigma - $deltaSigma);
        return $d;
    }

/* NOTE Main function for testing purposes */
$chk=array('lat'=>49.058346944,'lon'=>8.9124975,'alt'=>292.21);
$ecef=wgs2ecef($chk);
$wga=ecef2wga($ecef);
$ecef2=wgs2ecef($wga);
$ax=($chk['lat']-$wga['lat'])*3600; $ay=($chk['lon']-$wga['lon'])*3600; $az=($chk['alt']-$wga['alt'])*1000;
echo sprintf("WGS  ax=%9.6f[\"]  ay=%9.6f[\"]  az=%9.6f[mm]\r\n",$ax,$ay,$az);
$dx=($ecef['x']-$ecef2['x'])*1000; $dy=($ecef['y']-$ecef2['y'])*1000; $dz=($ecef['z']-$ecef2['z'])*1000;
echo sprintf("ECEF dx=%9.6f[mm] dy=%9.6f[mm] dz=%9.6f[mm] d=%9.6f[mm]\r\n",$dx,$dy,$dz,0.5*sqrt($dx*$dx+$dy*$dy+$dz*$dz));

#B = 49.058921   / 8.791156   / 204,700m
#K = 49.01124241 / 8.41125530 ?
$dist=haversine(49.058921,8.791156,49.01124241,8.41125530);
echo sprintf("Baseline B35OD - KARL00DEU0 = %.3f[m]\n",$dist);
$dist2=vincenty(49.058921,8.791156,49.01124241,8.41125530);
echo sprintf("Baseline B35OD - KARL00DEU0 = %.3f[m]\n",$dist2);
?>
