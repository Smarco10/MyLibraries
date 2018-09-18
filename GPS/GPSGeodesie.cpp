//
//  GPSGeodesie.cpp
//  test
//
//  Created by Marc-Antoine MARTIN on 01/04/2016.
//  Copyright © 2016 Marc-Antoine MARTIN. All rights reserved.
//

#include "GPSGeodesie.h"
#include "../Matrice.h"

using namespace std;

const double GRS_a = 6378137;
const double GRS_f = 1/298.257222101;
const double GRS_b = GRS_a*(1-GRS_f);
const double GRS_e = sqrt((pow(GRS_a,2) - pow(GRS_b,2)) / pow(GRS_a,2));
const double GRS_e_2 = pow(GRS_e, 2);

void GPSGeodesie::ECEF_2_Geographique(double X, double Y, double Z, double& longitude, double& latitude, double& he)
{
	double p = sqrt(X*X + Y*Y); // utilitaire de calcul
	longitude = 2 * atan2(Y, (X + p)); // longitude directement
	
	double phiOld = atan2(Z, p); // latitude a iterer
	double w = sqrt(1 - GRS_e_2 * pow(sin(phiOld), 2));
	double N = GRS_a / w;
	he = p * cos(phiOld) + Z * sin(phiOld) - GRS_a * w;
	double phiNew = atan2((Z / p), ((1 - N * GRS_e_2 /  (N + he))));
	
	while(abs(phiOld - phiNew) > 1e-10){
		phiOld = phiNew;
		w = sqrt(1 - GRS_e_2 * pow(sin(phiOld), 2));
		N = GRS_a / w;
		he = p * cos(phiOld) + Z * sin(phiOld) - GRS_a * w;
		phiNew = atan2((Z / p), ((1 - N * GRS_e_2 / (N + he))));
	}
	
	latitude = phiNew;
}

void GPSGeodesie::Geographique_2_ECEF(double lng, double lat, double alt, double& x, double& y, double& z){
	double w = sqrt(1 - GRS_e_2 * pow(sin(lat), 2)); //TODO: optimiser avec fastInvSqrt
	double N = GRS_a / w;
	
	x = (N + alt) * cos(lat) * cos(lng);
	y = (N + alt) * cos(lat) * sin(lng);
	z = (N * (1 - GRS_e_2) + alt) * sin(lat);
}

void GPSGeodesie::ENU_2_ECEF(double e, double n, double u, double& x, double& y, double& z, double lon0, double lat0, double he0)
{
	double slat = sin(lat0);
	double clat = cos(lat0);
	double slon = sin(lon0);
	double clon = cos(lon0);
	
	MatriceSqr<double, 3> C;
	
	C[0][0] = -slon;
	C[1][0] = clon;
	C[2][0] = 0;
	
	C[0][1] = -clon * slat;
	C[1][1] = -slon * slat;
	C[2][1] = clat;
	
	C[0][2] = clon * clat;
	C[1][2] = slon * clat;
	C[2][2] = slat;
	
	Matrice<double, 3, 1> mat, res;
	
	mat[0][0] = e;
	mat[1][0] = n;
	mat[2][0] = u;
	
	res = C * mat;
	
	x = res[0][0];
	y = res[1][0];
	z = res[2][0];
	
	double x0, y0, z0;
	Geographique_2_ECEF(lon0, lat0, he0, x0, y0, z0);
	
	x += x0;
	y += y0;
	z += z0;
}

void GPSGeodesie::ECEF_2_ENU(double x, double y, double z, double& e, double& n, double& u, double lon0, double lat0, double he0){
	double x0, y0, z0;
	Geographique_2_ECEF(lon0, lat0, he0, x0, y0, z0);
	
	x -= x0;
	y -= y0;
	z -= z0;
	
	double slat = sin(lat0);
	double clat = cos(lat0);
	double slon = sin(lon0);
	double clon = cos(lon0);
	
	Matrice<double, 3, 3> C;
	
	C[0][0] = -slon;
	C[0][1] = clon;
	C[0][2] = 0;
	
	C[1][0] = -clon * slat;
	C[1][1] = -slon * slat;
	C[1][2] = clat;
	
	C[2][0] = clon * clat;
	C[2][1] = slon * clat;
	C[2][2] = slat;
	
	Matrice<double, 3, 1> mat, res;
	
	mat[0][0] = x;
	mat[1][0] = y;
	mat[2][0] = z;
	
	res = C * mat;
	
	e = res[0][0];
	n = res[1][0];
	u = res[2][0];
}

void GPSGeodesie::XYZ2LatLng(double X, double Y, double Z, double lat0, double lng0, double alt0, double& lat, double& lng, double& alt){
	/// Convert from ECEF two ENU.
	/// @param[in] lon0 Longitude of the origin in radian.
	/// @param[in] lat0 Latitude of the origin in radian.
	/// @param[in] he0 Height of the origin in radian.
	double x = 0;
	double y = 0;
	double z = 0;
	ENU_2_ECEF(X, Y, Z, x, y, z, toRadians(lng0), toRadians(lat0), alt0);
	
	/// Convert from geographique to ECEF.
	/// @param[in] longitude Longitude in radian.
	/// @param[in] latitude Latitude in radian.
	/// @param[in] he Height in meter.
	ECEF_2_Geographique(x, y, z, lng, lat, alt);
	
	lng = toDegrees(lng);
	lat = toDegrees(lat);
}

////////////////////////////////////////////////////////////////////////

void GPSGeodesie::LatLng2XYZ(double lat, double lng, double alt, double lat0, double lng0, double alt0, double& X, double& Y, double& Z){
	/// Convert from geographique to ECEF.
	/// @param[in] longitude Longitude in radians.
	/// @param[in] latitude Latitude in radians.
	/// @param[in] he Height in meters.
	double x = 0;
	double y = 0;
	double z = 0;
	Geographique_2_ECEF(toRadians(lng), toRadians(lat), alt, x, y, z);
	
	/// Convert from ECEF two ENU.
	/// @param[in] lon0 Longitude of the origin in radians.
	/// @param[in] lat0 Latitude of the origin in radians.
	/// @param[in] he0 Height of the origin in meters.
	ECEF_2_ENU(x, y, z, X, Y, Z, toRadians(lng0), toRadians(lat0), alt0);
}
