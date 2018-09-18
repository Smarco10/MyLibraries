//
//  GPSGeodesie.hpp
//  test
//
//  Created by Marc-Antoine MARTIN on 01/04/2016.
//  Copyright © 2016 Marc-Antoine MARTIN. All rights reserved.
//

#ifndef GPSGeodesie_hpp
#define GPSGeodesie_hpp

#include <cstdio>
#include <cmath>

namespace GPSGeodesie {
#define toRadians(deg) ((double)(deg) * 0.01745329251994329576923690768489) // M_PI / 180.0
#define toDegrees(rad) ((double)(rad) * 57.295779513082320876798154814105) // 180.0 / M_PI
	
	/// Converts Cartesian (x, y) coordinates to polar coordinates (r, theta)
	template <typename _T1, typename _T2>
	void cartesianToPolar(const _T1 x, const _T1 y, _T2 & r, _T2 & theta) {
		r = std::sqrt(x*x + y*y);
		theta = std::atan2(x, y);
	}
	
	/// Converts polar coordinates (r, theta) to Cartesian (x, y) coordinates
	template <typename _T1, typename _T2>
	void polarToCartesian(const _T1 r, const _T1 theta, _T2 & x, _T2 & y) {
		x = r * std::cos(theta);
		y = r * std::sin(theta);
	}
	
	/// Converts Cartesian (x, y, z) coordinates to spherical coordinates (r, theta, phi)
	/// Angles expressed in radians.
	template <typename _T1, typename _T2>
	void cartesianToSpherical(const _T1 x, const _T1 y, const _T1 z, _T2 & r, _T2 & theta, _T2 & phi) {
		r = std::sqrt(x*x + y*y + z*z);
		theta = std::acos(z / r);
		phi = std::atan2(y, x);
	}
	
	/// Converts spherical coordinates (r, theta, phi) to Cartesian (x, y, z) coordinates.
	/// Angles expressed in radians.
	template <typename _T1, typename _T2>
	void sphericalToCartesian(const _T1 r, const _T1 theta, const _T1 phi, _T2 & x, _T2 & y, _T2 & z) {
		x = r * std::sin(theta) * std::cos(phi);
		y = r * std::sin(theta) * std::sin(phi);
		z = r * std::cos(theta);
	}
	
	void ECEF_2_Geographique(double X, double Y, double Z, double& longitude, double& latitude, double& he);
	void Geographique_2_ECEF(double lng, double lat, double alt, double& x, double& y, double& z);
	
	void ENU_2_ECEF(double e, double n, double u, double& x, double& y, double& z, double lon0, double lat0, double he0);
	void ECEF_2_ENU(double x, double y, double z, double& e, double& n, double& u, double lon0, double lat0, double he0);
	
	void XYZ2LatLng(double X, double Y, double Z, double lat0, double lng0, double alt0, double& lat, double& lng, double& alt);
	void LatLng2XYZ(double lat, double lng, double alt, double lat0, double lng0, double alt0, double& X, double& Y, double& Z);

}

#endif /* GPSGeodesie_hpp */
