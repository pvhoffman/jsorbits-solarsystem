// jsorbits-solarsystem.js
//
// Copyright (c) 2011, Paul V. Hoffman.
// All rights reserved.

var solarSystemOrbits = {
	_eccentricAnomaly : function(M, e) {
		var E0 = M + e * Math.sin(M) * (1.0 + e * Math.cos(M));
		var E1 = E0 - (E0 - e * Math.sin(E0) - M) / (1.0 - e * Math.cos(E0));
		var D;
		var I = 0;

		do {
			E0 = E1;
			E1 = E0 - (E0 - e * Math.sin(E0) - M) / (1 - e * Math.cos(E0));
			D = Math.abs((E1 - E0));
			I = I + 1;
		} while(D > 0.00001 && I < 10);

		if(I == 10 && D > 0.00001) {
			alert("Integration failed, life sucks...");
		}
		return E1;
	},
	_truncate : function(d) {
		if(d < 0) return Math.ceil(d);
		return Math.floor(d);
	},
	_normalizeAngle : function(a){
		return a - Math.floor(a / (Math.PI * 2.0)) * (Math.PI * 2.0);
	},
	_gmst0Hours : function(epoch){
		var M = 6.214192441848251 + 0.01720196961933223 * epoch;
		var w = 4.938241566909764 + 8.219366312879496e-07 * epoch;
		var L = this._normalizeAngle((M + w));
		L = L / ((Math.PI * 2.0) / 24.0);
		L = L + 12.0;
		return L;
	},
	_sidtimeHours : function(epoch, utchours, longitude){ //coordinates) {
		var g = this._gmst0Hours(epoch);
		var r = g + utchours + longitude / 15.0;
		if (r < 0.0) r = r + 24.0;
		if (r > 24.0) r = r - 24.0;
		return r;
	},
	_hourAngleDegrees : function (epoch, utchours, longitude, ra){ //coordinates, longitude) {
		var s = this._sidtimeHours(epoch, utchours, longitude);
		s = s * ((Math.PI * 2.0) / 24.0);
		var r = s - ra;
		r = this._normalizeAngle(r);
		return r;
	},
	_rectangularToSpherical : function(x, y, z){
		var e = 0.0001;
		var r = Math.sqrt(x*x + y*y + z*z);
		var ra = 0.0;
		var decl = 0.0;
		if (Math.abs(x) < e &&  Math.abs(y) < e) {
			decl = Math.atan2(z, Math.sqrt(x*x + y*y));
		} else {
			ra = Math.atan2(y, x);
			decl = Math.asin(z / r);
		}
		return {'r' : r, 'ra' : ra, 'decl' : decl};
	},
	epochFromDate : function(date) {
		var Y = date.getFullYear();
		var M = date.getMonth() + 1;
		var D = date.getDate();
		var h = date.getHours();
		var m = date.getMinutes();
		var s = date.getSeconds();
		var d = 367 * Y - this._truncate((7 * (Y + this._truncate(((M + 9) / 12)))) / 4) + this._truncate((275 * M) / 9) + D - 730530;
		d = d + h / 24.0 + m / 1440.0 + s / 86400.0;
		//return -3543;
		return d;
	},
	orbitalElementsMercury : {
		id : 'mercury',
		longitudeOfAscendingNode : function(epoch) {
			return 0.84354031676913532 + 5.6651118591708344e-07 * epoch;
		},
		inclination              : function(epoch) {
			return 0.1222550781144468 + 8.7266462599716477e-10 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 0.5083114366800805 + 1.7705318063931279e-07 * epoch;
		},
		semiMajorAxis            : function(epoch) {
			return 0.387098;
		},
		eccentricity             : function(epoch) {
			return 0.205635 + 5.59E-10 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 2.9436059939020605 + 0.0714247100149078 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			return 0.0;
		},
		latitudePerturbations    : function(epoch) {
			return 0.0;
		},
		distancePerturbations : function(epoch) {
			return 0.0;
		},
		parallax : function(r) {
			var par = 4.26345151167726e-05 / r;
			return par;
		},
		magnitude : function(r, R, phaseangle) {
			//log10 =  Math.log(arg) / Math.LN10;
			var FV   =  phaseangle * (180.0/Math.PI);
			var res  = -0.36 + 5.0 * (Math.log(r*R) / Math.LN10) + 0.027 * FV + 2.2E-13 * Math.pow(FV, 6);
			return res;
		},
		data : function(epoch, utchours, geoLat, geoLon) {
			var self = solarSystemOrbits;
			var c = self._planetaryCoordinates(epoch, this, geoLat, geoLon, utchours);
			var e = self._planetaryElongationPhaseAngleMagnitude(this, c, epoch);
			var res = {
				'c' : c,
				'e' : e
			};
			return res;
		}
	},
	orbitalElementsVenus : {
		id : 'venus',
		longitudeOfAscendingNode : function(epoch) {
			return 1.33832 + 4.30381e-07 * epoch;
		},
		inclination              : function(epoch) {
			return 0.0592469 + 4.79966e-10 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 0.958029 + 2.41508e-07 * epoch;
		},
		semiMajorAxis            : function(epoch) {
			return 0.723330;
		},
		eccentricity             : function(epoch) {
			return 0.006773 - 1.302E-9 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 0.837849 + 0.0279624 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			return 0.0;
		},
		latitudePerturbations    : function(epoch) {
			return 0.0;
		},
		distancePerturbations : function(epoch) {
			return 0.0;
		},
		parallax : function(r) {
			var par = 4.26345151167726e-05 / r;
			return par;
		},
		magnitude : function(r, R, phaseangle) {
			var FV  =  phaseangle * (180.0/Math.PI);
			var res = -4.34 + 5 * (Math.log(r*R) / Math.LN10) + 0.013 * FV + 4.2E-7  * Math.pow(FV, 3);
			return res;
		},
		data : function(epoch, utchours, geoLat, geoLon) {
			var self = solarSystemOrbits;
			var c = self._planetaryCoordinates(epoch, this, geoLat, geoLon, utchours);
			var e = self._planetaryElongationPhaseAngleMagnitude(this, c, epoch);
			var res = {
				'c' : c,
				'e' : e
			};
			return res;
		}
	},
	orbitalElementsSun : {
		id : 'sun',
		longitudeOfAscendingNode : function(epoch) {
			return 0.0;
		},
		inclination              : function(epoch) {
			return 0.0;
		},
		argumentOfPerihelion     : function(epoch) {
			return 4.9382415669097641 + 8.2193663128794961e-07 * epoch;
		},
		semiMajorAxis            : function(epoch){
			return 1.0;
		},
		eccentricity             : function(epoch) {
			return 0.016709 - 1.151E-9 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 6.214192441848251 + 0.017201969619332229 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			return 0.0;
		},
		latitudePerturbations    : function(epoch) {
			return 0.0;
		},
		distancePerturbations : function(epoch) {
			return 0.0;
		},
		parallax : function(r) {
			var par = 4.26345151167726e-05 / r;
			return par;
		},
		magnitude : function(r, R, FV) {
			return 0.0;
		},
		data : function(epoch, utchours, geoLat, geoLon) {
			var self = solarSystemOrbits;
			var N = self._normalizeAngle(this.longitudeOfAscendingNode(epoch));
			var i = self._normalizeAngle(this.inclination(epoch));
			var w = self._normalizeAngle(this.argumentOfPerihelion(epoch));
			var a = this.semiMajorAxis(epoch);
			var e = this.eccentricity(epoch);
			var M = self._normalizeAngle(this.meanAnomaly(epoch));
			var E = self._eccentricAnomaly(M, e);
			var x = a * (Math.cos(E) - e);
			var y = a * (Math.sin(E) * Math.sqrt(1.0 - e*e));
			var z = 0;
			var r = Math.sqrt(x*x + y*y);
			var v = self._normalizeAngle(Math.atan2(y, x));

			x = r * (Math.cos(N) * Math.cos(v+w) - Math.sin(N) * Math.sin(v+w) * Math.cos(i));
			y = r * (Math.sin(N) * Math.cos(v+w) + Math.cos(N) * Math.sin(v+w) * Math.cos(i));
			z = r * Math.sin(v+w) * Math.sin(i);

			var lon = self._normalizeAngle(Math.atan2(y, x));
			var lat = self._normalizeAngle(Math.asin(z / r));

			var heliocentric = {
				'x' : r * Math.cos(lon) * Math.cos(lat),
				'y' : r * Math.sin(lon) * Math.cos(lat),
				'z' : r * Math.sin(lat)
			};

			var oblecl = 0.4090929593627069 - 6.218608124855796e-09 * epoch;

			var geocentric = {
				'x' : heliocentric.x,
				'y' : heliocentric.y * Math.cos(oblecl) - heliocentric.z * Math.sin(oblecl),
				'z' : heliocentric.y * Math.sin(oblecl) + heliocentric.z * Math.cos(oblecl)
			};

			var spherical = self._rectangularToSpherical(geocentric.x, geocentric.y, geocentric.z);
			spherical.ra  = self._normalizeAngle(spherical.ra);
			var ha = self._hourAngleDegrees(epoch, utchours, geoLon, spherical.ra);
			var par = this.parallax(spherical.r);
			// latitude adjustment
			var gclat = geoLat - 0.1924 * Math.sin((2.0*geoLat) * (Math.PI / 180.0));
			gclat = gclat * (Math.PI / 180.0);
			// distance from the center of the earth
			var rho = 0.99833 + 0.00167 * Math.cos((2.0*geoLat) * (Math.PI / 180.0));
			// aux angle
			var g = Math.atan( Math.tan(gclat) / Math.cos(ha));
			g = self._normalizeAngle(g);
			// RA adjustment
			var topRA   = par * rho * Math.cos(gclat) * Math.sin(ha) / Math.cos(spherical.decl);
			// DECL adjustment
			var topDecl = par * rho * Math.sin(gclat) * Math.sin(g - spherical.decl) / Math.sin(g);
			// adjust the geocentric spherical coordinates to topocentric
			spherical.ra = spherical.ra - topRA;
			spherical.decl = spherical.decl - topDecl;
			spherical.ra = self._normalizeAngle(spherical.ra);
			// convert to rectnagle.
			// x axis points to the celestial equator in the south
			x = Math.cos(ha) * Math.cos(spherical.decl);
			// y axis points to the horizon in the west
			y = Math.sin(ha) * Math.cos(spherical.decl);
			// z axis points to the north celestial pole
			z = Math.sin(spherical.decl);
			// rotate alon the y axis (pointing east/west) so z points at the zenith
			var topocentric = {
				'x' : x * Math.sin(gclat) - z * Math.cos(gclat),
				'y' : y,
				'z' : x * Math.cos(gclat) + z * Math.sin(gclat)
			};

			var c = {
				'heliocentric' : heliocentric,
				'geocentric'   : geocentric,
				'topocentric'  : topocentric,
				'hourangle'    : ha
			};
			var e = self._planetaryElongationPhaseAngleMagnitudeZero(c, epoch);
			var res = {
				'c' : c,
				'e' : e
			};
			return res;

		}
	},
	orbitalElementsMars : {
		id : 'mars',
		longitudeOfAscendingNode : function(epoch) {
			return 0.86494 + 3.68406e-07 * epoch;
		},
		inclination              : function(epoch) {
			return 0.0322834 - 3.10669e-10 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 5.0004 + 5.11313e-07 * epoch;
		},
		semiMajorAxis            : function(epoch) {
			return 1.523688;
		},
		eccentricity             : function(epoch) {
			return 0.093405 + 2.516E-9 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 0.324668 + 0.00914589 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			return 0.0;
		},
		latitudePerturbations    : function(epoch) {
			return 0.0;
		},
		distancePerturbations : function(epoch) {
			return 0.0;
		},
		parallax : function(r) {
			var par = 4.26345151167726e-05 / r;
			return par;
		},
		magnitude : function(r, R, phaseangle) {
			var FV  =  phaseangle * (180.0/Math.PI);
			var res = -1.51 + 5.0 * (Math.log(r*R) / Math.LN10) + 0.016 * FV;
			return res;
		},
		data : function(epoch, utchours, geoLat, geoLon) {
			var self = solarSystemOrbits;
			var c = self._planetaryCoordinates(epoch, this, geoLat, geoLon, utchours);
			var e = self._planetaryElongationPhaseAngleMagnitude(this, c, epoch);
			var res = {
				'c' : c,
				'e' : e
			};
			return res;
		}
	},
	orbitalElementsJupiter : {
		id : 'jupiter',
		longitudeOfAscendingNode : function(epoch) {
			return 1.75326 + 4.83201e-07 * epoch;
		},
		inclination              : function(epoch) {
			return 0.0227416 - 2.71748e-09 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 4.78007 + 2.87115e-07 * epoch;
		},
		semiMajorAxis            : function(epoch){
			return 5.20256;
		},
		eccentricity             : function(epoch) {
			return 0.048498 + 4.469E-9 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 0.347233 + 0.00145011 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			var self = solarSystemOrbits;
			var Mj = self._planetaryPerturbationsMj(epoch);
			var Ms = self._planetaryPerturbationsMs(epoch);
			var p = 0.0;
			p = p + -0.00579449311662 * Math.sin(2*Mj - 5*Ms - 1.17984257435);
			p = p + -0.000977384381117 * Math.sin(2*Mj - 2*Ms + 0.366519142919);
			p = p + 0.000733038285838 * Math.sin(3*Mj - 5*Ms + 0.366519142919);
			p = p + -0.000628318530718 * Math.sin(Mj - 2*Ms);
			p = p + 0.000383972435439 * Math.cos(Mj - Ms);
			p = p + 0.000401425727959 * Math.sin(2*Mj - 3*Ms + 0.907571211037);
			p = p + -0.000279252680319 * Math.sin(Mj - 5*Ms - 1.20427718388);
			return p;
		},
		latitudePerturbations    : function(epoch) {
			return 0.0;
		},
		distancePerturbations : function(epoch) {
			return 0.0;
		},
		parallax : function(r) {
			var par = 4.26345151167726e-05 / r;
			return par;
		},
		magnitude : function(r, R, phaseangle) {
			var FV  =  phaseangle * (180.0/Math.PI);
			var res = -9.25 + 5.0 * (Math.log(r*R) / Math.LN10) + 0.014 * FV;
			return res;
		},
		data : function(epoch, utchours, geoLat, geoLon) {
			var self = solarSystemOrbits;
			var c = self._planetaryCoordinates(epoch, this, geoLat, geoLon, utchours);
			var e = self._planetaryElongationPhaseAngleMagnitude(this, c, epoch);
			var res = {
				'c' : c,
				'e' : e
			};
			return res;
		}
	},
	orbitalElementsSaturn : {
		id : 'saturn',
		longitudeOfAscendingNode : function(epoch) {
			return 1.9838 + 4.17099e-07 * epoch;
		},
		inclination              : function(epoch) {
			return 0.0434343 - 1.8867e-09 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 5.92354 + 5.19516e-07 * epoch;
		},
		semiMajorAxis            : function(epoch) {
			return 9.55475;
		},
		eccentricity             : function(epoch) {
			return 0.055546 - 9.499E-9 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 5.53212 + 0.000583712 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			var self = solarSystemOrbits;
			var Mj = self._planetaryPerturbationsMj(epoch);
			var Ms = self._planetaryPerturbationsMs(epoch);
			var p = 0.0;
			p = p + 0.0141720735262 * Math.sin(2*Mj - 5*Ms - 1.17984257435);
			p = p + -0.00399680398707 * Math.cos(2*Mj - 4*Ms - 0.0349065850399);
			p = p + 0.00207694180987 * Math.sin(Mj - 2*Ms - 0.0523598775598);
			p = p + 0.000802851455917 * Math.sin(2*Mj - 6*Ms - 1.20427718388);
			p = p + 0.000244346095279 * Math.sin(Mj - 3*Ms + 0.558505360638);
			return p;
		},
		latitudePerturbations    : function(epoch) {
			var self = solarSystemOrbits;
			var Mj = self._planetaryPerturbationsMj(epoch);
			var Ms = self._planetaryPerturbationsMs(epoch);
			var p = 0.0;
			p = p + -0.000349065850399 * Math.cos(2*Mj - 4*Ms - 0.0349065850399);
			p = p + 0.000314159265359 * Math.sin(2*Mj - 6*Ms - 0.855211333477);
			return p;

		},
		distancePerturbations : function(epoch) {
			return 0.0;
		},
		parallax : function(r) {
			var par = 4.26345151167726e-05 / r;
			return par;
		},
		magnitude : function(r, R, phaseangle) {
			var FV  = phaseangle * (180.0/Math.PI);
			var ring_magn = 0.0;
			var res = -9.0  + 5.0 *(Math.log(r*R) / Math.LN10) + 0.044 * FV + ring_magn;
			return res;
		},
		data : function(epoch, utchours, geoLat, geoLon) {
			var self = solarSystemOrbits;
			var c = self._planetaryCoordinates(epoch, this, geoLat, geoLon, utchours);
			var e = self._planetaryElongationPhaseAngleMagnitude(this, c, epoch);
			var res = {
				'c' : c,
				'e' : e
			};
			return res;
		}

	},
	orbitalElementsUranus : {
		id : 'uranus',
		longitudeOfAscendingNode : function(epoch) {
			return 1.29155 + 2.43962e-07 * epoch;
		},
		inclination              : function(epoch) {
			return 0.0134966 + 3.31613e-10 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 1.68706 + 5.3346e-07 * epoch;
		},
		semiMajorAxis            : function(epoch){
			return 19.18171 - 1.55E-8 * epoch;
		},
		eccentricity             : function(epoch) {
			return 0.047318 + 7.45E-9 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 2.48867 + 0.000204654 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			var self = solarSystemOrbits;
			var Mj = self._planetaryPerturbationsMj(epoch);
			var Ms = self._planetaryPerturbationsMs(epoch);
			var Mu = self._planetaryPerturbationsMu(epoch);
			var p = 0.0;
			p = p + 0.000698131700798 *  Math.sin(Ms - 2*Mu + 0.10471975512);
			p = p + 0.000610865238198 *  Math.sin(Ms - 3*Mu + 0.575958653158);
			p = p + -0.000261799387799 * Math.sin(Mj - Mu + 0.349065850399);
			return p;
		},
		latitudePerturbations    : function(epoch) {
			return 0.0;
		},
		distancePerturbations : function(epoch) {
			return 0.0;
		},
		parallax : function(r) {
			var par = 4.26345151167726e-05 / r;
			return par;
		},
		magnitude : function(r, R, phaseangle) {
			var FV  = phaseangle * (180.0/Math.PI);
			var res = -7.15 + 5.0 * (Math.log(r*R) / Math.LN10) + 0.001 * FV;
			return res;
		},
		data : function(epoch, utchours, geoLat, geoLon) {
			var self = solarSystemOrbits;
			var c = self._planetaryCoordinates(epoch, this, geoLat, geoLon, utchours);
			var e = self._planetaryElongationPhaseAngleMagnitude(this, c, epoch);
			var res = {
				'c' : c,
				'e' : e
			};
			return res;
		}
	},
	orbitalElementsNeptune : {
		id : 'neptune',
		longitudeOfAscendingNode : function(epoch) {
			return 2.30001 + 5.26618e-07 * epoch;
		},
		inclination              : function(epoch) {
			return 0.0308923 - 4.45059e-09 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 4.76206 - 1.05191e-07 * epoch;
		},
		semiMajorAxis            : function(epoch){
			return 30.05826 + 3.313E-8 * epoch;
		},
		eccentricity             : function(epoch) {
			return 0.008606 + 2.15E-9 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 4.54217 + 0.000104635 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			return 0.0;
		},
		latitudePerturbations    : function(epoch) {
			return 0.0;
		},
		distancePerturbations : function(epoch) {
			return 0.0;
		},
		parallax : function(r) {
			var par = 4.26345151167726e-05 / r;
			return par;
		},
		magnitude : function(r, R, phaseangle) {
			var FV  = phaseangle * (180.0/Math.PI);
			var res = -6.90 + 5.0 * (Math.log(r*R) / Math.LN10) + 0.001 * FV;
			return res;
		},
		data : function(epoch, utchours, geoLat, geoLon) {
			var self = solarSystemOrbits;
			var c = self._planetaryCoordinates(epoch, this, geoLat, geoLon, utchours);
			var e = self._planetaryElongationPhaseAngleMagnitude(this, c, epoch);
			var res = {
				'c' : c,
				'e'   : e
			};
			return res;
		}
	},
	orbitalElementsEarthMoon : {
		id : 'moon-earth',
		longitudeOfAscendingNode : function(epoch) {
			return 2.183804829314361 - 0.0009242183063049012 * epoch;
		},
		inclination              : function(epoch) {
			return 0.08980417133211624;
		},
		argumentOfPerihelion     : function(epoch) {
			return 5.551253560087733 + 0.0028685764238964994 * epoch;
		},
		semiMajorAxis            : function(epoch){
			return 60.2666;
		},
		eccentricity             : function(epoch) {
			return 0.054900;
		},
		meanAnomaly              : function(epoch) {
			return 2.0135060728802663 + 0.22802714374305486 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			var self = solarSystemOrbits;
			//Sun's  mean longitude
			var Ls = self.orbitalElementsSun.argumentOfPerihelion(epoch) + self.orbitalElementsSun.meanAnomaly(epoch);
			//Moon's mean longitude
			var Lm = this.longitudeOfAscendingNode(epoch) + this.argumentOfPerihelion(epoch) + this.meanAnomaly(epoch);
			//Moon's mean anomaly
			var Mm = this.meanAnomaly(epoch);
			// Sun's  mean anomaly
			var Ms = self.orbitalElementsSun.meanAnomaly(epoch);
			//Moon's mean elongation
			var D = Lm - Ls;
			//Moon's argument of latitude
			var F = Lm - this.longitudeOfAscendingNode(epoch);
			var p = 0.0;

			Ls = self._normalizeAngle(Ls);
			Lm = self._normalizeAngle(Lm);
			Mm = self._normalizeAngle(Mm);
			Ms = self._normalizeAngle(Ms);
			D  = self._normalizeAngle(D);
			F  = self._normalizeAngle(F);

			p = p + -0.0222354946704   * Math.sin(Mm - 2*D);
			p = p + 0.0114842664781    * Math.sin(2*D);
			p = p + -0.00324631240871  * Math.sin(Ms);
			p = p + -0.00102974425868  * Math.sin(2*Mm - 2*D);
			p = p + -0.000994837673637 * Math.sin(Mm - 2*D + Ms);
			p = p + 0.000925024503557  * Math.sin(Mm + 2*D);
			p = p + 0.000802851455917  * Math.sin(2*D - Ms);
			p = p + 0.000715584993318  * Math.sin(Mm - Ms);
			p = p + -0.000610865238198 * Math.sin(D);
			p = p + -0.000541052068118 * Math.sin(Mm + Ms);
			p = p + -0.000261799387799 * Math.sin(2*F - 2*D);
			p = p + 0.000191986217719  * Math.sin(Mm - 4*D);
			return p;
		},
		latitudePerturbations    : function(epoch) {
			var self = solarSystemOrbits;
			//Sun's  mean longitude
			var Ls = self.orbitalElementsSun.argumentOfPerihelion(epoch) + self.orbitalElementsSun.meanAnomaly(epoch);
			//Moon's mean longitude
			var Lm = this.longitudeOfAscendingNode(epoch) + this.argumentOfPerihelion(epoch) + this.meanAnomaly(epoch);
			//Moon's mean anomaly
			var Mm = this.meanAnomaly(epoch);
			// Sun's  mean anomaly
			var Ms = self.orbitalElementsSun.meanAnomaly(epoch);
			//Moon's mean elongation
			var D = Lm - Ls;
			//Moon's argument of latitude
			var F = Lm - this.longitudeOfAscendingNode(epoch);
			var p = 0.0;

			Ls = self._normalizeAngle(Ls);
			Lm = self._normalizeAngle(Lm);
			Mm = self._normalizeAngle(Mm);
			Ms = self._normalizeAngle(Ms);
			D  = self._normalizeAngle(D);
			F  = self._normalizeAngle(F);

			p = p + -0.00301941960595 * Math.sin(F - 2*D);
			p = p + -0.000959931088597 * Math.sin(Mm - F - 2*D);
			p = p + -0.000802851455917 * Math.sin(Mm + F - 2*D);
			p = p + 0.000575958653158 * Math.sin(F + 2*D);
			p = p + 0.000296705972839 * Math.sin(2*Mm + F);
			return p;
		},
		distancePerturbations : function(epoch) {
			var self = solarSystemOrbits;
			//Sun's  mean longitude
			var Ls = self.orbitalElementsSun.argumentOfPerihelion(epoch) + self.orbitalElementsSun.meanAnomaly(epoch);
			//Moon's mean longitude
			var Lm = this.longitudeOfAscendingNode(epoch) + this.argumentOfPerihelion(epoch) + this.meanAnomaly(epoch);
			//Moon's mean anomaly
			var Mm = this.meanAnomaly(epoch);
			// Sun's  mean anomaly
			var Ms = self.orbitalElementsSun.meanAnomaly(epoch);
			//Moon's mean elongation
			var D = Lm - Ls;
			//Moon's argument of latitude
			var F = Lm - this.longitudeOfAscendingNode(epoch);
			var p = 0.0;

			Ls = self._normalizeAngle(Ls);
			Lm = self._normalizeAngle(Lm);
			Mm = self._normalizeAngle(Mm);
			Ms = self._normalizeAngle(Ms);
			D  = self._normalizeAngle(D);
			F  = self._normalizeAngle(F);

			p = p + -0.58 * Math.cos(Mm - 2*D);
			p = p + -0.46 * Math.cos(2*D);

			return p;
		},
		parallax : function(r) {
			var par = Math.asin( 1.0 / r );
			return par;
		},
		magnitude : function(r, R, FV) {
			return 0.0;
		},
		data : function(epoch, utchours, geoLat, geoLon) {
			var self = solarSystemOrbits;
			var N = self._normalizeAngle(this.longitudeOfAscendingNode(epoch));
			var i = self._normalizeAngle(this.inclination(epoch));
			var w = self._normalizeAngle(this.argumentOfPerihelion(epoch));
			var a = this.semiMajorAxis(epoch);
			var e = this.eccentricity(epoch);
			var M = self._normalizeAngle(this.meanAnomaly(epoch));
			var E = self._eccentricAnomaly(M, e);
			var x = a * (Math.cos(E) - e);
			var y = a * (Math.sin(E) * Math.sqrt(1.0 - e*e));
			var z = 0;
			var r = Math.sqrt(x*x + y*y);
			var v = self._normalizeAngle(Math.atan2(y, x));

			x = r * (Math.cos(N) * Math.cos(v+w) - Math.sin(N) * Math.sin(v+w) * Math.cos(i));
			y = r * (Math.sin(N) * Math.cos(v+w) + Math.cos(N) * Math.sin(v+w) * Math.cos(i));
			z = r * Math.sin(v+w) * Math.sin(i);

			var lon = self._normalizeAngle(Math.atan2(y, x));
			var lat = self._normalizeAngle(Math.asin(z / r));

			var heliocentric = {
				'x' : r * Math.cos(lon) * Math.cos(lat),
				'y' : r * Math.sin(lon) * Math.cos(lat),
				'z' : r * Math.sin(lat)
			};

			var oblecl = 0.4090929593627069 - 6.218608124855796e-09 * epoch;

			var geocentric = {
				'x' : heliocentric.x,
				'y' : heliocentric.y * Math.cos(oblecl) - heliocentric.z * Math.sin(oblecl),
				'z' : heliocentric.y * Math.sin(oblecl) + heliocentric.z * Math.cos(oblecl)
			};

			var spherical = self._rectangularToSpherical(geocentric.x, geocentric.y, geocentric.z);
			spherical.ra  = self._normalizeAngle(spherical.ra);
			var ha = self._hourAngleDegrees(epoch, utchours, geoLon, spherical.ra);
			var par = this.parallax(spherical.r);
			// latitude adjustment
			var gclat = geoLat - 0.1924 * Math.sin((2.0*geoLat) * (Math.PI / 180.0));
			gclat = gclat * (Math.PI / 180.0);
			// distance from the center of the earth
			var rho = 0.99833 + 0.00167 * Math.cos((2.0*geoLat) * (Math.PI / 180.0));
			// aux angle
			var g = Math.atan( Math.tan(gclat) / Math.cos(ha));
			g = self._normalizeAngle(g);
			// RA adjustment
			var topRA   = par * rho * Math.cos(gclat) * Math.sin(ha) / Math.cos(spherical.decl);
			// DECL adjustment
			var topDecl = par * rho * Math.sin(gclat) * Math.sin(g - spherical.decl) / Math.sin(g);
			// adjust the geocentric spherical coordinates to topocentric
			spherical.ra = spherical.ra - topRA;
			spherical.decl = spherical.decl - topDecl;
			spherical.ra = self._normalizeAngle(spherical.ra);
			// convert to rectnagle.
			// x axis points to the celestial equator in the south
			x = Math.cos(ha) * Math.cos(spherical.decl);
			// y axis points to the horizon in the west
			y = Math.sin(ha) * Math.cos(spherical.decl);
			// z axis points to the north celestial pole
			z = Math.sin(spherical.decl);
			// rotate alon the y axis (pointing east/west) so z points at the zenith
			var topocentric = {
				'x' : x * Math.sin(gclat) - z * Math.cos(gclat),
				'y' : y,
				'z' : x * Math.cos(gclat) + z * Math.sin(gclat)
			};

			var c = {
				'heliocentric' : heliocentric,
				'geocentric'   : geocentric,
				'topocentric'  : topocentric,
				'hourangle'    : ha
			};
			var e = self._lunarElongationPhaseAngleMagnitude(this, c, epoch);
			var res = {
				'c' : c,
				'e' : e
			};
			return res;

		}
	},
	_planetaryCoordinates : function(epoch, body, geoLat, geoLon, utchours) {
		// return the rectangular heliocentric and geocentric coordinates
		var N = this._normalizeAngle(body.longitudeOfAscendingNode(epoch));
		var i = this._normalizeAngle(body.inclination(epoch));
		var w = this._normalizeAngle(body.argumentOfPerihelion(epoch));
		var a = body.semiMajorAxis(epoch);
		var e = body.eccentricity(epoch);
		var M = this._normalizeAngle(body.meanAnomaly(epoch));
		var E = this._eccentricAnomaly(M, e);


		var perturbationsLat  = body.latitudePerturbations(epoch);
		var perturbationsLon  = body.longitudePerturbations(epoch);
		var perturbationsDist = body.distancePerturbations(epoch);

		var x = a * (Math.cos(E) - e);
		var y = a * (Math.sin(E) * Math.sqrt(1.0 - e*e));
		var z = 0;

		var r = Math.sqrt(x*x + y*y);
		var v = Math.atan2(y, x);
		v = this._normalizeAngle(v);

		x = r * (Math.cos(N) * Math.cos(v+w) - Math.sin(N) * Math.sin(v+w) * Math.cos(i));
		y = r * (Math.sin(N) * Math.cos(v+w) + Math.cos(N) * Math.sin(v+w) * Math.cos(i));
		z = r * Math.sin(v+w) * Math.sin(i);

		var lon = Math.atan2(y, x);
		var lat = Math.asin(z / r);

		lon = lon + perturbationsLon;
		lat = lat + perturbationsLat;
		r   = r + perturbationsDist;

		lat = this._normalizeAngle(lat);
		lon = this._normalizeAngle(lon);

		var hx = r * Math.cos(lon) * Math.cos(lat);
		var hy = r * Math.sin(lon) * Math.cos(lat);
		var hz = r * Math.sin(lat);

		x = hx;
		y = hy;
		z = hz;

		// now get the geocentric
		var sun = this.orbitalElementsSun.data(epoch, utchours, geoLat, geoLon);
		x = sun.coords.heliocentric.x + x;
		y = sun.coords.heliocentric.y + y;
		z = sun.coords.heliocentric.z + z;

		var oblecl = 0.4090929593627069 - 6.218608124855796e-09 * epoch;

		var gx = x;
		var gy = y * Math.cos(oblecl) - z * Math.sin(oblecl);
		var gz = y * Math.sin(oblecl) + z * Math.cos(oblecl);

		// now get the topocentric
		var spherical =  this._rectangularToSpherical(gx, gy, gz);
		spherical.ra = this._normalizeAngle(spherical.ra);
		var ha = this._hourAngleDegrees(epoch, utchours, geoLon, spherical.ra);
		// parallax
		var par = body.parallax(spherical.r);
		// latitude adjustment
		var gclat = geoLat - 0.1924 * Math.sin((2.0*geoLat) * (Math.PI / 180.0));
		gclat = gclat * (Math.PI / 180.0);
		// distance from the center of the earth
		var rho = 0.99833 + 0.00167 * Math.cos((2.0*geoLat) * (Math.PI / 180.0));
		// aux angle
		var g = Math.atan( Math.tan(gclat) / Math.cos(ha));
		g = this._normalizeAngle(g);
		// RA adjustment
		var topRA   = par * rho * Math.cos(gclat) * Math.sin(ha) / Math.cos(spherical.decl);
		// DECL adjustment
		var topDecl = par * rho * Math.sin(gclat) * Math.sin(g - spherical.decl) / Math.sin(g);
		// adjust the geocentric spherical coordinates to topocentric
		spherical.ra = spherical.ra - topRA;
		spherical.decl = spherical.decl - topDecl;
		spherical.ra = this._normalizeAngle(spherical.ra);
		// convert to rectnagle.
		// x axis points to the celestial equator in the south
		x = Math.cos(ha) * Math.cos(spherical.decl);
		// y axis points to the horizon in the west
		y = Math.sin(ha) * Math.cos(spherical.decl);
		// z axis points to the north celestial pole
		z = Math.sin(spherical.decl);
		// rotate alon the y axis (pointing east/west) so z points at the zenith
		var tx = x * Math.sin(gclat) - z * Math.cos(gclat);
		var ty = y;
		var tz = x * Math.cos(gclat) + z * Math.sin(gclat);

		var res = {
			'heliocentric' : {'x' : hx, 'y' : hy, 'z' : hz},
			'geocentric'   : {'x' : gx, 'y' : gy, 'z' : gz},
			'topocentric'  : {'x' : tx, 'y' : ty, 'z' : tz},
			'hourangle'    : ha,
		};
		return res;
	},
	_lunarElongationPhaseAngleMagnitude : function(body, c, epoch, utc, geoLat, geoLon) {
		var sun = this.orbitalElementsSun.data(epoch, utc, geoLat, geoLon);
		var s = sun.coords;
		var slon = Math.atan2(s.geocentric.y, s.geocentric.x);
		var mlon = Math.atan2(c.geocentric.y, c.geocentric.x);
		var mlat = Math.asin(c.geocentric.z / Math.sqrt( c.geocentric.x * c.geocentric.x + c.geocentric.y * c.geocentric.y  + c.geocentric.z * c.geocentric.z) );

		var elong = Math.acos( Math.cos(slon - mlon) * Math.cos(mlat) );
		var FV = Math.PI - elong;
		var phase = ( 1.0 + Math.cos(FV) ) / 2.0;
		var res = {
			'elongation' : elong,
			'phase'      : phase,
			'magnitude'  : 0.0
		};
		return res;

	},
	_planetaryElongationPhaseAngleMagnitudeZero : function(body, coords, epoch) {
		var res = {
			'elongation' : 0.0,
			'phase'      : 0.0,
			'magnitude'  : 0.0
		};
		return res;
	},
	_planetaryElongationPhaseAngleMagnitude: function(body, c, epoch) {
		var r = Math.sqrt(c.heliocentric.x * c.heliocentric.x + c.heliocentric.y * c.heliocentric.y + c.heliocentric.z * c.heliocentric.z);
		var R = Math.sqrt(c.geocentric.x   * c.geocentric.x   + c.geocentric.y   * c.geocentric.y   + c.geocentric.z   * c.geocentric.z);
		var s = this._solarDistanceFromEarth(epoch);

		var elong = Math.acos( ( s*s + R*R - r*r ) / (2*s*R) );
		var FV    = Math.acos( ( r*r + R*R - s*s ) / (2*r*R) );
		var mag   = body.magnitude(r, R, FV);
		var phase = ( 1.0 + Math.cos(FV) ) / 2.0;
		var res   = {
			'elongation' : elong,
			'phase'      : phase,
			'magnitude'  : mag
		};
		return res;
	},
	_solarDistanceFromEarth : function(epoch){
		var body = this.orbitalElementsSun;
		var a = body.semiMajorAxis(epoch);
		var e = body.eccentricity(epoch);
		var M = this._normalizeAngle(body.meanAnomaly(epoch));
		var E = this._eccentricAnomaly(M, e);

		var x = a * (Math.cos(E) - e);
		var y = a * (Math.sin(E) * Math.sqrt(1.0 - e*e));
		var r = Math.sqrt(x*x + y*y);
		return r;


	},
	_planetaryPerturbationsMj : function(epoch) {
		return this.orbitalElementsJupiter.meanAnomaly(epoch);
	},
	_planetaryPerturbationsMs : function(epoch) {
		return this.orbitalElementsSaturn.meanAnomaly(epoch);
	},
	_planetaryPerturbationsMu : function(epoch) {
		return this.orbitalElementsUranus.meanAnomaly(epoch);
	},
};



