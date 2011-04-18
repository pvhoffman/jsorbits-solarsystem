var solarSystemOrbits = {
	_eccentricAnomaly : function(M, e) {
		var E0 = M + (180.0/Math.PI) * e * this._sind(M) * (1.0 + e * this._cosd(M));
		var E1 = E0 - (E0 - (180.0/Math.PI) * e * this._sind(E0) - M) / (1.0 - e * this._cosd(E0));
		var D;
		var I = 0;

		do {
			E0 = E1;
			E1 = E0 - (E0 - (180.0/Math.PI) * e * this._sind(E0) - M) / (1 - e * this._cosd(E0));
			D = Math.abs((E1 - E0));
			I = I + 1;
		} while(D > 0.0005 && I < 10);

		if(I == 10 && D > 0.0005) {
			alert("Integration failed, life sucks...");
		}

		return E1;


	},
	_truncate : function(d) {
		if(d < 0) return Math.ceil(d);
		return Math.floor(d);
	},
	_normalizeAngle : function(a){
		return a - Math.floor(a / 360.0) * 360.0;
	},

	_sind : function(x) {
		return Math.sin(x * (Math.PI / 180.0));
	},
	_cosd : function(x) {
		return Math.cos(x * (Math.PI / 180.0));
	},
	_tand : function(x) {
		return Math.tan(x * (Math.PI / 180.0));
	},
	_asind : function(x) {
		return (180.0 / Math.PI) * Math.asin(x);
	},
	_acosd : function(x) {
		return (180.0 / Math.PI) * Math.acos(x);
	},
	_atand : function(x) {
		return (180.0 / Math.PI) * Math.atan(x);
	},
	_atan2d : function(y, x) {
		return (180.0 / Math.PI) * Math.atan2(y, x);
	},
	epochFromDate : function(date) {
		var Y = date.getFullYear();
		var M = date.getMonth() + 1;
		var D = date.getDate();
		var h = date.getHours();
		var m = date.getMinutes();
		var s = date.getSeconds();
/*
                var Y = date.getUTCFullYear();
                var M = date.getUTCMonth() + 1;
                var D = date.getUTCDate();
                var h = date.getUTCHours();
                var m = date.getUTCMinutes();
                var s = date.getUTCSeconds();
 */
		var d = 367 * Y - this._truncate((7 * (Y + this._truncate(((M + 9) / 12)))) / 4) + this._truncate((275 * M) / 9) + D - 730530;
		d = d + ((h + (m / 60.0) + (s / 3600.0)) / 24.0);
		//return -3543;
		return d;


	},
	orbitalElementsMercury : {
		longitudeOfAscendingNode : function(epoch) {
			return 48.3313 + 3.24587E-5 * epoch;
		},
		inclination              : function(epoch) {
			return 7.0047 + 5.00E-8 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 29.1241 + 1.01444E-5 * epoch;
		},
		semiMajorAxis            : 0.387098,
		eccentricity             : function(epoch) {
			return 0.205635 + 5.59E-10 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 168.6562 + 4.0923344368 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			return 0.0;
		},
		latitudePerturbations    : function(epoch) {
			return 0.0;
		}
	},
	orbitalElementsVenus : {
		longitudeOfAscendingNode : function(epoch) {
			return 76.6799 + 2.46590E-5 * epoch;
		},
		inclination              : function(epoch) {
			return 3.3946 + 2.75E-8 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 54.8910 + 1.38374E-5 * epoch;
		},
		semiMajorAxis            :0.723330,
		eccentricity             : function(epoch) {
			return 0.006773 - 1.302E-9 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 48.0052 + 1.6021302244 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			return 0.0;
		},
		latitudePerturbations    : function(epoch) {
			return 0.0;
		}
	},
	orbitalElementsSun : {
		longitudeOfAscendingNode : function(epoch) {
			return 0.0;
		},
		inclination              : function(epoch) {
			return 0.0;
		},
		argumentOfPerihelion     : function(epoch) {
			return 282.9404 + 4.70935E-5 * epoch;
		},
		semiMajorAxis            : 1.0,
		eccentricity             : function(epoch) {
			return 0.016709 - 1.151E-9 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 356.0470 + 0.9856002585 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			return 0.0;
		},
		latitudePerturbations    : function(epoch) {
			return 0.0;
		}
	},
	orbitalElementsMars : {
		longitudeOfAscendingNode : function(epoch) {
			return 49.5574 + 2.11081E-5 * epoch;
		},
		inclination              : function(epoch) {
			return 1.8497 - 1.78E-8 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 286.5016 + 2.92961E-5 * epoch;
		},
		semiMajorAxis            : 1.523688,
		eccentricity             : function(epoch) {
			return 0.093405     + 2.516E-9 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 18.6021 + 0.5240207766 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			return 0.0;
		},
		latitudePerturbations    : function(epoch) {
			return 0.0;
		}
	},
	orbitalElementsJupiter : {
		longitudeOfAscendingNode : function(epoch) {
			return 100.4542 + 2.76854E-5 * epoch;
		},
		inclination              : function(epoch) {
			return 1.3030 - 1.557E-7 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 273.8777 + 1.64505E-5 * epoch;
		},
		semiMajorAxis            : 5.20256,
		eccentricity             : function(epoch) {
			return 0.048498 + 4.469E-9 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 19.8950 + 0.0830853001 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			return 0.0;
		},
		latitudePerturbations    : function(epoch) {
			return 0.0;
		}
	},
	orbitalCoordinatesForBodyOnDate: function(date, body) {
		var epoch = this.epochFromDate(date);

		var N = this._normalizeAngle(body.longitudeOfAscendingNode(epoch));
		var i = this._normalizeAngle(body.inclination(epoch));
		var w = this._normalizeAngle(body.argumentOfPerihelion(epoch));
		var a = body.semiMajorAxis;
		var e = body.eccentricity(epoch);
		var M = this._normalizeAngle(body.meanAnomaly(epoch));

		var E = this._eccentricAnomaly(M, e);

		var x = a * (this._cosd(E) - e);
		var y = a * (this._sind(E) * Math.sqrt(1.0 - e*e));
		var z = 0;

		var r = Math.sqrt(x*x + y*y);
		var v = this._atan2d(y, x);

		x = r * (this._cosd(N) * this._cosd(v+w) - this._sind(N) * this._sind(v+w) * this._cosd(i));
		y = r * (this._sind(N) * this._cosd(v+w) + this._cosd(N) * this._sind(v+w) * this._cosd(i));
		z = r * this._sind(v+w) * this._sind(i);

		// convert to spherical coordinates
		r = Math.sqrt(x*x + y*y + z*z);

		var lon = this._normalizeAngle(this._atan2d(y, x));
		var lat = this._asind(z / r);
		var res = {elements : { N : N, i : i, w : w, a : a, e : e, M : M}, spherical : {r : r, lat : lat, lon : lon}, rectangle : {x : x, y : y, z : z}, epoch : epoch, date : date};
		return res;
	},

	convertToGeocentricCoordinates : function(coordinates){
		var xg;
		var yg;
		var zg;
		var oblecl = 23.4393 - 3.563E-7 * coordinates.epoch;

		if(coordinates.rectangle.z == 0) {
			xg = coordinates.rectangle.x;
			yg = coordinates.rectangle.y;
			zg = coordinates.rectangle.z;
		} else {
			var x = this.orbitalCoordinatesForBodyOnDate(coordinates.date, this.orbitalElementsSun);
			xg = x.rectangle.x + coordinates.rectangle.x;
			yg = x.rectangle.y + coordinates.rectangle.y;
			zg = x.rectangle.z + coordinates.rectangle.z;
		}

		var xe = xg;
		var ye = yg * this._cosd(oblecl) - zg * this._sind(oblecl);
		var ze = yg * this._sind(oblecl) + zg * this._cosd(oblecl);

		r = Math.sqrt(xe*xe + ye*ye + ze*ze);
		var lon = this._normalizeAngle(this._atan2d(ye, xe));
		var lat = this._asind(ze / r);
		var res = {elements :
			   { N : coordinates.elements.N
			     , i : coordinates.elements.i
			     , w : coordinates.elements.w
			     , a : coordinates.elements.a
			     , e : coordinates.elements.e
			     , M : coordinates.elements.M}
			   , spherical :
			   {r : r
			    , lat : lat
			    , lon : lon}, rectangle :
			   {x : xe
			    , y : ye
			    , z : ze}
			   , epoch : coordinates.epoch
			   , date : coordinates.date};
		return res;
	},
	siderealTimeHourAngleAltitudeAndAzimuth : function(coordinates, geoLat, geoLon) {
		var L = coordinates.elements.M + coordinates.elements.w;
		var gmst0 = this._normalizeAngle(L) / 15 + 12;
		var ut    = coordinates.date.getUTCHours() + (coordinates.date.getUTCMinutes() / 60.0);
		var sidtime = gmst0 + ut + geoLon / 15;

		if(sidtime < 0) sidtime = sidtime + 24;
		if(sidtime > 24) sidtime = sidtime - 24;

		var ha = (sidtime * 15) - coordinates.spherical.lon;

		var x = this._cosd(ha) * this._cosd(coordinates.spherical.lat);
		var y = this._sind(ha) * this._cosd(coordinates.spherical.lat);
		var z = this._sind(coordinates.spherical.lat);

		var xh = x * this._sind(geoLat) - z * this._cosd(geoLat);
		var yh = y;
		var zh = x * this._cosd(geoLat) + z * this._sind(geoLat);

		var azimuth  = this._atan2d(yh, xh) + 180;
		var altitude = this._atan2d(zh, Math.sqrt(xh*xh + yh*yh));

		var res = {elements :
			   { N : coordinates.elements.N
			     , i : coordinates.elements.i
			     , w : coordinates.elements.w
			     , a : coordinates.elements.a
			     , e : coordinates.elements.e
			     , M : coordinates.elements.M}
			   , gmst0 : gmst0
			   , ut : ut
			   , ha : ha
			   , sidtime : sidtime
			   , altitude : altitude
			   , azimuth : azimuth
			   , epoch : coordinates.epoch
			   , date : coordinates.date};
		return res;
	},

	_perturbationsMj : function(epoch) {
		return 0;      //this.orbitalElementsJupiter.meanAnomaly(epoch);
	},
	_perturbationsMs : function(epoch) {
		return 0;      //this.orbitalElementsSaturn.meanAnomaly(epoch);
	},
	_perturbationsMu : function(epoch) {
		return 0;      //this.orbitalElementsUranus.meanAnomaly(epoch);
	},
};
