
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
	_gmst0Hours : function(epoch){
		var M = 356.0470 + 0.9856002585 * epoch;
		var w = 282.9404 + 4.70935E-5 * epoch;
		var L = this._normalizeAngle(M + w);
		return L / 15.0 + 12.0;
	},
	_sidtimeHours : function(coordinates) {
		var ut = coordinates.date.getUTCHours() + (coordinates.date.getUTCMinutes() / 60.0) + (coordinates.date.getUTCSeconds() / 3600.0);
		var st = this._gmst0Hours(coordinates.epoch) + ut + coordinates.location.longitude / 15.0;
		if(st < 0) st = st + 24;
		if(st > 24) st = st - 24;
		return st;
	},
	_hourAngleDegrees : function (coordinates, longitude) {
		var ha = (this._sidtimeHours(coordinates) * 15.0) - longitude;
		ha = this._normalizeAngle(ha);
		return ha;
	},
	_resultValue : function(params) {
		var res = {
			elements : ((params && params.elements) ? params.elements : {N : 0, i : 0, w : 0, a : 0, e : 0, M : 0}),
			coordinates : {spherical : ((params && params.coordinates && params.coordinates.spherical) ? params.coordinates.spherical : {r : 0, latitude : 0, longitude : 0}),
				       rectangle : ((params && params.coordinates && params.coordinates.rectangle) ? params.coordinates.rectangle : {x : 0, y : 0, z : 0}),
				       horizontal : ((params && params.coordinates && params.coordinates.horizontal) ? params.coordinates.horizontal : {azimuth : 0, altitude : 0})},
			epoch : ((params && params.epoch) ? params.epoch : 0),
			date  : ((params && params.date)  ? params.date  : 0),
			location : ((params && params.location) ? params.location : {latitude : 0, longitude : 0}),
			id : ((params && params.id) ? params.id : "unknown"),
			hourangle : ((params && params.hourangle) ? params.hourangle : 0.0)

		};
		return res;
	},
	degreesToDegreesMinutesSeconds : function(degrees) {
		var b = degrees;

		var d = this._truncate(degrees);

		var b1 = (b - d) * 60.0;

		var m = this._truncate(b1);

		var s = (b1 - m) * 60;

		var res = {degrees : d, minutes : m, seconds : s};

		return res;
	},
	degreesToHoursMinutesSeconds : function(degrees) {
		var b = degrees / 15.0;
		var h = this._truncate(b);

		var b1 = (b - h) * 60.0;
		var m = this._truncate(b1);

		var s = (b1 - m) * 60.0;

		var res = {hours : h, minutes : m, seconds : s};

		return res;
	},
	epochFromDate : function(date) {
		var Y = date.getFullYear();
		var M = date.getMonth() + 1;
		var D = date.getDate();
		var h = date.getHours();
		var m = date.getMinutes();
		var s = date.getSeconds();
		var d = 367 * Y - this._truncate((7 * (Y + this._truncate(((M + 9) / 12)))) / 4) + this._truncate((275 * M) / 9) + D - 730530;
		d = d + ((h + (m / 60.0) + (s / 3600.0)) / 24.0);
		//return -3543;
		return d;
	},
	orbitalElementsMercury : {
		id : 'Mercury',
		longitudeOfAscendingNode : function(epoch) {
			return 48.3313 + 3.24587E-5 * epoch;
		},
		inclination              : function(epoch) {
			return 7.0047 + 5.00E-8 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 29.1241 + 1.01444E-5 * epoch;
		},
		semiMajorAxis            : function(epoch) {
			return 0.387098;
		},
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
		id : 'Venus',
		longitudeOfAscendingNode : function(epoch) {
			return 76.6799 + 2.46590E-5 * epoch;
		},
		inclination              : function(epoch) {
			return 3.3946 + 2.75E-8 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 54.8910 + 1.38374E-5 * epoch;
		},
		semiMajorAxis            : function(epoch) {
			return 0.723330;
		},
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
		id : 'Sun',
		longitudeOfAscendingNode : function(epoch) {
			return 0.0;
		},
		inclination              : function(epoch) {
			return 0.0;
		},
		argumentOfPerihelion     : function(epoch) {
			return 282.9404 + 4.70935E-5 * epoch;
		},
		semiMajorAxis            : function(epoch){
			return 1.0;
		},
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
		id : 'Mars',
		longitudeOfAscendingNode : function(epoch) {
			return 49.5574 + 2.11081E-5 * epoch;
		},
		inclination              : function(epoch) {
			return 1.8497 - 1.78E-8 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 286.5016 + 2.92961E-5 * epoch;
		},
		semiMajorAxis            : function(epoch) {
			return 1.523688;
		},
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
		id : 'Jupiter',
		longitudeOfAscendingNode : function(epoch) {
			return 100.4542 + 2.76854E-5 * epoch;
		},
		inclination              : function(epoch) {
			return 1.3030 - 1.557E-7 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 273.8777 + 1.64505E-5 * epoch;
		},
		semiMajorAxis            : function(epoch){
			return 5.20256;
		},
		eccentricity             : function(epoch) {
			return 0.048498 + 4.469E-9 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 19.8950 + 0.0830853001 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			var self = solarSystemOrbits;
			var res = 0.0;
			var Mj = self._perturbationsMj(epoch);
			var Ms = self._perturbationsMs(epoch);
			res = res + (-0.332 * self._sind(2.0 * Mj - 5.0 * Ms - 67.6));
			res = res + (-0.056 * self._sind(2.0 * Mj - 2.0 * Ms + 21.0));
			res = res + (0.0420 * self._sind(3.0 * Mj - 5.0 * Ms + 21.0));
			res = res + (-0.036 * self._sind(Mj - 2.0 * Ms));
			res = res + (0.0220 * self._cosd(Mj - Ms));
			res = res + (0.0230 * self._sind(2.0 * Mj - 3.0 * Ms + 52));
			res = res + (-0.016 * self._sind(Mj - 5.0 * Ms - 69));
			return res;
		},
		latitudePerturbations    : function(epoch) {
			return 0.0;
		}
	},
	orbitalElementsSaturn : {
		id : 'Saturn',
		longitudeOfAscendingNode : function(epoch) {
			return 113.6634 + 2.38980E-5 * epoch;
		},
		inclination              : function(epoch) {
			return 2.4886 - 1.081E-7 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 339.3939 + 2.97661E-5 * epoch;
		},
		semiMajorAxis            : function(epoch) {
			return 9.55475;
		},
		eccentricity             : function(epoch) {
			return 0.055546 - 9.499E-9 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 316.9670 + 0.0334442282 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			var self = solarSystemOrbits;
			var res = 0.0;
			var Mj = self._perturbationsMj(epoch);
			var Ms = self._perturbationsMs(epoch);
			res = res + ( 0.812 * self._sind(2*Mj - 5*Ms - 67.6));
			res = res + (-0.229 * self._cosd(2*Mj - 4*Ms - 2));
			res = res + ( 0.119 * self._sind(Mj - 2*Ms - 3));
			res = res + ( 0.046 * self._sind(2*Mj - 6*Ms - 69));
			res = res + ( 0.014 * self._sind(Mj - 3*Ms + 32));
			return res;
		},
		latitudePerturbations    : function(epoch) {
			var self = solarSystemOrbits;
			var res = 0.0;
			var Mj = self._perturbationsMj(epoch);
			var Ms = self._perturbationsMs(epoch);
			res = res + (-0.020 * self._cosd(2*Mj - 4*Ms - 2));
			res = res + ( 0.018 * self._sind(2*Mj - 6*Ms - 49));
			return res;
		}
	},
	orbitalElementsUranus : {
		id : 'Uranus',
		longitudeOfAscendingNode : function(epoch) {
			return 74.0005 + 1.3978E-5 * epoch;
		},
		inclination              : function(epoch) {
			return 0.7733 + 1.9E-8 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 96.6612 + 3.0565E-5 * epoch;
		},
		semiMajorAxis            : function(epoch){
			return 19.18171 - 1.55E-8 * epoch;
		},
		eccentricity             : function(epoch) {
			return 0.047318 + 7.45E-9 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 142.5905 + 0.011725806 * epoch;
		},
		longitudePerturbations   : function(epoch) {
			var self = solarSystemOrbits;
			var res = 0.0;
			var Mj = self._perturbationsMj(epoch);
			var Ms = self._perturbationsMs(epoch);
			var Mu = self._perturbationsMu(epoch);
			res = res + ( 0.040 * self._sind(Ms - 2*Mu + 6));
			res = res + ( 0.035 * self._sind(Ms - 3*Mu + 33));
			res = res + (-0.015 * self._sind(Mj - Mu + 20));
			return res;
		},
		latitudePerturbations    : function(epoch) {
			return 0.0;
		}
	},
	orbitalElementsNeptune : {
		id : 'Neptune',
		longitudeOfAscendingNode : function(epoch) {
			return 131.7806 + 3.0173E-5 * epoch;
		},
		inclination              : function(epoch) {
			return 1.7700 - 2.55E-7 * epoch;
		},
		argumentOfPerihelion     : function(epoch) {
			return 272.8461 - 6.027E-6 * epoch;
		},
		semiMajorAxis            : function(epoch){
			return 30.05826 + 3.313E-8 * epoch;
		},
		eccentricity             : function(epoch) {
			return 0.008606 + 2.15E-9 * epoch;
		},
		meanAnomaly              : function(epoch) {
			return 260.2471 + 0.005995147 * epoch;
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
		var a = body.semiMajorAxis(epoch);
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

		var lon = this._normalizeAngle(this._atan2d(y, x));
		var lat = this._asind(z / r);
		// it would be nice to apply perturbations before converting to spherical
		lon = lon + body.longitudePerturbations(epoch);
		lat = lat + body.latitudePerturbations(epoch);

		x = r * this._cosd(lon) * this._cosd(lat);
		y = r * this._sind(lon) * this._cosd(lat);
		z = r * this._sind(lat);


		return this._resultValue({
		                                 elements : { N : N, i : i, w : w, a : a, e : e, M : M},
		                                 coordinates : {
		                                         spherical : {r : r, latitude : lat, longitude : lon},
		                                         rectangle : {x : x, y : y, z : z}
						 },
		                                 epoch : epoch,
		                                 date : date,
		                                 id : body.id
					 });
	},
	coordinatesForBodyAtTimeAndPlace : function(body, date, geoLat, geoLon) {
		var c = this.orbitalCoordinatesForBodyOnDate(date, body);
		var sx = c.coordinates.rectangle.x;
		var sy = c.coordinates.rectangle.y;
		var sz = c.coordinates.rectangle.z;
		c.location  = {latitude : geoLat, longitude : geoLon};
		if(body.id != 'Sun') {
			var s = this.orbitalCoordinatesForBodyOnDate(date, this.orbitalElementsSun);
			sx = s.coordinates.rectangle.x + sx;
			sy = s.coordinates.rectangle.y + sy;
			sz = s.coordinates.rectangle.z + sz;
		}
		var oblecl = 23.4393 - 3.563E-7 * c.epoch;
		var xx = sx;
		var yy = sy;
		var zz = sz;
		var x = xx;
		var y = yy * this._cosd(oblecl) - zz * this._sind(oblecl);
		var z = yy * this._sind(oblecl) + zz * this._cosd(oblecl);
		var r    = Math.sqrt(x*x + y*y + z*z);
		var ra   = this._normalizeAngle(this._atan2d(y, x));
		var decl = this._asind(z / r);

		c.hourangle = this._hourAngleDegrees(c, ra);
		c.coordinates.rectangle.x = x;
		c.coordinates.rectangle.y = y;
		c.coordinates.rectangle.z = z;
		c.coordinates.spherical.r = r;
		c.coordinates.spherical.longitude = ra;
		c.coordinates.spherical.latitude  = decl;

		var ha    = c.hourangle;
		var ppar  = (8.794/3600) / r;
		var gclat = geoLat - 0.1924 * this._sind(2.0 * geoLat);
		var rho   = 0.99833 + 0.00167 * this._cosd(2.0 * geoLat);
		var g     = this._normalizeAngle((this._atand( this._tand(gclat) / this._cosd(ha))));

		var topRA   = ra  - ppar * rho * this._cosd(gclat) * this._sind(ha) / this._cosd(decl);
		var topDecl = decl - ppar * rho * this._sind(gclat) * this._sind(g - decl) / this._sind(g);

		ra   = this._normalizeAngle(topRA);
		decl = topDecl;

		xx = this._cosd(ha) * this._cosd(decl);
		yy = this._sind(ha) * this._cosd(decl);
		zz = this._sind(decl);

		x = xx * this._sind(geoLat) - zz * this._cosd(geoLat);
		y = yy;
		z = xx * this._cosd(geoLat) + zz * this._sind(geoLat);

		var azm = this._atan2d(y, x) + 180.0;
		var alt = this._atan2d(z, Math.sqrt(x*x + y*y));

		azm = this._normalizeAngle(azm);
		c.coordinates.horizontal = {azimuth : azm, altitude : alt};

		return this._resultValue(c);
	},
	_perturbationsMj : function(epoch) {
		return this.orbitalElementsJupiter.meanAnomaly(epoch);
	},
	_perturbationsMs : function(epoch) {
		return this.orbitalElementsSaturn.meanAnomaly(epoch);
	},
	_perturbationsMu : function(epoch) {
		return this.orbitalElementsUranus.meanAnomaly(epoch);
	},
};



