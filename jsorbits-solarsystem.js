
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
	_meanLongitude : function(coordinates){
		return this._normalizeAngle((coordinates.elements.M + coordinates.elements.w));
	},
	//_gmst0Hours : function(coordinates){
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
				       topocentric : ((params && params.coordinates && params.coordinates.topocentric) ? params.coordinates.topocentric : {ra : 0, decl : 0}),
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

		var lon = this._normalizeAngle(this._atan2d(y, x));
		var lat = this._asind(z / r);
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

	convertToGeocentricCoordinates : function(coordinates){
		var xg;
		var yg;
		var zg;
		var oblecl = 23.4393 - 3.563E-7 * coordinates.epoch;

		if(coordinates.coordinates.rectangle.z == 0) {
			xg = coordinates.coordinates.rectangle.x;
			yg = coordinates.coordinates.rectangle.y;
			zg = coordinates.coordinates.rectangle.z;
		} else {
			var x = this.orbitalCoordinatesForBodyOnDate(coordinates.date, this.orbitalElementsSun);
			xg = x.coordinates.rectangle.x + coordinates.coordinates.rectangle.x;
			yg = x.coordinates.rectangle.y + coordinates.coordinates.rectangle.y;
			zg = x.coordinates.rectangle.z + coordinates.coordinates.rectangle.z;
		}

		var xe = xg;
		var ye = yg * this._cosd(oblecl) - zg * this._sind(oblecl);
		var ze = yg * this._sind(oblecl) + zg * this._cosd(oblecl);

		var r = Math.sqrt(xe*xe + ye*ye + ze*ze);

		var lon = this._normalizeAngle(this._atan2d(ye, xe));
		var lat = this._asind(ze / r);
		return this._resultValue({
		                                 elements : { N : coordinates.elements.N,
		                                              i : coordinates.elements.i,
		                                              w : coordinates.elements.w,
		                                              a : coordinates.elements.a,
		                                              e : coordinates.elements.e,
		                                              M : coordinates.elements.M},
		                                 coordinates : {
		                                         spherical : {r : r, latitude : lat, longitude : lon},
		                                         rectangle : {x : xe, y : ye, z : ze}
						 },
		                                 epoch : coordinates.epoch,
		                                 date : coordinates.date,
		                                 id   : coordinates.id || "unknown geocentric"
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
		var ppar = (8.794/3600.0) / r
		           var gclat = geoLat - 0.1924 * this._sind(2.0 * geoLat);
		var rho   = 0.99833 + 0.00167 * this._cosd(2.0 * geoLat);
		var g     = this._atand( this._tand(gclat) / this._cosd(c.hourangle) );
		var tdecl = decl - ppar * rho * this._sind(gclat) * this._sind(g - decl) / this._sind(g);
		c.coordinates.rectangle.x = x;
		c.coordinates.rectangle.y = y;
		c.coordinates.rectangle.z = z;
		c.coordinates.spherical.r = r;
		c.coordinates.spherical.longitude = ra;
		c.coordinates.spherical.latitude  = decl;
		decl = tdecl;
		xx = this._cosd(c.hourangle) * this._cosd(decl);
		yy = this._sind(c.hourangle) * this._cosd(decl);
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
		return 0;      //this.orbitalElementsJupiter.meanAnomaly(epoch);
	},
	_perturbationsMs : function(epoch) {
		return 0;      //this.orbitalElementsSaturn.meanAnomaly(epoch);
	},
	_perturbationsMu : function(epoch) {
		return 0;      //this.orbitalElementsUranus.meanAnomaly(epoch);
	},
};



