var solarSystemOrbits = {
	_eccentricAnomaly : function(M, e) {
		var E0 = M + (180.0/Math.PI) * e * this._sind(M) * (1.0 + e * this._cosd(M));
		var E1 = E0 - (E0 - (180.0/Math.PI) * e * this._sind(E0) - M) / (1.0 - e * this._cosd(E0));
		var D;

		do {
			E0 = E1;
			E1 = E0 - (E0 - (180.0/Math.PI) * e * this._sind(E0) - M) / (1 - e * this._cosd(E0));
			D = Math.abs((E1 - E0));
		} while(D > 0.0005);

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
		var Y = date.getUTCFullYear();
		var M = date.getUTCMonth() + 1;
		var D = date.getUTCDate();
		var h = date.getUTCHours();
		var m = date.getUTCMinutes();
		var s = date.getUTCSeconds();
		var d = 367 * Y - this._truncate((7 * (Y + this._truncate(((M + 9) / 12)))) / 4) + this._truncate((275 * M) / 9) + D - 730530;
		d = d + ((h + (m / 60.0) + (s / 3600.0)) / 24.0);
		return d


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

	distanceAndTrueAnomalyForBody : function(date, body) {
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

		var r = Math.sqrt(x*x + y*y);
		var v = this._atan2d(y, x);

		var res = {distance : r, trueAnomaly : v}

		return res;
	}


};
