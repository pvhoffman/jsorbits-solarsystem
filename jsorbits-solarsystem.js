
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
                M = 6.214192441848251 + 0.01720196961933223 * epoch;
                w = 4.938241566909764 + 8.219366312879496e-07 * epoch;
                L = this._normalizeAngle((M + w));
                L = L / ((Math.PI * 2.0) / 24.0);
                L = L + 12.0;
                return L;
	},
	_sidtimeHours : function(epoch, utchours, longitude){//coordinates) {
                var g = this._gmst0Hours(epoch);
                var r = g + utchours + longitude / 15.0;
                if (r < 0.0)  r = r + 24.0;
                if (r > 24.0) r = r - 24.0;
                return r;
	},
	_hourAngleDegrees : function (epoch, utchours, longitude, ra){ //coordinates, longitude) {
                s = this._sidtimeHours(epoch, utchours, longitude);
                r = (s * ((Math.PI * 2.0) / 24.0)) - ra;
                r = this._normalizeAngle(r);
                return r;
	},
        _rectangularToSpherical : function(x, y, z){
            var e = 0.0001;
            var r = Math.sqrt(x*x + y*y + z*z);
            var ra = 0.0; 
            var decl = 0.0;
            if (Math.abs(x) < e &&  Math.abs(y) < e){
                decl = Math.atan2(z, math.sqrt(x*x + y*y));
            } else {
                ra = Math.atan2(y, x);
                decl = Math.asin(z / r);
            }
            return {'r' : r, 'ra' : ra, 'decl' : decl};
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
                coordinates : function(date, geoLat, geoLon) {
                        var self = solarSystemOrbits;
                        var epoch = self.epochFromDate(date);
                        var utchours = date.getUTCHours() + date.getUTCMinutes() / 60.0 + date.getUTCSeconds() / 3600.0;
                        return self._planetaryCoordinates(epoch, this, geoLat, geoLon, utchours);
                }
	},
	orbitalElementsVenus : {
		id : 'Venus',
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
                coordinates : function(date, geoLat, geoLon) {
                        var self = solarSystemOrbits;
                        var epoch = self.epochFromDate(date);
                        var utchours = date.getUTCHours() + date.getUTCMinutes() / 60.0 + date.getUTCSeconds() / 3600.0;
                        return self._planetaryCoordinates(epoch, this, geoLat, geoLon, utchours);
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
                coordinates : function(date, geoLat, geoLon) {
                        var self = solarSystemOrbits;
                        var epoch = self.epochFromDate(date);
                        var utchours = date.getUTCHours() + date.getUTCMinutes() / 60.0 + date.getUTCSeconds() / 3600.0;
                        return self._planetaryCoordinates(epoch, this, geoLat, geoLon, utchours);
                }
        },
	orbitalElementsMars : {
		id : 'Mars',
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
                coordinates : function(date, geoLat, geoLon) {
                        var self = solarSystemOrbits;
                        var epoch = self.epochFromDate(date);
                        var utchours = date.getUTCHours() + date.getUTCMinutes() / 60.0 + date.getUTCSeconds() / 3600.0;
                        return self._planetaryCoordinates(epoch, this, geoLat, geoLon, utchours);
                }
	},
	orbitalElementsJupiter : {
		id : 'Jupiter',
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
                coordinates : function(date, geoLat, geoLon) {
                        var self = solarSystemOrbits;
                        var epoch = self.epochFromDate(date);
                        var utchours = date.getUTCHours() + date.getUTCMinutes() / 60.0 + date.getUTCSeconds() / 3600.0;
                        return self._planetaryCoordinates(epoch, this, geoLat, geoLon, utchours);
                }
	},
	orbitalElementsSaturn : {
		id : 'Saturn',
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
                coordinates : function(date, geoLat, geoLon) {
                        var self = solarSystemOrbits;
                        var epoch = self.epochFromDate(date);
                        var utchours = date.getUTCHours() + date.getUTCMinutes() / 60.0 + date.getUTCSeconds() / 3600.0;
                        return self._planetaryCoordinates(epoch, this, geoLat, geoLon, utchours);
                }

	},
	orbitalElementsUranus : {
		id : 'Uranus',
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
                coordinates : function(date, geoLat, geoLon) {
                        var self = solarSystemOrbits;
                        var epoch = self.epochFromDate(date);
                        var utchours = date.getUTCHours() + date.getUTCMinutes() / 60.0 + date.getUTCSeconds() / 3600.0;
                        return self._planetaryCoordinates(epoch, this, geoLat, geoLon, utchours);
                }
	},
	orbitalElementsNeptune : {
		id : 'Neptune',
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
                coordinates : function(date, geoLat, geoLon) {
                        var self = solarSystemOrbits;
                        var epoch = self.epochFromDate(date);
                        var utchours = date.getUTCHours() + date.getUTCMinutes() / 60.0 + date.getUTCSeconds() / 3600.0;
                        return self._planetaryCoordinates(epoch, this, geoLat, geoLon, utchours);
                }
	},
	orbitalElementsEarthMoon : {
		id : 'Earth Moon',
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
			return 0.0;
		},
		latitudePerturbations    : function(epoch) {
			return 0.0;
		},
	        distancePerturbations : function(epoch) {
                        return 0.0;
                },
                parallax : function(r) {
                        var par = Math.asin( 1.0 / r );
                        return par;
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

		var lon = this._normalizeAngle(Math.atan2(y, x));
		var lat = Math.asin(z / r);

		lon = lon + perturbationsLon;
		lat = lat + perturbationsLat;
                r   = r + perturbationsDist;

		var hx = r * Math.cos(lon) * Math.cos(lat);
		var hy = r * Math.sin(lon) * Math.cos(lat);
		var hz = r * Math.sin(lat);

                var x = hx;
                var y = hy;
                var z = hz;

                // now get the geocentric
                if(body.id != 'Sun'){
                        var sun = this._planetaryCoordinates(epoch, this.orbitalElementsSun, geoLat, geoLon, utchours); 
                        x = sun.heliocentric.x + x;
                        y = sun.heliocentric.y + y;
                        z = sun.heliocentric.z + z;
                }
                var oblecl = 0.4090929593627069 - 6.218608124855796e-09 * epoch;

                var gx = x;
                var gy = y * Math.cos(oblecl) - z * Math.sin(oblecl);
                var gz = y * Math.sin(oblecl) + z * Math.cos(oblecl);

                // now get the topocentric
                var spherical =  this._rectangularToSpherical(gx, gy, gz);
                spherical.ra = this._normalizeAngle(spherical.ra);

                var ha = this._hourAngleDegrees(epoch, utchours, geoLon, spherical.ra); 
                ha = this._normalizeAngle(ha);
                // parallax
                var par = body.parallax(spherical.r);

                var gclat = geoLat - 0.1924 * Math.sin((2.0*geoLat) * (Math.PI / 180.0));
                gclat = gclat * (Math.PI / 180.0);

                var rho = 0.99833 + 0.00167 * Math.cos((2.0*geoLat) * (Math.PI / 180.0));

                var g = Math.atan( Math.tan(gclat) / Math.cos(ha));
                g = this._normalizeAngle(g);

                var topRA   = par * rho * Math.cos(gclat) * Math.sin(ha) / Math.cos(spherical.decl);
                var topDecl = par * rho * Math.sin(gclat) * Math.sin(g - spherical.decl) / Math.sin(g);

                spherical.ra = spherical.ra - topRA;
                spherical.decl = spherical.decl - topDecl;
                spherical.ra = this._normalizeAngle(spherical.ra)

                ha = this._hourAngleDegrees(epoch, utchours, geoLon, spherical.ra);
                ha = this._normalizeAngle(ha);

                x = Math.cos(ha) * Math.cos(spherical.decl);
                y = Math.sin(ha) * Math.cos(spherical.decl);
                z = Math.sin(spherical.decl);

                var tx = x * Math.sin(gclat) - z * Math.cos(gclat);
        	var ty = y;
	        var tz = x * Math.cos(gclat) + z * Math.sin(gclat);

                var res = {
                        'heliocentric' : {'x' : hx, 'y' : hy, 'z' : hz},
                        'geocentric'   : {'x' : gx, 'y' : gy, 'z' : gz},
                        'topocentric'  : {'x' : tx, 'y' : ty, 'z' : tz},
                        'hourangle'    : ha
                };
                return res;
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



