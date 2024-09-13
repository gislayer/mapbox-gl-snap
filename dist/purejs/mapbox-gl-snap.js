(function (global, factory) {
    typeof exports === 'object' && typeof module !== 'undefined' ? module.exports = factory() :
    typeof define === 'function' && define.amd ? define(factory) :
    (global = typeof globalThis !== 'undefined' ? globalThis : global || self, global.MapboxSnap = factory());
})(this, (function () { 'use strict';

    /******************************************************************************
    Copyright (c) Microsoft Corporation.

    Permission to use, copy, modify, and/or distribute this software for any
    purpose with or without fee is hereby granted.

    THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
    REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
    AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
    INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
    LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
    OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
    PERFORMANCE OF THIS SOFTWARE.
    ***************************************************************************** */
    /* global Reflect, Promise, SuppressedError, Symbol, Iterator */


    var __assign = function() {
        __assign = Object.assign || function __assign(t) {
            for (var s, i = 1, n = arguments.length; i < n; i++) {
                s = arguments[i];
                for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p)) t[p] = s[p];
            }
            return t;
        };
        return __assign.apply(this, arguments);
    };

    typeof SuppressedError === "function" ? SuppressedError : function (error, suppressed, message) {
        var e = new Error(message);
        return e.name = "SuppressedError", e.error = error, e.suppressed = suppressed, e;
    };

    /**
     * @module helpers
     */
    /**
     * Earth Radius used with the Harvesine formula and approximates using a spherical (non-ellipsoid) Earth.
     *
     * @memberof helpers
     * @type {number}
     */
    var earthRadius$2 = 6371008.8;
    /**
     * Unit of measurement factors using a spherical (non-ellipsoid) earth radius.
     *
     * @memberof helpers
     * @type {Object}
     */
    var factors$2 = {
        centimeters: earthRadius$2 * 100,
        centimetres: earthRadius$2 * 100,
        degrees: earthRadius$2 / 111325,
        feet: earthRadius$2 * 3.28084,
        inches: earthRadius$2 * 39.37,
        kilometers: earthRadius$2 / 1000,
        kilometres: earthRadius$2 / 1000,
        meters: earthRadius$2,
        metres: earthRadius$2,
        miles: earthRadius$2 / 1609.344,
        millimeters: earthRadius$2 * 1000,
        millimetres: earthRadius$2 * 1000,
        nauticalmiles: earthRadius$2 / 1852,
        radians: 1,
        yards: earthRadius$2 * 1.0936,
    };
    /**
     * Wraps a GeoJSON {@link Geometry} in a GeoJSON {@link Feature}.
     *
     * @name feature
     * @param {Geometry} geometry input geometry
     * @param {Object} [properties={}] an Object of key-value pairs to add as properties
     * @param {Object} [options={}] Optional Parameters
     * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
     * @param {string|number} [options.id] Identifier associated with the Feature
     * @returns {Feature} a GeoJSON Feature
     * @example
     * var geometry = {
     *   "type": "Point",
     *   "coordinates": [110, 50]
     * };
     *
     * var feature = turf.feature(geometry);
     *
     * //=feature
     */
    function feature$2(geom, properties, options) {
        if (options === void 0) { options = {}; }
        var feat = { type: "Feature" };
        if (options.id === 0 || options.id) {
            feat.id = options.id;
        }
        if (options.bbox) {
            feat.bbox = options.bbox;
        }
        feat.properties = properties || {};
        feat.geometry = geom;
        return feat;
    }
    /**
     * Creates a {@link Point} {@link Feature} from a Position.
     *
     * @name point
     * @param {Array<number>} coordinates longitude, latitude position (each in decimal degrees)
     * @param {Object} [properties={}] an Object of key-value pairs to add as properties
     * @param {Object} [options={}] Optional Parameters
     * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
     * @param {string|number} [options.id] Identifier associated with the Feature
     * @returns {Feature<Point>} a Point feature
     * @example
     * var point = turf.point([-75.343, 39.984]);
     *
     * //=point
     */
    function point$2(coordinates, properties, options) {
        if (options === void 0) { options = {}; }
        if (!coordinates) {
            throw new Error("coordinates is required");
        }
        if (!Array.isArray(coordinates)) {
            throw new Error("coordinates must be an Array");
        }
        if (coordinates.length < 2) {
            throw new Error("coordinates must be at least 2 numbers long");
        }
        if (!isNumber$2(coordinates[0]) || !isNumber$2(coordinates[1])) {
            throw new Error("coordinates must contain numbers");
        }
        var geom = {
            type: "Point",
            coordinates: coordinates,
        };
        return feature$2(geom, properties, options);
    }
    /**
     * Creates a {@link LineString} {@link Feature} from an Array of Positions.
     *
     * @name lineString
     * @param {Array<Array<number>>} coordinates an array of Positions
     * @param {Object} [properties={}] an Object of key-value pairs to add as properties
     * @param {Object} [options={}] Optional Parameters
     * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
     * @param {string|number} [options.id] Identifier associated with the Feature
     * @returns {Feature<LineString>} LineString Feature
     * @example
     * var linestring1 = turf.lineString([[-24, 63], [-23, 60], [-25, 65], [-20, 69]], {name: 'line 1'});
     * var linestring2 = turf.lineString([[-14, 43], [-13, 40], [-15, 45], [-10, 49]], {name: 'line 2'});
     *
     * //=linestring1
     * //=linestring2
     */
    function lineString$1(coordinates, properties, options) {
        if (options === void 0) { options = {}; }
        if (coordinates.length < 2) {
            throw new Error("coordinates must be an array of two or more positions");
        }
        var geom = {
            type: "LineString",
            coordinates: coordinates,
        };
        return feature$2(geom, properties, options);
    }
    /**
     * Takes one or more {@link Feature|Features} and creates a {@link FeatureCollection}.
     *
     * @name featureCollection
     * @param {Feature[]} features input features
     * @param {Object} [options={}] Optional Parameters
     * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
     * @param {string|number} [options.id] Identifier associated with the Feature
     * @returns {FeatureCollection} FeatureCollection of Features
     * @example
     * var locationA = turf.point([-75.343, 39.984], {name: 'Location A'});
     * var locationB = turf.point([-75.833, 39.284], {name: 'Location B'});
     * var locationC = turf.point([-75.534, 39.123], {name: 'Location C'});
     *
     * var collection = turf.featureCollection([
     *   locationA,
     *   locationB,
     *   locationC
     * ]);
     *
     * //=collection
     */
    function featureCollection$2(features, options) {
        if (options === void 0) { options = {}; }
        var fc = { type: "FeatureCollection" };
        if (options.id) {
            fc.id = options.id;
        }
        if (options.bbox) {
            fc.bbox = options.bbox;
        }
        fc.features = features;
        return fc;
    }
    /**
     * Creates a {@link Feature<MultiLineString>} based on a
     * coordinate array. Properties can be added optionally.
     *
     * @name multiLineString
     * @param {Array<Array<Array<number>>>} coordinates an array of LineStrings
     * @param {Object} [properties={}] an Object of key-value pairs to add as properties
     * @param {Object} [options={}] Optional Parameters
     * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
     * @param {string|number} [options.id] Identifier associated with the Feature
     * @returns {Feature<MultiLineString>} a MultiLineString feature
     * @throws {Error} if no coordinates are passed
     * @example
     * var multiLine = turf.multiLineString([[[0,0],[10,10]]]);
     *
     * //=multiLine
     */
    function multiLineString$1(coordinates, properties, options) {
        if (options === void 0) { options = {}; }
        var geom = {
            type: "MultiLineString",
            coordinates: coordinates,
        };
        return feature$2(geom, properties, options);
    }
    /**
     * Convert a distance measurement (assuming a spherical Earth) from radians to a more friendly unit.
     * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
     *
     * @name radiansToLength
     * @param {number} radians in radians across the sphere
     * @param {string} [units="kilometers"] can be degrees, radians, miles, inches, yards, metres,
     * meters, kilometres, kilometers.
     * @returns {number} distance
     */
    function radiansToLength$1(radians, units) {
        if (units === void 0) { units = "kilometers"; }
        var factor = factors$2[units];
        if (!factor) {
            throw new Error(units + " units is invalid");
        }
        return radians * factor;
    }
    /**
     * Convert a distance measurement (assuming a spherical Earth) from a real-world unit into radians
     * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
     *
     * @name lengthToRadians
     * @param {number} distance in real units
     * @param {string} [units="kilometers"] can be degrees, radians, miles, inches, yards, metres,
     * meters, kilometres, kilometers.
     * @returns {number} radians
     */
    function lengthToRadians$2(distance, units) {
        if (units === void 0) { units = "kilometers"; }
        var factor = factors$2[units];
        if (!factor) {
            throw new Error(units + " units is invalid");
        }
        return distance / factor;
    }
    /**
     * Converts an angle in radians to degrees
     *
     * @name radiansToDegrees
     * @param {number} radians angle in radians
     * @returns {number} degrees between 0 and 360 degrees
     */
    function radiansToDegrees$2(radians) {
        var degrees = radians % (2 * Math.PI);
        return (degrees * 180) / Math.PI;
    }
    /**
     * Converts an angle in degrees to radians
     *
     * @name degreesToRadians
     * @param {number} degrees angle between 0 and 360 degrees
     * @returns {number} angle in radians
     */
    function degreesToRadians$2(degrees) {
        var radians = degrees % 360;
        return (radians * Math.PI) / 180;
    }
    /**
     * isNumber
     *
     * @param {*} num Number to validate
     * @returns {boolean} true/false
     * @example
     * turf.isNumber(123)
     * //=true
     * turf.isNumber('foo')
     * //=false
     */
    function isNumber$2(num) {
        return !isNaN(num) && num !== null && !Array.isArray(num);
    }

    /**
     * Unwrap a coordinate from a Point Feature, Geometry or a single coordinate.
     *
     * @name getCoord
     * @param {Array<number>|Geometry<Point>|Feature<Point>} coord GeoJSON Point or an Array of numbers
     * @returns {Array<number>} coordinates
     * @example
     * var pt = turf.point([10, 10]);
     *
     * var coord = turf.getCoord(pt);
     * //= [10, 10]
     */
    function getCoord$1(coord) {
        if (!coord) {
            throw new Error("coord is required");
        }
        if (!Array.isArray(coord)) {
            if (coord.type === "Feature" &&
                coord.geometry !== null &&
                coord.geometry.type === "Point") {
                return coord.geometry.coordinates;
            }
            if (coord.type === "Point") {
                return coord.coordinates;
            }
        }
        if (Array.isArray(coord) &&
            coord.length >= 2 &&
            !Array.isArray(coord[0]) &&
            !Array.isArray(coord[1])) {
            return coord;
        }
        throw new Error("coord must be GeoJSON Point or an Array of numbers");
    }
    /**
     * Unwrap coordinates from a Feature, Geometry Object or an Array
     *
     * @name getCoords
     * @param {Array<any>|Geometry|Feature} coords Feature, Geometry Object or an Array
     * @returns {Array<any>} coordinates
     * @example
     * var poly = turf.polygon([[[119.32, -8.7], [119.55, -8.69], [119.51, -8.54], [119.32, -8.7]]]);
     *
     * var coords = turf.getCoords(poly);
     * //= [[[119.32, -8.7], [119.55, -8.69], [119.51, -8.54], [119.32, -8.7]]]
     */
    function getCoords(coords) {
        if (Array.isArray(coords)) {
            return coords;
        }
        // Feature
        if (coords.type === "Feature") {
            if (coords.geometry !== null) {
                return coords.geometry.coordinates;
            }
        }
        else {
            // Geometry
            if (coords.coordinates) {
                return coords.coordinates;
            }
        }
        throw new Error("coords must be GeoJSON Feature, Geometry Object or an Array");
    }
    /**
     * Get Geometry from Feature or Geometry Object
     *
     * @param {Feature|Geometry} geojson GeoJSON Feature or Geometry Object
     * @returns {Geometry|null} GeoJSON Geometry Object
     * @throws {Error} if geojson is not a Feature or Geometry Object
     * @example
     * var point = {
     *   "type": "Feature",
     *   "properties": {},
     *   "geometry": {
     *     "type": "Point",
     *     "coordinates": [110, 40]
     *   }
     * }
     * var geom = turf.getGeom(point)
     * //={"type": "Point", "coordinates": [110, 40]}
     */
    function getGeom(geojson) {
        if (geojson.type === "Feature") {
            return geojson.geometry;
        }
        return geojson;
    }

    //http://en.wikipedia.org/wiki/Haversine_formula
    //http://www.movable-type.co.uk/scripts/latlong.html
    /**
     * Calculates the distance between two {@link Point|points} in degrees, radians, miles, or kilometers.
     * This uses the [Haversine formula](http://en.wikipedia.org/wiki/Haversine_formula) to account for global curvature.
     *
     * @name distance
     * @param {Coord | Point} from origin point or coordinate
     * @param {Coord | Point} to destination point or coordinate
     * @param {Object} [options={}] Optional parameters
     * @param {string} [options.units='kilometers'] can be degrees, radians, miles, or kilometers
     * @returns {number} distance between the two points
     * @example
     * var from = turf.point([-75.343, 39.984]);
     * var to = turf.point([-75.534, 39.123]);
     * var options = {units: 'miles'};
     *
     * var distance = turf.distance(from, to, options);
     *
     * //addToMap
     * var addToMap = [from, to];
     * from.properties.distance = distance;
     * to.properties.distance = distance;
     */
    function distance(from, to, options) {
        if (options === void 0) { options = {}; }
        var coordinates1 = getCoord$1(from);
        var coordinates2 = getCoord$1(to);
        var dLat = degreesToRadians$2(coordinates2[1] - coordinates1[1]);
        var dLon = degreesToRadians$2(coordinates2[0] - coordinates1[0]);
        var lat1 = degreesToRadians$2(coordinates1[1]);
        var lat2 = degreesToRadians$2(coordinates2[1]);
        var a = Math.pow(Math.sin(dLat / 2), 2) +
            Math.pow(Math.sin(dLon / 2), 2) * Math.cos(lat1) * Math.cos(lat2);
        return radiansToLength$1(2 * Math.atan2(Math.sqrt(a), Math.sqrt(1 - a)), options.units);
    }

    // index.ts
    var earthRadius$1 = 63710088e-1;
    var factors$1 = {
      centimeters: earthRadius$1 * 100,
      centimetres: earthRadius$1 * 100,
      degrees: 360 / (2 * Math.PI),
      feet: earthRadius$1 * 3.28084,
      inches: earthRadius$1 * 39.37,
      kilometers: earthRadius$1 / 1e3,
      kilometres: earthRadius$1 / 1e3,
      meters: earthRadius$1,
      metres: earthRadius$1,
      miles: earthRadius$1 / 1609.344,
      millimeters: earthRadius$1 * 1e3,
      millimetres: earthRadius$1 * 1e3,
      nauticalmiles: earthRadius$1 / 1852,
      radians: 1,
      yards: earthRadius$1 * 1.0936
    };
    function feature$1(geom, properties, options = {}) {
      const feat = { type: "Feature" };
      if (options.id === 0 || options.id) {
        feat.id = options.id;
      }
      if (options.bbox) {
        feat.bbox = options.bbox;
      }
      feat.properties = properties || {};
      feat.geometry = geom;
      return feat;
    }
    function point$1(coordinates, properties, options = {}) {
      if (!coordinates) {
        throw new Error("coordinates is required");
      }
      if (!Array.isArray(coordinates)) {
        throw new Error("coordinates must be an Array");
      }
      if (coordinates.length < 2) {
        throw new Error("coordinates must be at least 2 numbers long");
      }
      if (!isNumber$1(coordinates[0]) || !isNumber$1(coordinates[1])) {
        throw new Error("coordinates must contain numbers");
      }
      const geom = {
        type: "Point",
        coordinates
      };
      return feature$1(geom, properties, options);
    }
    function polygon$1(coordinates, properties, options = {}) {
      for (const ring of coordinates) {
        if (ring.length < 4) {
          throw new Error(
            "Each LinearRing of a Polygon must have 4 or more Positions."
          );
        }
        if (ring[ring.length - 1].length !== ring[0].length) {
          throw new Error("First and last Position are not equivalent.");
        }
        for (let j = 0; j < ring[ring.length - 1].length; j++) {
          if (ring[ring.length - 1][j] !== ring[0][j]) {
            throw new Error("First and last Position are not equivalent.");
          }
        }
      }
      const geom = {
        type: "Polygon",
        coordinates
      };
      return feature$1(geom, properties, options);
    }
    function lengthToRadians$1(distance, units = "kilometers") {
      const factor = factors$1[units];
      if (!factor) {
        throw new Error(units + " units is invalid");
      }
      return distance / factor;
    }
    function radiansToDegrees$1(radians) {
      const degrees = radians % (2 * Math.PI);
      return degrees * 180 / Math.PI;
    }
    function degreesToRadians$1(degrees) {
      const radians = degrees % 360;
      return radians * Math.PI / 180;
    }
    function isNumber$1(num) {
      return !isNaN(num) && num !== null && !Array.isArray(num);
    }

    // index.ts
    function getCoord(coord) {
      if (!coord) {
        throw new Error("coord is required");
      }
      if (!Array.isArray(coord)) {
        if (coord.type === "Feature" && coord.geometry !== null && coord.geometry.type === "Point") {
          return [...coord.geometry.coordinates];
        }
        if (coord.type === "Point") {
          return [...coord.coordinates];
        }
      }
      if (Array.isArray(coord) && coord.length >= 2 && !Array.isArray(coord[0]) && !Array.isArray(coord[1])) {
        return [...coord];
      }
      throw new Error("coord must be GeoJSON Point or an Array of numbers");
    }

    // index.ts
    function destination$1(origin, distance, bearing, options = {}) {
      const coordinates1 = getCoord(origin);
      const longitude1 = degreesToRadians$1(coordinates1[0]);
      const latitude1 = degreesToRadians$1(coordinates1[1]);
      const bearingRad = degreesToRadians$1(bearing);
      const radians = lengthToRadians$1(distance, options.units);
      const latitude2 = Math.asin(
        Math.sin(latitude1) * Math.cos(radians) + Math.cos(latitude1) * Math.sin(radians) * Math.cos(bearingRad)
      );
      const longitude2 = longitude1 + Math.atan2(
        Math.sin(bearingRad) * Math.sin(radians) * Math.cos(latitude1),
        Math.cos(radians) - Math.sin(latitude1) * Math.sin(latitude2)
      );
      const lng = radiansToDegrees$1(longitude2);
      const lat = radiansToDegrees$1(latitude2);
      return point$1([lng, lat], options.properties);
    }

    // index.ts
    function circle(center, radius, options = {}) {
      const steps = options.steps || 64;
      const properties = options.properties ? options.properties : !Array.isArray(center) && center.type === "Feature" && center.properties ? center.properties : {};
      const coordinates = [];
      for (let i = 0; i < steps; i++) {
        coordinates.push(
          destination$1(center, radius, i * -360 / steps, options).geometry.coordinates
        );
      }
      coordinates.push(coordinates[0]);
      return polygon$1([coordinates], properties);
    }

    /**
     * Callback for coordEach
     *
     * @callback coordEachCallback
     * @param {Array<number>} currentCoord The current coordinate being processed.
     * @param {number} coordIndex The current index of the coordinate being processed.
     * @param {number} featureIndex The current index of the Feature being processed.
     * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
     * @param {number} geometryIndex The current index of the Geometry being processed.
     */

    /**
     * Iterate over coordinates in any GeoJSON object, similar to Array.forEach()
     *
     * @name coordEach
     * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
     * @param {Function} callback a method that takes (currentCoord, coordIndex, featureIndex, multiFeatureIndex)
     * @param {boolean} [excludeWrapCoord=false] whether or not to include the final coordinate of LinearRings that wraps the ring in its iteration.
     * @returns {void}
     * @example
     * var features = turf.featureCollection([
     *   turf.point([26, 37], {"foo": "bar"}),
     *   turf.point([36, 53], {"hello": "world"})
     * ]);
     *
     * turf.coordEach(features, function (currentCoord, coordIndex, featureIndex, multiFeatureIndex, geometryIndex) {
     *   //=currentCoord
     *   //=coordIndex
     *   //=featureIndex
     *   //=multiFeatureIndex
     *   //=geometryIndex
     * });
     */
    function coordEach$2(geojson, callback, excludeWrapCoord) {
      // Handles null Geometry -- Skips this GeoJSON
      if (geojson === null) return;
      var j,
        k,
        l,
        geometry,
        stopG,
        coords,
        geometryMaybeCollection,
        wrapShrink = 0,
        coordIndex = 0,
        isGeometryCollection,
        type = geojson.type,
        isFeatureCollection = type === "FeatureCollection",
        isFeature = type === "Feature",
        stop = isFeatureCollection ? geojson.features.length : 1;

      // This logic may look a little weird. The reason why it is that way
      // is because it's trying to be fast. GeoJSON supports multiple kinds
      // of objects at its root: FeatureCollection, Features, Geometries.
      // This function has the responsibility of handling all of them, and that
      // means that some of the `for` loops you see below actually just don't apply
      // to certain inputs. For instance, if you give this just a
      // Point geometry, then both loops are short-circuited and all we do
      // is gradually rename the input until it's called 'geometry'.
      //
      // This also aims to allocate as few resources as possible: just a
      // few numbers and booleans, rather than any temporary arrays as would
      // be required with the normalization approach.
      for (var featureIndex = 0; featureIndex < stop; featureIndex++) {
        geometryMaybeCollection = isFeatureCollection
          ? geojson.features[featureIndex].geometry
          : isFeature
          ? geojson.geometry
          : geojson;
        isGeometryCollection = geometryMaybeCollection
          ? geometryMaybeCollection.type === "GeometryCollection"
          : false;
        stopG = isGeometryCollection
          ? geometryMaybeCollection.geometries.length
          : 1;

        for (var geomIndex = 0; geomIndex < stopG; geomIndex++) {
          var multiFeatureIndex = 0;
          var geometryIndex = 0;
          geometry = isGeometryCollection
            ? geometryMaybeCollection.geometries[geomIndex]
            : geometryMaybeCollection;

          // Handles null Geometry -- Skips this geometry
          if (geometry === null) continue;
          coords = geometry.coordinates;
          var geomType = geometry.type;

          wrapShrink =
            0;

          switch (geomType) {
            case null:
              break;
            case "Point":
              if (
                callback(
                  coords,
                  coordIndex,
                  featureIndex,
                  multiFeatureIndex,
                  geometryIndex
                ) === false
              )
                return false;
              coordIndex++;
              multiFeatureIndex++;
              break;
            case "LineString":
            case "MultiPoint":
              for (j = 0; j < coords.length; j++) {
                if (
                  callback(
                    coords[j],
                    coordIndex,
                    featureIndex,
                    multiFeatureIndex,
                    geometryIndex
                  ) === false
                )
                  return false;
                coordIndex++;
                if (geomType === "MultiPoint") multiFeatureIndex++;
              }
              if (geomType === "LineString") multiFeatureIndex++;
              break;
            case "Polygon":
            case "MultiLineString":
              for (j = 0; j < coords.length; j++) {
                for (k = 0; k < coords[j].length - wrapShrink; k++) {
                  if (
                    callback(
                      coords[j][k],
                      coordIndex,
                      featureIndex,
                      multiFeatureIndex,
                      geometryIndex
                    ) === false
                  )
                    return false;
                  coordIndex++;
                }
                if (geomType === "MultiLineString") multiFeatureIndex++;
                if (geomType === "Polygon") geometryIndex++;
              }
              if (geomType === "Polygon") multiFeatureIndex++;
              break;
            case "MultiPolygon":
              for (j = 0; j < coords.length; j++) {
                geometryIndex = 0;
                for (k = 0; k < coords[j].length; k++) {
                  for (l = 0; l < coords[j][k].length - wrapShrink; l++) {
                    if (
                      callback(
                        coords[j][k][l],
                        coordIndex,
                        featureIndex,
                        multiFeatureIndex,
                        geometryIndex
                      ) === false
                    )
                      return false;
                    coordIndex++;
                  }
                  geometryIndex++;
                }
                multiFeatureIndex++;
              }
              break;
            case "GeometryCollection":
              for (j = 0; j < geometry.geometries.length; j++)
                if (
                  coordEach$2(geometry.geometries[j], callback) ===
                  false
                )
                  return false;
              break;
            default:
              throw new Error("Unknown Geometry Type");
          }
        }
      }
    }

    /**
     * Callback for featureEach
     *
     * @callback featureEachCallback
     * @param {Feature<any>} currentFeature The current Feature being processed.
     * @param {number} featureIndex The current index of the Feature being processed.
     */

    /**
     * Iterate over features in any GeoJSON object, similar to
     * Array.forEach.
     *
     * @name featureEach
     * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
     * @param {Function} callback a method that takes (currentFeature, featureIndex)
     * @returns {void}
     * @example
     * var features = turf.featureCollection([
     *   turf.point([26, 37], {foo: 'bar'}),
     *   turf.point([36, 53], {hello: 'world'})
     * ]);
     *
     * turf.featureEach(features, function (currentFeature, featureIndex) {
     *   //=currentFeature
     *   //=featureIndex
     * });
     */
    function featureEach$3(geojson, callback) {
      if (geojson.type === "Feature") {
        callback(geojson, 0);
      } else if (geojson.type === "FeatureCollection") {
        for (var i = 0; i < geojson.features.length; i++) {
          if (callback(geojson.features[i], i) === false) break;
        }
      }
    }

    /**
     * Get all coordinates from any GeoJSON object.
     *
     * @name coordAll
     * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
     * @returns {Array<Array<number>>} coordinate position array
     * @example
     * var features = turf.featureCollection([
     *   turf.point([26, 37], {foo: 'bar'}),
     *   turf.point([36, 53], {hello: 'world'})
     * ]);
     *
     * var coords = turf.coordAll(features);
     * //= [[26, 37], [36, 53]]
     */
    function coordAll$2(geojson) {
      var coords = [];
      coordEach$2(geojson, function (coord) {
        coords.push(coord);
      });
      return coords;
    }

    /**
     * Callback for geomEach
     *
     * @callback geomEachCallback
     * @param {Geometry} currentGeometry The current Geometry being processed.
     * @param {number} featureIndex The current index of the Feature being processed.
     * @param {Object} featureProperties The current Feature Properties being processed.
     * @param {Array<number>} featureBBox The current Feature BBox being processed.
     * @param {number|string} featureId The current Feature Id being processed.
     */

    /**
     * Iterate over each geometry in any GeoJSON object, similar to Array.forEach()
     *
     * @name geomEach
     * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
     * @param {Function} callback a method that takes (currentGeometry, featureIndex, featureProperties, featureBBox, featureId)
     * @returns {void}
     * @example
     * var features = turf.featureCollection([
     *     turf.point([26, 37], {foo: 'bar'}),
     *     turf.point([36, 53], {hello: 'world'})
     * ]);
     *
     * turf.geomEach(features, function (currentGeometry, featureIndex, featureProperties, featureBBox, featureId) {
     *   //=currentGeometry
     *   //=featureIndex
     *   //=featureProperties
     *   //=featureBBox
     *   //=featureId
     * });
     */
    function geomEach$2(geojson, callback) {
      var i,
        j,
        g,
        geometry,
        stopG,
        geometryMaybeCollection,
        isGeometryCollection,
        featureProperties,
        featureBBox,
        featureId,
        featureIndex = 0,
        isFeatureCollection = geojson.type === "FeatureCollection",
        isFeature = geojson.type === "Feature",
        stop = isFeatureCollection ? geojson.features.length : 1;

      // This logic may look a little weird. The reason why it is that way
      // is because it's trying to be fast. GeoJSON supports multiple kinds
      // of objects at its root: FeatureCollection, Features, Geometries.
      // This function has the responsibility of handling all of them, and that
      // means that some of the `for` loops you see below actually just don't apply
      // to certain inputs. For instance, if you give this just a
      // Point geometry, then both loops are short-circuited and all we do
      // is gradually rename the input until it's called 'geometry'.
      //
      // This also aims to allocate as few resources as possible: just a
      // few numbers and booleans, rather than any temporary arrays as would
      // be required with the normalization approach.
      for (i = 0; i < stop; i++) {
        geometryMaybeCollection = isFeatureCollection
          ? geojson.features[i].geometry
          : isFeature
          ? geojson.geometry
          : geojson;
        featureProperties = isFeatureCollection
          ? geojson.features[i].properties
          : isFeature
          ? geojson.properties
          : {};
        featureBBox = isFeatureCollection
          ? geojson.features[i].bbox
          : isFeature
          ? geojson.bbox
          : undefined;
        featureId = isFeatureCollection
          ? geojson.features[i].id
          : isFeature
          ? geojson.id
          : undefined;
        isGeometryCollection = geometryMaybeCollection
          ? geometryMaybeCollection.type === "GeometryCollection"
          : false;
        stopG = isGeometryCollection
          ? geometryMaybeCollection.geometries.length
          : 1;

        for (g = 0; g < stopG; g++) {
          geometry = isGeometryCollection
            ? geometryMaybeCollection.geometries[g]
            : geometryMaybeCollection;

          // Handle null Geometry
          if (geometry === null) {
            if (
              callback(
                null,
                featureIndex,
                featureProperties,
                featureBBox,
                featureId
              ) === false
            )
              return false;
            continue;
          }
          switch (geometry.type) {
            case "Point":
            case "LineString":
            case "MultiPoint":
            case "Polygon":
            case "MultiLineString":
            case "MultiPolygon": {
              if (
                callback(
                  geometry,
                  featureIndex,
                  featureProperties,
                  featureBBox,
                  featureId
                ) === false
              )
                return false;
              break;
            }
            case "GeometryCollection": {
              for (j = 0; j < geometry.geometries.length; j++) {
                if (
                  callback(
                    geometry.geometries[j],
                    featureIndex,
                    featureProperties,
                    featureBBox,
                    featureId
                  ) === false
                )
                  return false;
              }
              break;
            }
            default:
              throw new Error("Unknown Geometry Type");
          }
        }
        // Only increase `featureIndex` per each feature
        featureIndex++;
      }
    }

    /**
     * Callback for flattenEach
     *
     * @callback flattenEachCallback
     * @param {Feature} currentFeature The current flattened feature being processed.
     * @param {number} featureIndex The current index of the Feature being processed.
     * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
     */

    /**
     * Iterate over flattened features in any GeoJSON object, similar to
     * Array.forEach.
     *
     * @name flattenEach
     * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
     * @param {Function} callback a method that takes (currentFeature, featureIndex, multiFeatureIndex)
     * @example
     * var features = turf.featureCollection([
     *     turf.point([26, 37], {foo: 'bar'}),
     *     turf.multiPoint([[40, 30], [36, 53]], {hello: 'world'})
     * ]);
     *
     * turf.flattenEach(features, function (currentFeature, featureIndex, multiFeatureIndex) {
     *   //=currentFeature
     *   //=featureIndex
     *   //=multiFeatureIndex
     * });
     */
    function flattenEach$2(geojson, callback) {
      geomEach$2(geojson, function (geometry, featureIndex, properties, bbox, id) {
        // Callback for single geometry
        var type = geometry === null ? null : geometry.type;
        switch (type) {
          case null:
          case "Point":
          case "LineString":
          case "Polygon":
            if (
              callback(
                feature$2(geometry, properties, { bbox: bbox, id: id }),
                featureIndex,
                0
              ) === false
            )
              return false;
            return;
        }

        var geomType;

        // Callback for multi-geometry
        switch (type) {
          case "MultiPoint":
            geomType = "Point";
            break;
          case "MultiLineString":
            geomType = "LineString";
            break;
          case "MultiPolygon":
            geomType = "Polygon";
            break;
        }

        for (
          var multiFeatureIndex = 0;
          multiFeatureIndex < geometry.coordinates.length;
          multiFeatureIndex++
        ) {
          var coordinate = geometry.coordinates[multiFeatureIndex];
          var geom = {
            type: geomType,
            coordinates: coordinate,
          };
          if (
            callback(feature$2(geom, properties), featureIndex, multiFeatureIndex) ===
            false
          )
            return false;
        }
      });
    }

    /**
     * @private
     */
    function polygonToLine(poly, options) {
        if (options === void 0) { options = {}; }
        var geom = getGeom(poly);
        var coords = geom.coordinates;
        var properties = options.properties
            ? options.properties
            : poly.type === "Feature"
                ? poly.properties
                : {};
        return coordsToLine(coords, properties);
    }
    /**
     * @private
     */
    function coordsToLine(coords, properties) {
        if (coords.length > 1) {
            return multiLineString$1(coords, properties);
        }
        return lineString$1(coords[0], properties);
    }

    // http://en.wikipedia.org/wiki/Haversine_formula
    // http://www.movable-type.co.uk/scripts/latlong.html
    /**
     * Takes two {@link Point|points} and finds the geographic bearing between them,
     * i.e. the angle measured in degrees from the north line (0 degrees)
     *
     * @name bearing
     * @param {Coord} start starting Point
     * @param {Coord} end ending Point
     * @param {Object} [options={}] Optional parameters
     * @param {boolean} [options.final=false] calculates the final bearing if true
     * @returns {number} bearing in decimal degrees, between -180 and 180 degrees (positive clockwise)
     * @example
     * var point1 = turf.point([-75.343, 39.984]);
     * var point2 = turf.point([-75.534, 39.123]);
     *
     * var bearing = turf.bearing(point1, point2);
     *
     * //addToMap
     * var addToMap = [point1, point2]
     * point1.properties['marker-color'] = '#f00'
     * point2.properties['marker-color'] = '#0f0'
     * point1.properties.bearing = bearing
     */
    function bearing(start, end, options) {
        if (options === void 0) { options = {}; }
        // Reverse calculation
        if (options.final === true) {
            return calculateFinalBearing(start, end);
        }
        var coordinates1 = getCoord$1(start);
        var coordinates2 = getCoord$1(end);
        var lon1 = degreesToRadians$2(coordinates1[0]);
        var lon2 = degreesToRadians$2(coordinates2[0]);
        var lat1 = degreesToRadians$2(coordinates1[1]);
        var lat2 = degreesToRadians$2(coordinates2[1]);
        var a = Math.sin(lon2 - lon1) * Math.cos(lat2);
        var b = Math.cos(lat1) * Math.sin(lat2) -
            Math.sin(lat1) * Math.cos(lat2) * Math.cos(lon2 - lon1);
        return radiansToDegrees$2(Math.atan2(a, b));
    }
    /**
     * Calculates Final Bearing
     *
     * @private
     * @param {Coord} start starting Point
     * @param {Coord} end ending Point
     * @returns {number} bearing
     */
    function calculateFinalBearing(start, end) {
        // Swap start & end
        var bear = bearing(end, start);
        bear = (bear + 180) % 360;
        return bear;
    }

    // http://en.wikipedia.org/wiki/Haversine_formula
    // http://www.movable-type.co.uk/scripts/latlong.html
    /**
     * Takes a {@link Point} and calculates the location of a destination point given a distance in
     * degrees, radians, miles, or kilometers; and bearing in degrees.
     * This uses the [Haversine formula](http://en.wikipedia.org/wiki/Haversine_formula) to account for global curvature.
     *
     * @name destination
     * @param {Coord} origin starting point
     * @param {number} distance distance from the origin point
     * @param {number} bearing ranging from -180 to 180
     * @param {Object} [options={}] Optional parameters
     * @param {string} [options.units='kilometers'] miles, kilometers, degrees, or radians
     * @param {Object} [options.properties={}] Translate properties to Point
     * @returns {Feature<Point>} destination point
     * @example
     * var point = turf.point([-75.343, 39.984]);
     * var distance = 50;
     * var bearing = 90;
     * var options = {units: 'miles'};
     *
     * var destination = turf.destination(point, distance, bearing, options);
     *
     * //addToMap
     * var addToMap = [point, destination]
     * destination.properties['marker-color'] = '#f00';
     * point.properties['marker-color'] = '#0f0';
     */
    function destination(origin, distance, bearing, options) {
        if (options === void 0) { options = {}; }
        // Handle input
        var coordinates1 = getCoord$1(origin);
        var longitude1 = degreesToRadians$2(coordinates1[0]);
        var latitude1 = degreesToRadians$2(coordinates1[1]);
        var bearingRad = degreesToRadians$2(bearing);
        var radians = lengthToRadians$2(distance, options.units);
        // Main
        var latitude2 = Math.asin(Math.sin(latitude1) * Math.cos(radians) +
            Math.cos(latitude1) * Math.sin(radians) * Math.cos(bearingRad));
        var longitude2 = longitude1 +
            Math.atan2(Math.sin(bearingRad) * Math.sin(radians) * Math.cos(latitude1), Math.cos(radians) - Math.sin(latitude1) * Math.sin(latitude2));
        var lng = radiansToDegrees$2(longitude2);
        var lat = radiansToDegrees$2(latitude2);
        return point$2([lng, lat], options.properties);
    }

    /**
     * Creates a {@link FeatureCollection} of 2-vertex {@link LineString} segments from a
     * {@link LineString|(Multi)LineString} or {@link Polygon|(Multi)Polygon}.
     *
     * @name lineSegment
     * @param {GeoJSON} geojson GeoJSON Polygon or LineString
     * @returns {FeatureCollection<LineString>} 2-vertex line segments
     * @example
     * var polygon = turf.polygon([[[-50, 5], [-40, -10], [-50, -10], [-40, 5], [-50, 5]]]);
     * var segments = turf.lineSegment(polygon);
     *
     * //addToMap
     * var addToMap = [polygon, segments]
     */
    function lineSegment(geojson) {
        if (!geojson) {
            throw new Error("geojson is required");
        }
        var results = [];
        flattenEach$2(geojson, function (feature) {
            lineSegmentFeature(feature, results);
        });
        return featureCollection$2(results);
    }
    /**
     * Line Segment
     *
     * @private
     * @param {Feature<LineString|Polygon>} geojson Line or polygon feature
     * @param {Array} results push to results
     * @returns {void}
     */
    function lineSegmentFeature(geojson, results) {
        var coords = [];
        var geometry = geojson.geometry;
        if (geometry !== null) {
            switch (geometry.type) {
                case "Polygon":
                    coords = getCoords(geometry);
                    break;
                case "LineString":
                    coords = [getCoords(geometry)];
            }
            coords.forEach(function (coord) {
                var segments = createSegments(coord, geojson.properties);
                segments.forEach(function (segment) {
                    segment.id = results.length;
                    results.push(segment);
                });
            });
        }
    }
    /**
     * Create Segments from LineString coordinates
     *
     * @private
     * @param {Array<Array<number>>} coords LineString coordinates
     * @param {*} properties GeoJSON properties
     * @returns {Array<Feature<LineString>>} line segments
     */
    function createSegments(coords, properties) {
        var segments = [];
        coords.reduce(function (previousCoords, currentCoords) {
            var segment = lineString$1([previousCoords, currentCoords], properties);
            segment.bbox = bbox$1(previousCoords, currentCoords);
            segments.push(segment);
            return currentCoords;
        });
        return segments;
    }
    /**
     * Create BBox between two coordinates (faster than @turf/bbox)
     *
     * @private
     * @param {Array<number>} coords1 Point coordinate
     * @param {Array<number>} coords2 Point coordinate
     * @returns {BBox} [west, south, east, north]
     */
    function bbox$1(coords1, coords2) {
        var x1 = coords1[0];
        var y1 = coords1[1];
        var x2 = coords2[0];
        var y2 = coords2[1];
        var west = x1 < x2 ? x1 : x2;
        var south = y1 < y2 ? y1 : y2;
        var east = x1 > x2 ? x1 : x2;
        var north = y1 > y2 ? y1 : y2;
        return [west, south, east, north];
    }

    function getDefaultExportFromCjs (x) {
    	return x && x.__esModule && Object.prototype.hasOwnProperty.call(x, 'default') ? x['default'] : x;
    }

    function getAugmentedNamespace(n) {
      if (n.__esModule) return n;
      var f = n.default;
    	if (typeof f == "function") {
    		var a = function a () {
    			if (this instanceof a) {
            return Reflect.construct(f, arguments, this.constructor);
    			}
    			return f.apply(this, arguments);
    		};
    		a.prototype = f.prototype;
      } else a = {};
      Object.defineProperty(a, '__esModule', {value: true});
    	Object.keys(n).forEach(function (k) {
    		var d = Object.getOwnPropertyDescriptor(n, k);
    		Object.defineProperty(a, k, d.get ? d : {
    			enumerable: true,
    			get: function () {
    				return n[k];
    			}
    		});
    	});
    	return a;
    }

    var geojsonRbush$1 = {exports: {}};

    function quickselect(arr, k, left, right, compare) {
        quickselectStep(arr, k, left || 0, right || (arr.length - 1), compare || defaultCompare);
    }

    function quickselectStep(arr, k, left, right, compare) {

        while (right > left) {
            if (right - left > 600) {
                var n = right - left + 1;
                var m = k - left + 1;
                var z = Math.log(n);
                var s = 0.5 * Math.exp(2 * z / 3);
                var sd = 0.5 * Math.sqrt(z * s * (n - s) / n) * (m - n / 2 < 0 ? -1 : 1);
                var newLeft = Math.max(left, Math.floor(k - m * s / n + sd));
                var newRight = Math.min(right, Math.floor(k + (n - m) * s / n + sd));
                quickselectStep(arr, k, newLeft, newRight, compare);
            }

            var t = arr[k];
            var i = left;
            var j = right;

            swap(arr, left, k);
            if (compare(arr[right], t) > 0) swap(arr, left, right);

            while (i < j) {
                swap(arr, i, j);
                i++;
                j--;
                while (compare(arr[i], t) < 0) i++;
                while (compare(arr[j], t) > 0) j--;
            }

            if (compare(arr[left], t) === 0) swap(arr, left, j);
            else {
                j++;
                swap(arr, j, right);
            }

            if (j <= k) left = j + 1;
            if (k <= j) right = j - 1;
        }
    }

    function swap(arr, i, j) {
        var tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }

    function defaultCompare(a, b) {
        return a < b ? -1 : a > b ? 1 : 0;
    }

    class RBush {
        constructor(maxEntries = 9) {
            // max entries in a node is 9 by default; min node fill is 40% for best performance
            this._maxEntries = Math.max(4, maxEntries);
            this._minEntries = Math.max(2, Math.ceil(this._maxEntries * 0.4));
            this.clear();
        }

        all() {
            return this._all(this.data, []);
        }

        search(bbox) {
            let node = this.data;
            const result = [];

            if (!intersects$1(bbox, node)) return result;

            const toBBox = this.toBBox;
            const nodesToSearch = [];

            while (node) {
                for (let i = 0; i < node.children.length; i++) {
                    const child = node.children[i];
                    const childBBox = node.leaf ? toBBox(child) : child;

                    if (intersects$1(bbox, childBBox)) {
                        if (node.leaf) result.push(child);
                        else if (contains(bbox, childBBox)) this._all(child, result);
                        else nodesToSearch.push(child);
                    }
                }
                node = nodesToSearch.pop();
            }

            return result;
        }

        collides(bbox) {
            let node = this.data;

            if (!intersects$1(bbox, node)) return false;

            const nodesToSearch = [];
            while (node) {
                for (let i = 0; i < node.children.length; i++) {
                    const child = node.children[i];
                    const childBBox = node.leaf ? this.toBBox(child) : child;

                    if (intersects$1(bbox, childBBox)) {
                        if (node.leaf || contains(bbox, childBBox)) return true;
                        nodesToSearch.push(child);
                    }
                }
                node = nodesToSearch.pop();
            }

            return false;
        }

        load(data) {
            if (!(data && data.length)) return this;

            if (data.length < this._minEntries) {
                for (let i = 0; i < data.length; i++) {
                    this.insert(data[i]);
                }
                return this;
            }

            // recursively build the tree with the given data from scratch using OMT algorithm
            let node = this._build(data.slice(), 0, data.length - 1, 0);

            if (!this.data.children.length) {
                // save as is if tree is empty
                this.data = node;

            } else if (this.data.height === node.height) {
                // split root if trees have the same height
                this._splitRoot(this.data, node);

            } else {
                if (this.data.height < node.height) {
                    // swap trees if inserted one is bigger
                    const tmpNode = this.data;
                    this.data = node;
                    node = tmpNode;
                }

                // insert the small tree into the large tree at appropriate level
                this._insert(node, this.data.height - node.height - 1, true);
            }

            return this;
        }

        insert(item) {
            if (item) this._insert(item, this.data.height - 1);
            return this;
        }

        clear() {
            this.data = createNode([]);
            return this;
        }

        remove(item, equalsFn) {
            if (!item) return this;

            let node = this.data;
            const bbox = this.toBBox(item);
            const path = [];
            const indexes = [];
            let i, parent, goingUp;

            // depth-first iterative tree traversal
            while (node || path.length) {

                if (!node) { // go up
                    node = path.pop();
                    parent = path[path.length - 1];
                    i = indexes.pop();
                    goingUp = true;
                }

                if (node.leaf) { // check current node
                    const index = findItem(item, node.children, equalsFn);

                    if (index !== -1) {
                        // item found, remove the item and condense tree upwards
                        node.children.splice(index, 1);
                        path.push(node);
                        this._condense(path);
                        return this;
                    }
                }

                if (!goingUp && !node.leaf && contains(node, bbox)) { // go down
                    path.push(node);
                    indexes.push(i);
                    i = 0;
                    parent = node;
                    node = node.children[0];

                } else if (parent) { // go right
                    i++;
                    node = parent.children[i];
                    goingUp = false;

                } else node = null; // nothing found
            }

            return this;
        }

        toBBox(item) { return item; }

        compareMinX(a, b) { return a.minX - b.minX; }
        compareMinY(a, b) { return a.minY - b.minY; }

        toJSON() { return this.data; }

        fromJSON(data) {
            this.data = data;
            return this;
        }

        _all(node, result) {
            const nodesToSearch = [];
            while (node) {
                if (node.leaf) result.push(...node.children);
                else nodesToSearch.push(...node.children);

                node = nodesToSearch.pop();
            }
            return result;
        }

        _build(items, left, right, height) {

            const N = right - left + 1;
            let M = this._maxEntries;
            let node;

            if (N <= M) {
                // reached leaf level; return leaf
                node = createNode(items.slice(left, right + 1));
                calcBBox(node, this.toBBox);
                return node;
            }

            if (!height) {
                // target height of the bulk-loaded tree
                height = Math.ceil(Math.log(N) / Math.log(M));

                // target number of root entries to maximize storage utilization
                M = Math.ceil(N / Math.pow(M, height - 1));
            }

            node = createNode([]);
            node.leaf = false;
            node.height = height;

            // split the items into M mostly square tiles

            const N2 = Math.ceil(N / M);
            const N1 = N2 * Math.ceil(Math.sqrt(M));

            multiSelect(items, left, right, N1, this.compareMinX);

            for (let i = left; i <= right; i += N1) {

                const right2 = Math.min(i + N1 - 1, right);

                multiSelect(items, i, right2, N2, this.compareMinY);

                for (let j = i; j <= right2; j += N2) {

                    const right3 = Math.min(j + N2 - 1, right2);

                    // pack each entry recursively
                    node.children.push(this._build(items, j, right3, height - 1));
                }
            }

            calcBBox(node, this.toBBox);

            return node;
        }

        _chooseSubtree(bbox, node, level, path) {
            while (true) {
                path.push(node);

                if (node.leaf || path.length - 1 === level) break;

                let minArea = Infinity;
                let minEnlargement = Infinity;
                let targetNode;

                for (let i = 0; i < node.children.length; i++) {
                    const child = node.children[i];
                    const area = bboxArea(child);
                    const enlargement = enlargedArea(bbox, child) - area;

                    // choose entry with the least area enlargement
                    if (enlargement < minEnlargement) {
                        minEnlargement = enlargement;
                        minArea = area < minArea ? area : minArea;
                        targetNode = child;

                    } else if (enlargement === minEnlargement) {
                        // otherwise choose one with the smallest area
                        if (area < minArea) {
                            minArea = area;
                            targetNode = child;
                        }
                    }
                }

                node = targetNode || node.children[0];
            }

            return node;
        }

        _insert(item, level, isNode) {
            const bbox = isNode ? item : this.toBBox(item);
            const insertPath = [];

            // find the best node for accommodating the item, saving all nodes along the path too
            const node = this._chooseSubtree(bbox, this.data, level, insertPath);

            // put the item into the node
            node.children.push(item);
            extend(node, bbox);

            // split on node overflow; propagate upwards if necessary
            while (level >= 0) {
                if (insertPath[level].children.length > this._maxEntries) {
                    this._split(insertPath, level);
                    level--;
                } else break;
            }

            // adjust bboxes along the insertion path
            this._adjustParentBBoxes(bbox, insertPath, level);
        }

        // split overflowed node into two
        _split(insertPath, level) {
            const node = insertPath[level];
            const M = node.children.length;
            const m = this._minEntries;

            this._chooseSplitAxis(node, m, M);

            const splitIndex = this._chooseSplitIndex(node, m, M);

            const newNode = createNode(node.children.splice(splitIndex, node.children.length - splitIndex));
            newNode.height = node.height;
            newNode.leaf = node.leaf;

            calcBBox(node, this.toBBox);
            calcBBox(newNode, this.toBBox);

            if (level) insertPath[level - 1].children.push(newNode);
            else this._splitRoot(node, newNode);
        }

        _splitRoot(node, newNode) {
            // split root node
            this.data = createNode([node, newNode]);
            this.data.height = node.height + 1;
            this.data.leaf = false;
            calcBBox(this.data, this.toBBox);
        }

        _chooseSplitIndex(node, m, M) {
            let index;
            let minOverlap = Infinity;
            let minArea = Infinity;

            for (let i = m; i <= M - m; i++) {
                const bbox1 = distBBox(node, 0, i, this.toBBox);
                const bbox2 = distBBox(node, i, M, this.toBBox);

                const overlap = intersectionArea(bbox1, bbox2);
                const area = bboxArea(bbox1) + bboxArea(bbox2);

                // choose distribution with minimum overlap
                if (overlap < minOverlap) {
                    minOverlap = overlap;
                    index = i;

                    minArea = area < minArea ? area : minArea;

                } else if (overlap === minOverlap) {
                    // otherwise choose distribution with minimum area
                    if (area < minArea) {
                        minArea = area;
                        index = i;
                    }
                }
            }

            return index || M - m;
        }

        // sorts node children by the best axis for split
        _chooseSplitAxis(node, m, M) {
            const compareMinX = node.leaf ? this.compareMinX : compareNodeMinX;
            const compareMinY = node.leaf ? this.compareMinY : compareNodeMinY;
            const xMargin = this._allDistMargin(node, m, M, compareMinX);
            const yMargin = this._allDistMargin(node, m, M, compareMinY);

            // if total distributions margin value is minimal for x, sort by minX,
            // otherwise it's already sorted by minY
            if (xMargin < yMargin) node.children.sort(compareMinX);
        }

        // total margin of all possible split distributions where each node is at least m full
        _allDistMargin(node, m, M, compare) {
            node.children.sort(compare);

            const toBBox = this.toBBox;
            const leftBBox = distBBox(node, 0, m, toBBox);
            const rightBBox = distBBox(node, M - m, M, toBBox);
            let margin = bboxMargin(leftBBox) + bboxMargin(rightBBox);

            for (let i = m; i < M - m; i++) {
                const child = node.children[i];
                extend(leftBBox, node.leaf ? toBBox(child) : child);
                margin += bboxMargin(leftBBox);
            }

            for (let i = M - m - 1; i >= m; i--) {
                const child = node.children[i];
                extend(rightBBox, node.leaf ? toBBox(child) : child);
                margin += bboxMargin(rightBBox);
            }

            return margin;
        }

        _adjustParentBBoxes(bbox, path, level) {
            // adjust bboxes along the given tree path
            for (let i = level; i >= 0; i--) {
                extend(path[i], bbox);
            }
        }

        _condense(path) {
            // go through the path, removing empty nodes and updating bboxes
            for (let i = path.length - 1, siblings; i >= 0; i--) {
                if (path[i].children.length === 0) {
                    if (i > 0) {
                        siblings = path[i - 1].children;
                        siblings.splice(siblings.indexOf(path[i]), 1);

                    } else this.clear();

                } else calcBBox(path[i], this.toBBox);
            }
        }
    }

    function findItem(item, items, equalsFn) {
        if (!equalsFn) return items.indexOf(item);

        for (let i = 0; i < items.length; i++) {
            if (equalsFn(item, items[i])) return i;
        }
        return -1;
    }

    // calculate node's bbox from bboxes of its children
    function calcBBox(node, toBBox) {
        distBBox(node, 0, node.children.length, toBBox, node);
    }

    // min bounding rectangle of node children from k to p-1
    function distBBox(node, k, p, toBBox, destNode) {
        if (!destNode) destNode = createNode(null);
        destNode.minX = Infinity;
        destNode.minY = Infinity;
        destNode.maxX = -Infinity;
        destNode.maxY = -Infinity;

        for (let i = k; i < p; i++) {
            const child = node.children[i];
            extend(destNode, node.leaf ? toBBox(child) : child);
        }

        return destNode;
    }

    function extend(a, b) {
        a.minX = Math.min(a.minX, b.minX);
        a.minY = Math.min(a.minY, b.minY);
        a.maxX = Math.max(a.maxX, b.maxX);
        a.maxY = Math.max(a.maxY, b.maxY);
        return a;
    }

    function compareNodeMinX(a, b) { return a.minX - b.minX; }
    function compareNodeMinY(a, b) { return a.minY - b.minY; }

    function bboxArea(a)   { return (a.maxX - a.minX) * (a.maxY - a.minY); }
    function bboxMargin(a) { return (a.maxX - a.minX) + (a.maxY - a.minY); }

    function enlargedArea(a, b) {
        return (Math.max(b.maxX, a.maxX) - Math.min(b.minX, a.minX)) *
               (Math.max(b.maxY, a.maxY) - Math.min(b.minY, a.minY));
    }

    function intersectionArea(a, b) {
        const minX = Math.max(a.minX, b.minX);
        const minY = Math.max(a.minY, b.minY);
        const maxX = Math.min(a.maxX, b.maxX);
        const maxY = Math.min(a.maxY, b.maxY);

        return Math.max(0, maxX - minX) *
               Math.max(0, maxY - minY);
    }

    function contains(a, b) {
        return a.minX <= b.minX &&
               a.minY <= b.minY &&
               b.maxX <= a.maxX &&
               b.maxY <= a.maxY;
    }

    function intersects$1(a, b) {
        return b.minX <= a.maxX &&
               b.minY <= a.maxY &&
               b.maxX >= a.minX &&
               b.maxY >= a.minY;
    }

    function createNode(children) {
        return {
            children,
            height: 1,
            leaf: true,
            minX: Infinity,
            minY: Infinity,
            maxX: -Infinity,
            maxY: -Infinity
        };
    }

    // sort an array so that items come in groups of n unsorted items, with groups sorted between each other;
    // combines selection algorithm with binary divide & conquer approach

    function multiSelect(arr, left, right, n, compare) {
        const stack = [left, right];

        while (stack.length) {
            right = stack.pop();
            left = stack.pop();

            if (right - left <= n) continue;

            const mid = left + Math.ceil((right - left) / n / 2) * n;
            quickselect(arr, mid, left, right, compare);

            stack.push(left, mid, mid, right);
        }
    }

    var rbush$2 = /*#__PURE__*/Object.freeze({
        __proto__: null,
        default: RBush
    });

    var require$$0 = /*@__PURE__*/getAugmentedNamespace(rbush$2);

    var js$1 = {};

    (function (exports) {
    	Object.defineProperty(exports, "__esModule", { value: true });
    	/**
    	 * @module helpers
    	 */
    	/**
    	 * Earth Radius used with the Harvesine formula and approximates using a spherical (non-ellipsoid) Earth.
    	 *
    	 * @memberof helpers
    	 * @type {number}
    	 */
    	exports.earthRadius = 6371008.8;
    	/**
    	 * Unit of measurement factors using a spherical (non-ellipsoid) earth radius.
    	 *
    	 * @memberof helpers
    	 * @type {Object}
    	 */
    	exports.factors = {
    	    centimeters: exports.earthRadius * 100,
    	    centimetres: exports.earthRadius * 100,
    	    degrees: exports.earthRadius / 111325,
    	    feet: exports.earthRadius * 3.28084,
    	    inches: exports.earthRadius * 39.37,
    	    kilometers: exports.earthRadius / 1000,
    	    kilometres: exports.earthRadius / 1000,
    	    meters: exports.earthRadius,
    	    metres: exports.earthRadius,
    	    miles: exports.earthRadius / 1609.344,
    	    millimeters: exports.earthRadius * 1000,
    	    millimetres: exports.earthRadius * 1000,
    	    nauticalmiles: exports.earthRadius / 1852,
    	    radians: 1,
    	    yards: exports.earthRadius * 1.0936,
    	};
    	/**
    	 * Units of measurement factors based on 1 meter.
    	 *
    	 * @memberof helpers
    	 * @type {Object}
    	 */
    	exports.unitsFactors = {
    	    centimeters: 100,
    	    centimetres: 100,
    	    degrees: 1 / 111325,
    	    feet: 3.28084,
    	    inches: 39.37,
    	    kilometers: 1 / 1000,
    	    kilometres: 1 / 1000,
    	    meters: 1,
    	    metres: 1,
    	    miles: 1 / 1609.344,
    	    millimeters: 1000,
    	    millimetres: 1000,
    	    nauticalmiles: 1 / 1852,
    	    radians: 1 / exports.earthRadius,
    	    yards: 1.0936133,
    	};
    	/**
    	 * Area of measurement factors based on 1 square meter.
    	 *
    	 * @memberof helpers
    	 * @type {Object}
    	 */
    	exports.areaFactors = {
    	    acres: 0.000247105,
    	    centimeters: 10000,
    	    centimetres: 10000,
    	    feet: 10.763910417,
    	    hectares: 0.0001,
    	    inches: 1550.003100006,
    	    kilometers: 0.000001,
    	    kilometres: 0.000001,
    	    meters: 1,
    	    metres: 1,
    	    miles: 3.86e-7,
    	    millimeters: 1000000,
    	    millimetres: 1000000,
    	    yards: 1.195990046,
    	};
    	/**
    	 * Wraps a GeoJSON {@link Geometry} in a GeoJSON {@link Feature}.
    	 *
    	 * @name feature
    	 * @param {Geometry} geometry input geometry
    	 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
    	 * @param {Object} [options={}] Optional Parameters
    	 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
    	 * @param {string|number} [options.id] Identifier associated with the Feature
    	 * @returns {Feature} a GeoJSON Feature
    	 * @example
    	 * var geometry = {
    	 *   "type": "Point",
    	 *   "coordinates": [110, 50]
    	 * };
    	 *
    	 * var feature = turf.feature(geometry);
    	 *
    	 * //=feature
    	 */
    	function feature(geom, properties, options) {
    	    if (options === void 0) { options = {}; }
    	    var feat = { type: "Feature" };
    	    if (options.id === 0 || options.id) {
    	        feat.id = options.id;
    	    }
    	    if (options.bbox) {
    	        feat.bbox = options.bbox;
    	    }
    	    feat.properties = properties || {};
    	    feat.geometry = geom;
    	    return feat;
    	}
    	exports.feature = feature;
    	/**
    	 * Creates a GeoJSON {@link Geometry} from a Geometry string type & coordinates.
    	 * For GeometryCollection type use `helpers.geometryCollection`
    	 *
    	 * @name geometry
    	 * @param {string} type Geometry Type
    	 * @param {Array<any>} coordinates Coordinates
    	 * @param {Object} [options={}] Optional Parameters
    	 * @returns {Geometry} a GeoJSON Geometry
    	 * @example
    	 * var type = "Point";
    	 * var coordinates = [110, 50];
    	 * var geometry = turf.geometry(type, coordinates);
    	 * // => geometry
    	 */
    	function geometry(type, coordinates, _options) {
    	    switch (type) {
    	        case "Point":
    	            return point(coordinates).geometry;
    	        case "LineString":
    	            return lineString(coordinates).geometry;
    	        case "Polygon":
    	            return polygon(coordinates).geometry;
    	        case "MultiPoint":
    	            return multiPoint(coordinates).geometry;
    	        case "MultiLineString":
    	            return multiLineString(coordinates).geometry;
    	        case "MultiPolygon":
    	            return multiPolygon(coordinates).geometry;
    	        default:
    	            throw new Error(type + " is invalid");
    	    }
    	}
    	exports.geometry = geometry;
    	/**
    	 * Creates a {@link Point} {@link Feature} from a Position.
    	 *
    	 * @name point
    	 * @param {Array<number>} coordinates longitude, latitude position (each in decimal degrees)
    	 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
    	 * @param {Object} [options={}] Optional Parameters
    	 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
    	 * @param {string|number} [options.id] Identifier associated with the Feature
    	 * @returns {Feature<Point>} a Point feature
    	 * @example
    	 * var point = turf.point([-75.343, 39.984]);
    	 *
    	 * //=point
    	 */
    	function point(coordinates, properties, options) {
    	    if (options === void 0) { options = {}; }
    	    if (!coordinates) {
    	        throw new Error("coordinates is required");
    	    }
    	    if (!Array.isArray(coordinates)) {
    	        throw new Error("coordinates must be an Array");
    	    }
    	    if (coordinates.length < 2) {
    	        throw new Error("coordinates must be at least 2 numbers long");
    	    }
    	    if (!isNumber(coordinates[0]) || !isNumber(coordinates[1])) {
    	        throw new Error("coordinates must contain numbers");
    	    }
    	    var geom = {
    	        type: "Point",
    	        coordinates: coordinates,
    	    };
    	    return feature(geom, properties, options);
    	}
    	exports.point = point;
    	/**
    	 * Creates a {@link Point} {@link FeatureCollection} from an Array of Point coordinates.
    	 *
    	 * @name points
    	 * @param {Array<Array<number>>} coordinates an array of Points
    	 * @param {Object} [properties={}] Translate these properties to each Feature
    	 * @param {Object} [options={}] Optional Parameters
    	 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north]
    	 * associated with the FeatureCollection
    	 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
    	 * @returns {FeatureCollection<Point>} Point Feature
    	 * @example
    	 * var points = turf.points([
    	 *   [-75, 39],
    	 *   [-80, 45],
    	 *   [-78, 50]
    	 * ]);
    	 *
    	 * //=points
    	 */
    	function points(coordinates, properties, options) {
    	    if (options === void 0) { options = {}; }
    	    return featureCollection(coordinates.map(function (coords) {
    	        return point(coords, properties);
    	    }), options);
    	}
    	exports.points = points;
    	/**
    	 * Creates a {@link Polygon} {@link Feature} from an Array of LinearRings.
    	 *
    	 * @name polygon
    	 * @param {Array<Array<Array<number>>>} coordinates an array of LinearRings
    	 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
    	 * @param {Object} [options={}] Optional Parameters
    	 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
    	 * @param {string|number} [options.id] Identifier associated with the Feature
    	 * @returns {Feature<Polygon>} Polygon Feature
    	 * @example
    	 * var polygon = turf.polygon([[[-5, 52], [-4, 56], [-2, 51], [-7, 54], [-5, 52]]], { name: 'poly1' });
    	 *
    	 * //=polygon
    	 */
    	function polygon(coordinates, properties, options) {
    	    if (options === void 0) { options = {}; }
    	    for (var _i = 0, coordinates_1 = coordinates; _i < coordinates_1.length; _i++) {
    	        var ring = coordinates_1[_i];
    	        if (ring.length < 4) {
    	            throw new Error("Each LinearRing of a Polygon must have 4 or more Positions.");
    	        }
    	        for (var j = 0; j < ring[ring.length - 1].length; j++) {
    	            // Check if first point of Polygon contains two numbers
    	            if (ring[ring.length - 1][j] !== ring[0][j]) {
    	                throw new Error("First and last Position are not equivalent.");
    	            }
    	        }
    	    }
    	    var geom = {
    	        type: "Polygon",
    	        coordinates: coordinates,
    	    };
    	    return feature(geom, properties, options);
    	}
    	exports.polygon = polygon;
    	/**
    	 * Creates a {@link Polygon} {@link FeatureCollection} from an Array of Polygon coordinates.
    	 *
    	 * @name polygons
    	 * @param {Array<Array<Array<Array<number>>>>} coordinates an array of Polygon coordinates
    	 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
    	 * @param {Object} [options={}] Optional Parameters
    	 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
    	 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
    	 * @returns {FeatureCollection<Polygon>} Polygon FeatureCollection
    	 * @example
    	 * var polygons = turf.polygons([
    	 *   [[[-5, 52], [-4, 56], [-2, 51], [-7, 54], [-5, 52]]],
    	 *   [[[-15, 42], [-14, 46], [-12, 41], [-17, 44], [-15, 42]]],
    	 * ]);
    	 *
    	 * //=polygons
    	 */
    	function polygons(coordinates, properties, options) {
    	    if (options === void 0) { options = {}; }
    	    return featureCollection(coordinates.map(function (coords) {
    	        return polygon(coords, properties);
    	    }), options);
    	}
    	exports.polygons = polygons;
    	/**
    	 * Creates a {@link LineString} {@link Feature} from an Array of Positions.
    	 *
    	 * @name lineString
    	 * @param {Array<Array<number>>} coordinates an array of Positions
    	 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
    	 * @param {Object} [options={}] Optional Parameters
    	 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
    	 * @param {string|number} [options.id] Identifier associated with the Feature
    	 * @returns {Feature<LineString>} LineString Feature
    	 * @example
    	 * var linestring1 = turf.lineString([[-24, 63], [-23, 60], [-25, 65], [-20, 69]], {name: 'line 1'});
    	 * var linestring2 = turf.lineString([[-14, 43], [-13, 40], [-15, 45], [-10, 49]], {name: 'line 2'});
    	 *
    	 * //=linestring1
    	 * //=linestring2
    	 */
    	function lineString(coordinates, properties, options) {
    	    if (options === void 0) { options = {}; }
    	    if (coordinates.length < 2) {
    	        throw new Error("coordinates must be an array of two or more positions");
    	    }
    	    var geom = {
    	        type: "LineString",
    	        coordinates: coordinates,
    	    };
    	    return feature(geom, properties, options);
    	}
    	exports.lineString = lineString;
    	/**
    	 * Creates a {@link LineString} {@link FeatureCollection} from an Array of LineString coordinates.
    	 *
    	 * @name lineStrings
    	 * @param {Array<Array<Array<number>>>} coordinates an array of LinearRings
    	 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
    	 * @param {Object} [options={}] Optional Parameters
    	 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north]
    	 * associated with the FeatureCollection
    	 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
    	 * @returns {FeatureCollection<LineString>} LineString FeatureCollection
    	 * @example
    	 * var linestrings = turf.lineStrings([
    	 *   [[-24, 63], [-23, 60], [-25, 65], [-20, 69]],
    	 *   [[-14, 43], [-13, 40], [-15, 45], [-10, 49]]
    	 * ]);
    	 *
    	 * //=linestrings
    	 */
    	function lineStrings(coordinates, properties, options) {
    	    if (options === void 0) { options = {}; }
    	    return featureCollection(coordinates.map(function (coords) {
    	        return lineString(coords, properties);
    	    }), options);
    	}
    	exports.lineStrings = lineStrings;
    	/**
    	 * Takes one or more {@link Feature|Features} and creates a {@link FeatureCollection}.
    	 *
    	 * @name featureCollection
    	 * @param {Feature[]} features input features
    	 * @param {Object} [options={}] Optional Parameters
    	 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
    	 * @param {string|number} [options.id] Identifier associated with the Feature
    	 * @returns {FeatureCollection} FeatureCollection of Features
    	 * @example
    	 * var locationA = turf.point([-75.343, 39.984], {name: 'Location A'});
    	 * var locationB = turf.point([-75.833, 39.284], {name: 'Location B'});
    	 * var locationC = turf.point([-75.534, 39.123], {name: 'Location C'});
    	 *
    	 * var collection = turf.featureCollection([
    	 *   locationA,
    	 *   locationB,
    	 *   locationC
    	 * ]);
    	 *
    	 * //=collection
    	 */
    	function featureCollection(features, options) {
    	    if (options === void 0) { options = {}; }
    	    var fc = { type: "FeatureCollection" };
    	    if (options.id) {
    	        fc.id = options.id;
    	    }
    	    if (options.bbox) {
    	        fc.bbox = options.bbox;
    	    }
    	    fc.features = features;
    	    return fc;
    	}
    	exports.featureCollection = featureCollection;
    	/**
    	 * Creates a {@link Feature<MultiLineString>} based on a
    	 * coordinate array. Properties can be added optionally.
    	 *
    	 * @name multiLineString
    	 * @param {Array<Array<Array<number>>>} coordinates an array of LineStrings
    	 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
    	 * @param {Object} [options={}] Optional Parameters
    	 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
    	 * @param {string|number} [options.id] Identifier associated with the Feature
    	 * @returns {Feature<MultiLineString>} a MultiLineString feature
    	 * @throws {Error} if no coordinates are passed
    	 * @example
    	 * var multiLine = turf.multiLineString([[[0,0],[10,10]]]);
    	 *
    	 * //=multiLine
    	 */
    	function multiLineString(coordinates, properties, options) {
    	    if (options === void 0) { options = {}; }
    	    var geom = {
    	        type: "MultiLineString",
    	        coordinates: coordinates,
    	    };
    	    return feature(geom, properties, options);
    	}
    	exports.multiLineString = multiLineString;
    	/**
    	 * Creates a {@link Feature<MultiPoint>} based on a
    	 * coordinate array. Properties can be added optionally.
    	 *
    	 * @name multiPoint
    	 * @param {Array<Array<number>>} coordinates an array of Positions
    	 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
    	 * @param {Object} [options={}] Optional Parameters
    	 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
    	 * @param {string|number} [options.id] Identifier associated with the Feature
    	 * @returns {Feature<MultiPoint>} a MultiPoint feature
    	 * @throws {Error} if no coordinates are passed
    	 * @example
    	 * var multiPt = turf.multiPoint([[0,0],[10,10]]);
    	 *
    	 * //=multiPt
    	 */
    	function multiPoint(coordinates, properties, options) {
    	    if (options === void 0) { options = {}; }
    	    var geom = {
    	        type: "MultiPoint",
    	        coordinates: coordinates,
    	    };
    	    return feature(geom, properties, options);
    	}
    	exports.multiPoint = multiPoint;
    	/**
    	 * Creates a {@link Feature<MultiPolygon>} based on a
    	 * coordinate array. Properties can be added optionally.
    	 *
    	 * @name multiPolygon
    	 * @param {Array<Array<Array<Array<number>>>>} coordinates an array of Polygons
    	 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
    	 * @param {Object} [options={}] Optional Parameters
    	 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
    	 * @param {string|number} [options.id] Identifier associated with the Feature
    	 * @returns {Feature<MultiPolygon>} a multipolygon feature
    	 * @throws {Error} if no coordinates are passed
    	 * @example
    	 * var multiPoly = turf.multiPolygon([[[[0,0],[0,10],[10,10],[10,0],[0,0]]]]);
    	 *
    	 * //=multiPoly
    	 *
    	 */
    	function multiPolygon(coordinates, properties, options) {
    	    if (options === void 0) { options = {}; }
    	    var geom = {
    	        type: "MultiPolygon",
    	        coordinates: coordinates,
    	    };
    	    return feature(geom, properties, options);
    	}
    	exports.multiPolygon = multiPolygon;
    	/**
    	 * Creates a {@link Feature<GeometryCollection>} based on a
    	 * coordinate array. Properties can be added optionally.
    	 *
    	 * @name geometryCollection
    	 * @param {Array<Geometry>} geometries an array of GeoJSON Geometries
    	 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
    	 * @param {Object} [options={}] Optional Parameters
    	 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
    	 * @param {string|number} [options.id] Identifier associated with the Feature
    	 * @returns {Feature<GeometryCollection>} a GeoJSON GeometryCollection Feature
    	 * @example
    	 * var pt = turf.geometry("Point", [100, 0]);
    	 * var line = turf.geometry("LineString", [[101, 0], [102, 1]]);
    	 * var collection = turf.geometryCollection([pt, line]);
    	 *
    	 * // => collection
    	 */
    	function geometryCollection(geometries, properties, options) {
    	    if (options === void 0) { options = {}; }
    	    var geom = {
    	        type: "GeometryCollection",
    	        geometries: geometries,
    	    };
    	    return feature(geom, properties, options);
    	}
    	exports.geometryCollection = geometryCollection;
    	/**
    	 * Round number to precision
    	 *
    	 * @param {number} num Number
    	 * @param {number} [precision=0] Precision
    	 * @returns {number} rounded number
    	 * @example
    	 * turf.round(120.4321)
    	 * //=120
    	 *
    	 * turf.round(120.4321, 2)
    	 * //=120.43
    	 */
    	function round(num, precision) {
    	    if (precision === void 0) { precision = 0; }
    	    if (precision && !(precision >= 0)) {
    	        throw new Error("precision must be a positive number");
    	    }
    	    var multiplier = Math.pow(10, precision || 0);
    	    return Math.round(num * multiplier) / multiplier;
    	}
    	exports.round = round;
    	/**
    	 * Convert a distance measurement (assuming a spherical Earth) from radians to a more friendly unit.
    	 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
    	 *
    	 * @name radiansToLength
    	 * @param {number} radians in radians across the sphere
    	 * @param {string} [units="kilometers"] can be degrees, radians, miles, inches, yards, metres,
    	 * meters, kilometres, kilometers.
    	 * @returns {number} distance
    	 */
    	function radiansToLength(radians, units) {
    	    if (units === void 0) { units = "kilometers"; }
    	    var factor = exports.factors[units];
    	    if (!factor) {
    	        throw new Error(units + " units is invalid");
    	    }
    	    return radians * factor;
    	}
    	exports.radiansToLength = radiansToLength;
    	/**
    	 * Convert a distance measurement (assuming a spherical Earth) from a real-world unit into radians
    	 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
    	 *
    	 * @name lengthToRadians
    	 * @param {number} distance in real units
    	 * @param {string} [units="kilometers"] can be degrees, radians, miles, inches, yards, metres,
    	 * meters, kilometres, kilometers.
    	 * @returns {number} radians
    	 */
    	function lengthToRadians(distance, units) {
    	    if (units === void 0) { units = "kilometers"; }
    	    var factor = exports.factors[units];
    	    if (!factor) {
    	        throw new Error(units + " units is invalid");
    	    }
    	    return distance / factor;
    	}
    	exports.lengthToRadians = lengthToRadians;
    	/**
    	 * Convert a distance measurement (assuming a spherical Earth) from a real-world unit into degrees
    	 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, centimeters, kilometres, feet
    	 *
    	 * @name lengthToDegrees
    	 * @param {number} distance in real units
    	 * @param {string} [units="kilometers"] can be degrees, radians, miles, inches, yards, metres,
    	 * meters, kilometres, kilometers.
    	 * @returns {number} degrees
    	 */
    	function lengthToDegrees(distance, units) {
    	    return radiansToDegrees(lengthToRadians(distance, units));
    	}
    	exports.lengthToDegrees = lengthToDegrees;
    	/**
    	 * Converts any bearing angle from the north line direction (positive clockwise)
    	 * and returns an angle between 0-360 degrees (positive clockwise), 0 being the north line
    	 *
    	 * @name bearingToAzimuth
    	 * @param {number} bearing angle, between -180 and +180 degrees
    	 * @returns {number} angle between 0 and 360 degrees
    	 */
    	function bearingToAzimuth(bearing) {
    	    var angle = bearing % 360;
    	    if (angle < 0) {
    	        angle += 360;
    	    }
    	    return angle;
    	}
    	exports.bearingToAzimuth = bearingToAzimuth;
    	/**
    	 * Converts an angle in radians to degrees
    	 *
    	 * @name radiansToDegrees
    	 * @param {number} radians angle in radians
    	 * @returns {number} degrees between 0 and 360 degrees
    	 */
    	function radiansToDegrees(radians) {
    	    var degrees = radians % (2 * Math.PI);
    	    return (degrees * 180) / Math.PI;
    	}
    	exports.radiansToDegrees = radiansToDegrees;
    	/**
    	 * Converts an angle in degrees to radians
    	 *
    	 * @name degreesToRadians
    	 * @param {number} degrees angle between 0 and 360 degrees
    	 * @returns {number} angle in radians
    	 */
    	function degreesToRadians(degrees) {
    	    var radians = degrees % 360;
    	    return (radians * Math.PI) / 180;
    	}
    	exports.degreesToRadians = degreesToRadians;
    	/**
    	 * Converts a length to the requested unit.
    	 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
    	 *
    	 * @param {number} length to be converted
    	 * @param {Units} [originalUnit="kilometers"] of the length
    	 * @param {Units} [finalUnit="kilometers"] returned unit
    	 * @returns {number} the converted length
    	 */
    	function convertLength(length, originalUnit, finalUnit) {
    	    if (originalUnit === void 0) { originalUnit = "kilometers"; }
    	    if (finalUnit === void 0) { finalUnit = "kilometers"; }
    	    if (!(length >= 0)) {
    	        throw new Error("length must be a positive number");
    	    }
    	    return radiansToLength(lengthToRadians(length, originalUnit), finalUnit);
    	}
    	exports.convertLength = convertLength;
    	/**
    	 * Converts a area to the requested unit.
    	 * Valid units: kilometers, kilometres, meters, metres, centimetres, millimeters, acres, miles, yards, feet, inches, hectares
    	 * @param {number} area to be converted
    	 * @param {Units} [originalUnit="meters"] of the distance
    	 * @param {Units} [finalUnit="kilometers"] returned unit
    	 * @returns {number} the converted area
    	 */
    	function convertArea(area, originalUnit, finalUnit) {
    	    if (originalUnit === void 0) { originalUnit = "meters"; }
    	    if (finalUnit === void 0) { finalUnit = "kilometers"; }
    	    if (!(area >= 0)) {
    	        throw new Error("area must be a positive number");
    	    }
    	    var startFactor = exports.areaFactors[originalUnit];
    	    if (!startFactor) {
    	        throw new Error("invalid original units");
    	    }
    	    var finalFactor = exports.areaFactors[finalUnit];
    	    if (!finalFactor) {
    	        throw new Error("invalid final units");
    	    }
    	    return (area / startFactor) * finalFactor;
    	}
    	exports.convertArea = convertArea;
    	/**
    	 * isNumber
    	 *
    	 * @param {*} num Number to validate
    	 * @returns {boolean} true/false
    	 * @example
    	 * turf.isNumber(123)
    	 * //=true
    	 * turf.isNumber('foo')
    	 * //=false
    	 */
    	function isNumber(num) {
    	    return !isNaN(num) && num !== null && !Array.isArray(num);
    	}
    	exports.isNumber = isNumber;
    	/**
    	 * isObject
    	 *
    	 * @param {*} input variable to validate
    	 * @returns {boolean} true/false
    	 * @example
    	 * turf.isObject({elevation: 10})
    	 * //=true
    	 * turf.isObject('foo')
    	 * //=false
    	 */
    	function isObject(input) {
    	    return !!input && input.constructor === Object;
    	}
    	exports.isObject = isObject;
    	/**
    	 * Validate BBox
    	 *
    	 * @private
    	 * @param {Array<number>} bbox BBox to validate
    	 * @returns {void}
    	 * @throws Error if BBox is not valid
    	 * @example
    	 * validateBBox([-180, -40, 110, 50])
    	 * //=OK
    	 * validateBBox([-180, -40])
    	 * //=Error
    	 * validateBBox('Foo')
    	 * //=Error
    	 * validateBBox(5)
    	 * //=Error
    	 * validateBBox(null)
    	 * //=Error
    	 * validateBBox(undefined)
    	 * //=Error
    	 */
    	function validateBBox(bbox) {
    	    if (!bbox) {
    	        throw new Error("bbox is required");
    	    }
    	    if (!Array.isArray(bbox)) {
    	        throw new Error("bbox must be an Array");
    	    }
    	    if (bbox.length !== 4 && bbox.length !== 6) {
    	        throw new Error("bbox must be an Array of 4 or 6 numbers");
    	    }
    	    bbox.forEach(function (num) {
    	        if (!isNumber(num)) {
    	            throw new Error("bbox must only contain numbers");
    	        }
    	    });
    	}
    	exports.validateBBox = validateBBox;
    	/**
    	 * Validate Id
    	 *
    	 * @private
    	 * @param {string|number} id Id to validate
    	 * @returns {void}
    	 * @throws Error if Id is not valid
    	 * @example
    	 * validateId([-180, -40, 110, 50])
    	 * //=Error
    	 * validateId([-180, -40])
    	 * //=Error
    	 * validateId('Foo')
    	 * //=OK
    	 * validateId(5)
    	 * //=OK
    	 * validateId(null)
    	 * //=Error
    	 * validateId(undefined)
    	 * //=Error
    	 */
    	function validateId(id) {
    	    if (!id) {
    	        throw new Error("id is required");
    	    }
    	    if (["string", "number"].indexOf(typeof id) === -1) {
    	        throw new Error("id must be a number or a string");
    	    }
    	}
    	exports.validateId = validateId; 
    } (js$1));

    var js = {};

    Object.defineProperty(js, '__esModule', { value: true });

    var helpers$1 = js$1;

    /**
     * Callback for coordEach
     *
     * @callback coordEachCallback
     * @param {Array<number>} currentCoord The current coordinate being processed.
     * @param {number} coordIndex The current index of the coordinate being processed.
     * @param {number} featureIndex The current index of the Feature being processed.
     * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
     * @param {number} geometryIndex The current index of the Geometry being processed.
     */

    /**
     * Iterate over coordinates in any GeoJSON object, similar to Array.forEach()
     *
     * @name coordEach
     * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
     * @param {Function} callback a method that takes (currentCoord, coordIndex, featureIndex, multiFeatureIndex)
     * @param {boolean} [excludeWrapCoord=false] whether or not to include the final coordinate of LinearRings that wraps the ring in its iteration.
     * @returns {void}
     * @example
     * var features = turf.featureCollection([
     *   turf.point([26, 37], {"foo": "bar"}),
     *   turf.point([36, 53], {"hello": "world"})
     * ]);
     *
     * turf.coordEach(features, function (currentCoord, coordIndex, featureIndex, multiFeatureIndex, geometryIndex) {
     *   //=currentCoord
     *   //=coordIndex
     *   //=featureIndex
     *   //=multiFeatureIndex
     *   //=geometryIndex
     * });
     */
    function coordEach$1(geojson, callback, excludeWrapCoord) {
      // Handles null Geometry -- Skips this GeoJSON
      if (geojson === null) return;
      var j,
        k,
        l,
        geometry,
        stopG,
        coords,
        geometryMaybeCollection,
        wrapShrink = 0,
        coordIndex = 0,
        isGeometryCollection,
        type = geojson.type,
        isFeatureCollection = type === "FeatureCollection",
        isFeature = type === "Feature",
        stop = isFeatureCollection ? geojson.features.length : 1;

      // This logic may look a little weird. The reason why it is that way
      // is because it's trying to be fast. GeoJSON supports multiple kinds
      // of objects at its root: FeatureCollection, Features, Geometries.
      // This function has the responsibility of handling all of them, and that
      // means that some of the `for` loops you see below actually just don't apply
      // to certain inputs. For instance, if you give this just a
      // Point geometry, then both loops are short-circuited and all we do
      // is gradually rename the input until it's called 'geometry'.
      //
      // This also aims to allocate as few resources as possible: just a
      // few numbers and booleans, rather than any temporary arrays as would
      // be required with the normalization approach.
      for (var featureIndex = 0; featureIndex < stop; featureIndex++) {
        geometryMaybeCollection = isFeatureCollection
          ? geojson.features[featureIndex].geometry
          : isFeature
          ? geojson.geometry
          : geojson;
        isGeometryCollection = geometryMaybeCollection
          ? geometryMaybeCollection.type === "GeometryCollection"
          : false;
        stopG = isGeometryCollection
          ? geometryMaybeCollection.geometries.length
          : 1;

        for (var geomIndex = 0; geomIndex < stopG; geomIndex++) {
          var multiFeatureIndex = 0;
          var geometryIndex = 0;
          geometry = isGeometryCollection
            ? geometryMaybeCollection.geometries[geomIndex]
            : geometryMaybeCollection;

          // Handles null Geometry -- Skips this geometry
          if (geometry === null) continue;
          coords = geometry.coordinates;
          var geomType = geometry.type;

          wrapShrink =
            excludeWrapCoord &&
            (geomType === "Polygon" || geomType === "MultiPolygon")
              ? 1
              : 0;

          switch (geomType) {
            case null:
              break;
            case "Point":
              if (
                callback(
                  coords,
                  coordIndex,
                  featureIndex,
                  multiFeatureIndex,
                  geometryIndex
                ) === false
              )
                return false;
              coordIndex++;
              multiFeatureIndex++;
              break;
            case "LineString":
            case "MultiPoint":
              for (j = 0; j < coords.length; j++) {
                if (
                  callback(
                    coords[j],
                    coordIndex,
                    featureIndex,
                    multiFeatureIndex,
                    geometryIndex
                  ) === false
                )
                  return false;
                coordIndex++;
                if (geomType === "MultiPoint") multiFeatureIndex++;
              }
              if (geomType === "LineString") multiFeatureIndex++;
              break;
            case "Polygon":
            case "MultiLineString":
              for (j = 0; j < coords.length; j++) {
                for (k = 0; k < coords[j].length - wrapShrink; k++) {
                  if (
                    callback(
                      coords[j][k],
                      coordIndex,
                      featureIndex,
                      multiFeatureIndex,
                      geometryIndex
                    ) === false
                  )
                    return false;
                  coordIndex++;
                }
                if (geomType === "MultiLineString") multiFeatureIndex++;
                if (geomType === "Polygon") geometryIndex++;
              }
              if (geomType === "Polygon") multiFeatureIndex++;
              break;
            case "MultiPolygon":
              for (j = 0; j < coords.length; j++) {
                geometryIndex = 0;
                for (k = 0; k < coords[j].length; k++) {
                  for (l = 0; l < coords[j][k].length - wrapShrink; l++) {
                    if (
                      callback(
                        coords[j][k][l],
                        coordIndex,
                        featureIndex,
                        multiFeatureIndex,
                        geometryIndex
                      ) === false
                    )
                      return false;
                    coordIndex++;
                  }
                  geometryIndex++;
                }
                multiFeatureIndex++;
              }
              break;
            case "GeometryCollection":
              for (j = 0; j < geometry.geometries.length; j++)
                if (
                  coordEach$1(geometry.geometries[j], callback, excludeWrapCoord) ===
                  false
                )
                  return false;
              break;
            default:
              throw new Error("Unknown Geometry Type");
          }
        }
      }
    }

    /**
     * Callback for coordReduce
     *
     * The first time the callback function is called, the values provided as arguments depend
     * on whether the reduce method has an initialValue argument.
     *
     * If an initialValue is provided to the reduce method:
     *  - The previousValue argument is initialValue.
     *  - The currentValue argument is the value of the first element present in the array.
     *
     * If an initialValue is not provided:
     *  - The previousValue argument is the value of the first element present in the array.
     *  - The currentValue argument is the value of the second element present in the array.
     *
     * @callback coordReduceCallback
     * @param {*} previousValue The accumulated value previously returned in the last invocation
     * of the callback, or initialValue, if supplied.
     * @param {Array<number>} currentCoord The current coordinate being processed.
     * @param {number} coordIndex The current index of the coordinate being processed.
     * Starts at index 0, if an initialValue is provided, and at index 1 otherwise.
     * @param {number} featureIndex The current index of the Feature being processed.
     * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
     * @param {number} geometryIndex The current index of the Geometry being processed.
     */

    /**
     * Reduce coordinates in any GeoJSON object, similar to Array.reduce()
     *
     * @name coordReduce
     * @param {FeatureCollection|Geometry|Feature} geojson any GeoJSON object
     * @param {Function} callback a method that takes (previousValue, currentCoord, coordIndex)
     * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
     * @param {boolean} [excludeWrapCoord=false] whether or not to include the final coordinate of LinearRings that wraps the ring in its iteration.
     * @returns {*} The value that results from the reduction.
     * @example
     * var features = turf.featureCollection([
     *   turf.point([26, 37], {"foo": "bar"}),
     *   turf.point([36, 53], {"hello": "world"})
     * ]);
     *
     * turf.coordReduce(features, function (previousValue, currentCoord, coordIndex, featureIndex, multiFeatureIndex, geometryIndex) {
     *   //=previousValue
     *   //=currentCoord
     *   //=coordIndex
     *   //=featureIndex
     *   //=multiFeatureIndex
     *   //=geometryIndex
     *   return currentCoord;
     * });
     */
    function coordReduce$1(geojson, callback, initialValue, excludeWrapCoord) {
      var previousValue = initialValue;
      coordEach$1(
        geojson,
        function (
          currentCoord,
          coordIndex,
          featureIndex,
          multiFeatureIndex,
          geometryIndex
        ) {
          if (coordIndex === 0 && initialValue === undefined)
            previousValue = currentCoord;
          else
            previousValue = callback(
              previousValue,
              currentCoord,
              coordIndex,
              featureIndex,
              multiFeatureIndex,
              geometryIndex
            );
        },
        excludeWrapCoord
      );
      return previousValue;
    }

    /**
     * Callback for propEach
     *
     * @callback propEachCallback
     * @param {Object} currentProperties The current Properties being processed.
     * @param {number} featureIndex The current index of the Feature being processed.
     */

    /**
     * Iterate over properties in any GeoJSON object, similar to Array.forEach()
     *
     * @name propEach
     * @param {FeatureCollection|Feature} geojson any GeoJSON object
     * @param {Function} callback a method that takes (currentProperties, featureIndex)
     * @returns {void}
     * @example
     * var features = turf.featureCollection([
     *     turf.point([26, 37], {foo: 'bar'}),
     *     turf.point([36, 53], {hello: 'world'})
     * ]);
     *
     * turf.propEach(features, function (currentProperties, featureIndex) {
     *   //=currentProperties
     *   //=featureIndex
     * });
     */
    function propEach$1(geojson, callback) {
      var i;
      switch (geojson.type) {
        case "FeatureCollection":
          for (i = 0; i < geojson.features.length; i++) {
            if (callback(geojson.features[i].properties, i) === false) break;
          }
          break;
        case "Feature":
          callback(geojson.properties, 0);
          break;
      }
    }

    /**
     * Callback for propReduce
     *
     * The first time the callback function is called, the values provided as arguments depend
     * on whether the reduce method has an initialValue argument.
     *
     * If an initialValue is provided to the reduce method:
     *  - The previousValue argument is initialValue.
     *  - The currentValue argument is the value of the first element present in the array.
     *
     * If an initialValue is not provided:
     *  - The previousValue argument is the value of the first element present in the array.
     *  - The currentValue argument is the value of the second element present in the array.
     *
     * @callback propReduceCallback
     * @param {*} previousValue The accumulated value previously returned in the last invocation
     * of the callback, or initialValue, if supplied.
     * @param {*} currentProperties The current Properties being processed.
     * @param {number} featureIndex The current index of the Feature being processed.
     */

    /**
     * Reduce properties in any GeoJSON object into a single value,
     * similar to how Array.reduce works. However, in this case we lazily run
     * the reduction, so an array of all properties is unnecessary.
     *
     * @name propReduce
     * @param {FeatureCollection|Feature} geojson any GeoJSON object
     * @param {Function} callback a method that takes (previousValue, currentProperties, featureIndex)
     * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
     * @returns {*} The value that results from the reduction.
     * @example
     * var features = turf.featureCollection([
     *     turf.point([26, 37], {foo: 'bar'}),
     *     turf.point([36, 53], {hello: 'world'})
     * ]);
     *
     * turf.propReduce(features, function (previousValue, currentProperties, featureIndex) {
     *   //=previousValue
     *   //=currentProperties
     *   //=featureIndex
     *   return currentProperties
     * });
     */
    function propReduce$1(geojson, callback, initialValue) {
      var previousValue = initialValue;
      propEach$1(geojson, function (currentProperties, featureIndex) {
        if (featureIndex === 0 && initialValue === undefined)
          previousValue = currentProperties;
        else
          previousValue = callback(previousValue, currentProperties, featureIndex);
      });
      return previousValue;
    }

    /**
     * Callback for featureEach
     *
     * @callback featureEachCallback
     * @param {Feature<any>} currentFeature The current Feature being processed.
     * @param {number} featureIndex The current index of the Feature being processed.
     */

    /**
     * Iterate over features in any GeoJSON object, similar to
     * Array.forEach.
     *
     * @name featureEach
     * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
     * @param {Function} callback a method that takes (currentFeature, featureIndex)
     * @returns {void}
     * @example
     * var features = turf.featureCollection([
     *   turf.point([26, 37], {foo: 'bar'}),
     *   turf.point([36, 53], {hello: 'world'})
     * ]);
     *
     * turf.featureEach(features, function (currentFeature, featureIndex) {
     *   //=currentFeature
     *   //=featureIndex
     * });
     */
    function featureEach$2(geojson, callback) {
      if (geojson.type === "Feature") {
        callback(geojson, 0);
      } else if (geojson.type === "FeatureCollection") {
        for (var i = 0; i < geojson.features.length; i++) {
          if (callback(geojson.features[i], i) === false) break;
        }
      }
    }

    /**
     * Callback for featureReduce
     *
     * The first time the callback function is called, the values provided as arguments depend
     * on whether the reduce method has an initialValue argument.
     *
     * If an initialValue is provided to the reduce method:
     *  - The previousValue argument is initialValue.
     *  - The currentValue argument is the value of the first element present in the array.
     *
     * If an initialValue is not provided:
     *  - The previousValue argument is the value of the first element present in the array.
     *  - The currentValue argument is the value of the second element present in the array.
     *
     * @callback featureReduceCallback
     * @param {*} previousValue The accumulated value previously returned in the last invocation
     * of the callback, or initialValue, if supplied.
     * @param {Feature} currentFeature The current Feature being processed.
     * @param {number} featureIndex The current index of the Feature being processed.
     */

    /**
     * Reduce features in any GeoJSON object, similar to Array.reduce().
     *
     * @name featureReduce
     * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
     * @param {Function} callback a method that takes (previousValue, currentFeature, featureIndex)
     * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
     * @returns {*} The value that results from the reduction.
     * @example
     * var features = turf.featureCollection([
     *   turf.point([26, 37], {"foo": "bar"}),
     *   turf.point([36, 53], {"hello": "world"})
     * ]);
     *
     * turf.featureReduce(features, function (previousValue, currentFeature, featureIndex) {
     *   //=previousValue
     *   //=currentFeature
     *   //=featureIndex
     *   return currentFeature
     * });
     */
    function featureReduce$1(geojson, callback, initialValue) {
      var previousValue = initialValue;
      featureEach$2(geojson, function (currentFeature, featureIndex) {
        if (featureIndex === 0 && initialValue === undefined)
          previousValue = currentFeature;
        else previousValue = callback(previousValue, currentFeature, featureIndex);
      });
      return previousValue;
    }

    /**
     * Get all coordinates from any GeoJSON object.
     *
     * @name coordAll
     * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
     * @returns {Array<Array<number>>} coordinate position array
     * @example
     * var features = turf.featureCollection([
     *   turf.point([26, 37], {foo: 'bar'}),
     *   turf.point([36, 53], {hello: 'world'})
     * ]);
     *
     * var coords = turf.coordAll(features);
     * //= [[26, 37], [36, 53]]
     */
    function coordAll$1(geojson) {
      var coords = [];
      coordEach$1(geojson, function (coord) {
        coords.push(coord);
      });
      return coords;
    }

    /**
     * Callback for geomEach
     *
     * @callback geomEachCallback
     * @param {Geometry} currentGeometry The current Geometry being processed.
     * @param {number} featureIndex The current index of the Feature being processed.
     * @param {Object} featureProperties The current Feature Properties being processed.
     * @param {Array<number>} featureBBox The current Feature BBox being processed.
     * @param {number|string} featureId The current Feature Id being processed.
     */

    /**
     * Iterate over each geometry in any GeoJSON object, similar to Array.forEach()
     *
     * @name geomEach
     * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
     * @param {Function} callback a method that takes (currentGeometry, featureIndex, featureProperties, featureBBox, featureId)
     * @returns {void}
     * @example
     * var features = turf.featureCollection([
     *     turf.point([26, 37], {foo: 'bar'}),
     *     turf.point([36, 53], {hello: 'world'})
     * ]);
     *
     * turf.geomEach(features, function (currentGeometry, featureIndex, featureProperties, featureBBox, featureId) {
     *   //=currentGeometry
     *   //=featureIndex
     *   //=featureProperties
     *   //=featureBBox
     *   //=featureId
     * });
     */
    function geomEach$1(geojson, callback) {
      var i,
        j,
        g,
        geometry,
        stopG,
        geometryMaybeCollection,
        isGeometryCollection,
        featureProperties,
        featureBBox,
        featureId,
        featureIndex = 0,
        isFeatureCollection = geojson.type === "FeatureCollection",
        isFeature = geojson.type === "Feature",
        stop = isFeatureCollection ? geojson.features.length : 1;

      // This logic may look a little weird. The reason why it is that way
      // is because it's trying to be fast. GeoJSON supports multiple kinds
      // of objects at its root: FeatureCollection, Features, Geometries.
      // This function has the responsibility of handling all of them, and that
      // means that some of the `for` loops you see below actually just don't apply
      // to certain inputs. For instance, if you give this just a
      // Point geometry, then both loops are short-circuited and all we do
      // is gradually rename the input until it's called 'geometry'.
      //
      // This also aims to allocate as few resources as possible: just a
      // few numbers and booleans, rather than any temporary arrays as would
      // be required with the normalization approach.
      for (i = 0; i < stop; i++) {
        geometryMaybeCollection = isFeatureCollection
          ? geojson.features[i].geometry
          : isFeature
          ? geojson.geometry
          : geojson;
        featureProperties = isFeatureCollection
          ? geojson.features[i].properties
          : isFeature
          ? geojson.properties
          : {};
        featureBBox = isFeatureCollection
          ? geojson.features[i].bbox
          : isFeature
          ? geojson.bbox
          : undefined;
        featureId = isFeatureCollection
          ? geojson.features[i].id
          : isFeature
          ? geojson.id
          : undefined;
        isGeometryCollection = geometryMaybeCollection
          ? geometryMaybeCollection.type === "GeometryCollection"
          : false;
        stopG = isGeometryCollection
          ? geometryMaybeCollection.geometries.length
          : 1;

        for (g = 0; g < stopG; g++) {
          geometry = isGeometryCollection
            ? geometryMaybeCollection.geometries[g]
            : geometryMaybeCollection;

          // Handle null Geometry
          if (geometry === null) {
            if (
              callback(
                null,
                featureIndex,
                featureProperties,
                featureBBox,
                featureId
              ) === false
            )
              return false;
            continue;
          }
          switch (geometry.type) {
            case "Point":
            case "LineString":
            case "MultiPoint":
            case "Polygon":
            case "MultiLineString":
            case "MultiPolygon": {
              if (
                callback(
                  geometry,
                  featureIndex,
                  featureProperties,
                  featureBBox,
                  featureId
                ) === false
              )
                return false;
              break;
            }
            case "GeometryCollection": {
              for (j = 0; j < geometry.geometries.length; j++) {
                if (
                  callback(
                    geometry.geometries[j],
                    featureIndex,
                    featureProperties,
                    featureBBox,
                    featureId
                  ) === false
                )
                  return false;
              }
              break;
            }
            default:
              throw new Error("Unknown Geometry Type");
          }
        }
        // Only increase `featureIndex` per each feature
        featureIndex++;
      }
    }

    /**
     * Callback for geomReduce
     *
     * The first time the callback function is called, the values provided as arguments depend
     * on whether the reduce method has an initialValue argument.
     *
     * If an initialValue is provided to the reduce method:
     *  - The previousValue argument is initialValue.
     *  - The currentValue argument is the value of the first element present in the array.
     *
     * If an initialValue is not provided:
     *  - The previousValue argument is the value of the first element present in the array.
     *  - The currentValue argument is the value of the second element present in the array.
     *
     * @callback geomReduceCallback
     * @param {*} previousValue The accumulated value previously returned in the last invocation
     * of the callback, or initialValue, if supplied.
     * @param {Geometry} currentGeometry The current Geometry being processed.
     * @param {number} featureIndex The current index of the Feature being processed.
     * @param {Object} featureProperties The current Feature Properties being processed.
     * @param {Array<number>} featureBBox The current Feature BBox being processed.
     * @param {number|string} featureId The current Feature Id being processed.
     */

    /**
     * Reduce geometry in any GeoJSON object, similar to Array.reduce().
     *
     * @name geomReduce
     * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
     * @param {Function} callback a method that takes (previousValue, currentGeometry, featureIndex, featureProperties, featureBBox, featureId)
     * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
     * @returns {*} The value that results from the reduction.
     * @example
     * var features = turf.featureCollection([
     *     turf.point([26, 37], {foo: 'bar'}),
     *     turf.point([36, 53], {hello: 'world'})
     * ]);
     *
     * turf.geomReduce(features, function (previousValue, currentGeometry, featureIndex, featureProperties, featureBBox, featureId) {
     *   //=previousValue
     *   //=currentGeometry
     *   //=featureIndex
     *   //=featureProperties
     *   //=featureBBox
     *   //=featureId
     *   return currentGeometry
     * });
     */
    function geomReduce$1(geojson, callback, initialValue) {
      var previousValue = initialValue;
      geomEach$1(
        geojson,
        function (
          currentGeometry,
          featureIndex,
          featureProperties,
          featureBBox,
          featureId
        ) {
          if (featureIndex === 0 && initialValue === undefined)
            previousValue = currentGeometry;
          else
            previousValue = callback(
              previousValue,
              currentGeometry,
              featureIndex,
              featureProperties,
              featureBBox,
              featureId
            );
        }
      );
      return previousValue;
    }

    /**
     * Callback for flattenEach
     *
     * @callback flattenEachCallback
     * @param {Feature} currentFeature The current flattened feature being processed.
     * @param {number} featureIndex The current index of the Feature being processed.
     * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
     */

    /**
     * Iterate over flattened features in any GeoJSON object, similar to
     * Array.forEach.
     *
     * @name flattenEach
     * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
     * @param {Function} callback a method that takes (currentFeature, featureIndex, multiFeatureIndex)
     * @example
     * var features = turf.featureCollection([
     *     turf.point([26, 37], {foo: 'bar'}),
     *     turf.multiPoint([[40, 30], [36, 53]], {hello: 'world'})
     * ]);
     *
     * turf.flattenEach(features, function (currentFeature, featureIndex, multiFeatureIndex) {
     *   //=currentFeature
     *   //=featureIndex
     *   //=multiFeatureIndex
     * });
     */
    function flattenEach$1(geojson, callback) {
      geomEach$1(geojson, function (geometry, featureIndex, properties, bbox, id) {
        // Callback for single geometry
        var type = geometry === null ? null : geometry.type;
        switch (type) {
          case null:
          case "Point":
          case "LineString":
          case "Polygon":
            if (
              callback(
                helpers$1.feature(geometry, properties, { bbox: bbox, id: id }),
                featureIndex,
                0
              ) === false
            )
              return false;
            return;
        }

        var geomType;

        // Callback for multi-geometry
        switch (type) {
          case "MultiPoint":
            geomType = "Point";
            break;
          case "MultiLineString":
            geomType = "LineString";
            break;
          case "MultiPolygon":
            geomType = "Polygon";
            break;
        }

        for (
          var multiFeatureIndex = 0;
          multiFeatureIndex < geometry.coordinates.length;
          multiFeatureIndex++
        ) {
          var coordinate = geometry.coordinates[multiFeatureIndex];
          var geom = {
            type: geomType,
            coordinates: coordinate,
          };
          if (
            callback(helpers$1.feature(geom, properties), featureIndex, multiFeatureIndex) ===
            false
          )
            return false;
        }
      });
    }

    /**
     * Callback for flattenReduce
     *
     * The first time the callback function is called, the values provided as arguments depend
     * on whether the reduce method has an initialValue argument.
     *
     * If an initialValue is provided to the reduce method:
     *  - The previousValue argument is initialValue.
     *  - The currentValue argument is the value of the first element present in the array.
     *
     * If an initialValue is not provided:
     *  - The previousValue argument is the value of the first element present in the array.
     *  - The currentValue argument is the value of the second element present in the array.
     *
     * @callback flattenReduceCallback
     * @param {*} previousValue The accumulated value previously returned in the last invocation
     * of the callback, or initialValue, if supplied.
     * @param {Feature} currentFeature The current Feature being processed.
     * @param {number} featureIndex The current index of the Feature being processed.
     * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
     */

    /**
     * Reduce flattened features in any GeoJSON object, similar to Array.reduce().
     *
     * @name flattenReduce
     * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON object
     * @param {Function} callback a method that takes (previousValue, currentFeature, featureIndex, multiFeatureIndex)
     * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
     * @returns {*} The value that results from the reduction.
     * @example
     * var features = turf.featureCollection([
     *     turf.point([26, 37], {foo: 'bar'}),
     *     turf.multiPoint([[40, 30], [36, 53]], {hello: 'world'})
     * ]);
     *
     * turf.flattenReduce(features, function (previousValue, currentFeature, featureIndex, multiFeatureIndex) {
     *   //=previousValue
     *   //=currentFeature
     *   //=featureIndex
     *   //=multiFeatureIndex
     *   return currentFeature
     * });
     */
    function flattenReduce$1(geojson, callback, initialValue) {
      var previousValue = initialValue;
      flattenEach$1(
        geojson,
        function (currentFeature, featureIndex, multiFeatureIndex) {
          if (
            featureIndex === 0 &&
            multiFeatureIndex === 0 &&
            initialValue === undefined
          )
            previousValue = currentFeature;
          else
            previousValue = callback(
              previousValue,
              currentFeature,
              featureIndex,
              multiFeatureIndex
            );
        }
      );
      return previousValue;
    }

    /**
     * Callback for segmentEach
     *
     * @callback segmentEachCallback
     * @param {Feature<LineString>} currentSegment The current Segment being processed.
     * @param {number} featureIndex The current index of the Feature being processed.
     * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
     * @param {number} geometryIndex The current index of the Geometry being processed.
     * @param {number} segmentIndex The current index of the Segment being processed.
     * @returns {void}
     */

    /**
     * Iterate over 2-vertex line segment in any GeoJSON object, similar to Array.forEach()
     * (Multi)Point geometries do not contain segments therefore they are ignored during this operation.
     *
     * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON
     * @param {Function} callback a method that takes (currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex)
     * @returns {void}
     * @example
     * var polygon = turf.polygon([[[-50, 5], [-40, -10], [-50, -10], [-40, 5], [-50, 5]]]);
     *
     * // Iterate over GeoJSON by 2-vertex segments
     * turf.segmentEach(polygon, function (currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex) {
     *   //=currentSegment
     *   //=featureIndex
     *   //=multiFeatureIndex
     *   //=geometryIndex
     *   //=segmentIndex
     * });
     *
     * // Calculate the total number of segments
     * var total = 0;
     * turf.segmentEach(polygon, function () {
     *     total++;
     * });
     */
    function segmentEach$1(geojson, callback) {
      flattenEach$1(geojson, function (feature, featureIndex, multiFeatureIndex) {
        var segmentIndex = 0;

        // Exclude null Geometries
        if (!feature.geometry) return;
        // (Multi)Point geometries do not contain segments therefore they are ignored during this operation.
        var type = feature.geometry.type;
        if (type === "Point" || type === "MultiPoint") return;

        // Generate 2-vertex line segments
        var previousCoords;
        var previousFeatureIndex = 0;
        var previousMultiIndex = 0;
        var prevGeomIndex = 0;
        if (
          coordEach$1(
            feature,
            function (
              currentCoord,
              coordIndex,
              featureIndexCoord,
              multiPartIndexCoord,
              geometryIndex
            ) {
              // Simulating a meta.coordReduce() since `reduce` operations cannot be stopped by returning `false`
              if (
                previousCoords === undefined ||
                featureIndex > previousFeatureIndex ||
                multiPartIndexCoord > previousMultiIndex ||
                geometryIndex > prevGeomIndex
              ) {
                previousCoords = currentCoord;
                previousFeatureIndex = featureIndex;
                previousMultiIndex = multiPartIndexCoord;
                prevGeomIndex = geometryIndex;
                segmentIndex = 0;
                return;
              }
              var currentSegment = helpers$1.lineString(
                [previousCoords, currentCoord],
                feature.properties
              );
              if (
                callback(
                  currentSegment,
                  featureIndex,
                  multiFeatureIndex,
                  geometryIndex,
                  segmentIndex
                ) === false
              )
                return false;
              segmentIndex++;
              previousCoords = currentCoord;
            }
          ) === false
        )
          return false;
      });
    }

    /**
     * Callback for segmentReduce
     *
     * The first time the callback function is called, the values provided as arguments depend
     * on whether the reduce method has an initialValue argument.
     *
     * If an initialValue is provided to the reduce method:
     *  - The previousValue argument is initialValue.
     *  - The currentValue argument is the value of the first element present in the array.
     *
     * If an initialValue is not provided:
     *  - The previousValue argument is the value of the first element present in the array.
     *  - The currentValue argument is the value of the second element present in the array.
     *
     * @callback segmentReduceCallback
     * @param {*} previousValue The accumulated value previously returned in the last invocation
     * of the callback, or initialValue, if supplied.
     * @param {Feature<LineString>} currentSegment The current Segment being processed.
     * @param {number} featureIndex The current index of the Feature being processed.
     * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed.
     * @param {number} geometryIndex The current index of the Geometry being processed.
     * @param {number} segmentIndex The current index of the Segment being processed.
     */

    /**
     * Reduce 2-vertex line segment in any GeoJSON object, similar to Array.reduce()
     * (Multi)Point geometries do not contain segments therefore they are ignored during this operation.
     *
     * @param {FeatureCollection|Feature|Geometry} geojson any GeoJSON
     * @param {Function} callback a method that takes (previousValue, currentSegment, currentIndex)
     * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
     * @returns {void}
     * @example
     * var polygon = turf.polygon([[[-50, 5], [-40, -10], [-50, -10], [-40, 5], [-50, 5]]]);
     *
     * // Iterate over GeoJSON by 2-vertex segments
     * turf.segmentReduce(polygon, function (previousSegment, currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex) {
     *   //= previousSegment
     *   //= currentSegment
     *   //= featureIndex
     *   //= multiFeatureIndex
     *   //= geometryIndex
     *   //= segmentIndex
     *   return currentSegment
     * });
     *
     * // Calculate the total number of segments
     * var initialValue = 0
     * var total = turf.segmentReduce(polygon, function (previousValue) {
     *     previousValue++;
     *     return previousValue;
     * }, initialValue);
     */
    function segmentReduce$1(geojson, callback, initialValue) {
      var previousValue = initialValue;
      var started = false;
      segmentEach$1(
        geojson,
        function (
          currentSegment,
          featureIndex,
          multiFeatureIndex,
          geometryIndex,
          segmentIndex
        ) {
          if (started === false && initialValue === undefined)
            previousValue = currentSegment;
          else
            previousValue = callback(
              previousValue,
              currentSegment,
              featureIndex,
              multiFeatureIndex,
              geometryIndex,
              segmentIndex
            );
          started = true;
        }
      );
      return previousValue;
    }

    /**
     * Callback for lineEach
     *
     * @callback lineEachCallback
     * @param {Feature<LineString>} currentLine The current LineString|LinearRing being processed
     * @param {number} featureIndex The current index of the Feature being processed
     * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed
     * @param {number} geometryIndex The current index of the Geometry being processed
     */

    /**
     * Iterate over line or ring coordinates in LineString, Polygon, MultiLineString, MultiPolygon Features or Geometries,
     * similar to Array.forEach.
     *
     * @name lineEach
     * @param {Geometry|Feature<LineString|Polygon|MultiLineString|MultiPolygon>} geojson object
     * @param {Function} callback a method that takes (currentLine, featureIndex, multiFeatureIndex, geometryIndex)
     * @example
     * var multiLine = turf.multiLineString([
     *   [[26, 37], [35, 45]],
     *   [[36, 53], [38, 50], [41, 55]]
     * ]);
     *
     * turf.lineEach(multiLine, function (currentLine, featureIndex, multiFeatureIndex, geometryIndex) {
     *   //=currentLine
     *   //=featureIndex
     *   //=multiFeatureIndex
     *   //=geometryIndex
     * });
     */
    function lineEach$1(geojson, callback) {
      // validation
      if (!geojson) throw new Error("geojson is required");

      flattenEach$1(geojson, function (feature, featureIndex, multiFeatureIndex) {
        if (feature.geometry === null) return;
        var type = feature.geometry.type;
        var coords = feature.geometry.coordinates;
        switch (type) {
          case "LineString":
            if (callback(feature, featureIndex, multiFeatureIndex, 0, 0) === false)
              return false;
            break;
          case "Polygon":
            for (
              var geometryIndex = 0;
              geometryIndex < coords.length;
              geometryIndex++
            ) {
              if (
                callback(
                  helpers$1.lineString(coords[geometryIndex], feature.properties),
                  featureIndex,
                  multiFeatureIndex,
                  geometryIndex
                ) === false
              )
                return false;
            }
            break;
        }
      });
    }

    /**
     * Callback for lineReduce
     *
     * The first time the callback function is called, the values provided as arguments depend
     * on whether the reduce method has an initialValue argument.
     *
     * If an initialValue is provided to the reduce method:
     *  - The previousValue argument is initialValue.
     *  - The currentValue argument is the value of the first element present in the array.
     *
     * If an initialValue is not provided:
     *  - The previousValue argument is the value of the first element present in the array.
     *  - The currentValue argument is the value of the second element present in the array.
     *
     * @callback lineReduceCallback
     * @param {*} previousValue The accumulated value previously returned in the last invocation
     * of the callback, or initialValue, if supplied.
     * @param {Feature<LineString>} currentLine The current LineString|LinearRing being processed.
     * @param {number} featureIndex The current index of the Feature being processed
     * @param {number} multiFeatureIndex The current index of the Multi-Feature being processed
     * @param {number} geometryIndex The current index of the Geometry being processed
     */

    /**
     * Reduce features in any GeoJSON object, similar to Array.reduce().
     *
     * @name lineReduce
     * @param {Geometry|Feature<LineString|Polygon|MultiLineString|MultiPolygon>} geojson object
     * @param {Function} callback a method that takes (previousValue, currentLine, featureIndex, multiFeatureIndex, geometryIndex)
     * @param {*} [initialValue] Value to use as the first argument to the first call of the callback.
     * @returns {*} The value that results from the reduction.
     * @example
     * var multiPoly = turf.multiPolygon([
     *   turf.polygon([[[12,48],[2,41],[24,38],[12,48]], [[9,44],[13,41],[13,45],[9,44]]]),
     *   turf.polygon([[[5, 5], [0, 0], [2, 2], [4, 4], [5, 5]]])
     * ]);
     *
     * turf.lineReduce(multiPoly, function (previousValue, currentLine, featureIndex, multiFeatureIndex, geometryIndex) {
     *   //=previousValue
     *   //=currentLine
     *   //=featureIndex
     *   //=multiFeatureIndex
     *   //=geometryIndex
     *   return currentLine
     * });
     */
    function lineReduce$1(geojson, callback, initialValue) {
      var previousValue = initialValue;
      lineEach$1(
        geojson,
        function (currentLine, featureIndex, multiFeatureIndex, geometryIndex) {
          if (featureIndex === 0 && initialValue === undefined)
            previousValue = currentLine;
          else
            previousValue = callback(
              previousValue,
              currentLine,
              featureIndex,
              multiFeatureIndex,
              geometryIndex
            );
        }
      );
      return previousValue;
    }

    /**
     * Finds a particular 2-vertex LineString Segment from a GeoJSON using `@turf/meta` indexes.
     *
     * Negative indexes are permitted.
     * Point & MultiPoint will always return null.
     *
     * @param {FeatureCollection|Feature|Geometry} geojson Any GeoJSON Feature or Geometry
     * @param {Object} [options={}] Optional parameters
     * @param {number} [options.featureIndex=0] Feature Index
     * @param {number} [options.multiFeatureIndex=0] Multi-Feature Index
     * @param {number} [options.geometryIndex=0] Geometry Index
     * @param {number} [options.segmentIndex=0] Segment Index
     * @param {Object} [options.properties={}] Translate Properties to output LineString
     * @param {BBox} [options.bbox={}] Translate BBox to output LineString
     * @param {number|string} [options.id={}] Translate Id to output LineString
     * @returns {Feature<LineString>} 2-vertex GeoJSON Feature LineString
     * @example
     * var multiLine = turf.multiLineString([
     *     [[10, 10], [50, 30], [30, 40]],
     *     [[-10, -10], [-50, -30], [-30, -40]]
     * ]);
     *
     * // First Segment (defaults are 0)
     * turf.findSegment(multiLine);
     * // => Feature<LineString<[[10, 10], [50, 30]]>>
     *
     * // First Segment of 2nd Multi Feature
     * turf.findSegment(multiLine, {multiFeatureIndex: 1});
     * // => Feature<LineString<[[-10, -10], [-50, -30]]>>
     *
     * // Last Segment of Last Multi Feature
     * turf.findSegment(multiLine, {multiFeatureIndex: -1, segmentIndex: -1});
     * // => Feature<LineString<[[-50, -30], [-30, -40]]>>
     */
    function findSegment$1(geojson, options) {
      // Optional Parameters
      options = options || {};
      if (!helpers$1.isObject(options)) throw new Error("options is invalid");
      var featureIndex = options.featureIndex || 0;
      var multiFeatureIndex = options.multiFeatureIndex || 0;
      var geometryIndex = options.geometryIndex || 0;
      var segmentIndex = options.segmentIndex || 0;

      // Find FeatureIndex
      var properties = options.properties;
      var geometry;

      switch (geojson.type) {
        case "FeatureCollection":
          if (featureIndex < 0)
            featureIndex = geojson.features.length + featureIndex;
          properties = properties || geojson.features[featureIndex].properties;
          geometry = geojson.features[featureIndex].geometry;
          break;
        case "Feature":
          properties = properties || geojson.properties;
          geometry = geojson.geometry;
          break;
        case "Point":
        case "MultiPoint":
          return null;
        case "LineString":
        case "Polygon":
        case "MultiLineString":
        case "MultiPolygon":
          geometry = geojson;
          break;
        default:
          throw new Error("geojson is invalid");
      }

      // Find SegmentIndex
      if (geometry === null) return null;
      var coords = geometry.coordinates;
      switch (geometry.type) {
        case "Point":
        case "MultiPoint":
          return null;
        case "LineString":
          if (segmentIndex < 0) segmentIndex = coords.length + segmentIndex - 1;
          return helpers$1.lineString(
            [coords[segmentIndex], coords[segmentIndex + 1]],
            properties,
            options
          );
        case "Polygon":
          if (geometryIndex < 0) geometryIndex = coords.length + geometryIndex;
          if (segmentIndex < 0)
            segmentIndex = coords[geometryIndex].length + segmentIndex - 1;
          return helpers$1.lineString(
            [
              coords[geometryIndex][segmentIndex],
              coords[geometryIndex][segmentIndex + 1],
            ],
            properties,
            options
          );
        case "MultiLineString":
          if (multiFeatureIndex < 0)
            multiFeatureIndex = coords.length + multiFeatureIndex;
          if (segmentIndex < 0)
            segmentIndex = coords[multiFeatureIndex].length + segmentIndex - 1;
          return helpers$1.lineString(
            [
              coords[multiFeatureIndex][segmentIndex],
              coords[multiFeatureIndex][segmentIndex + 1],
            ],
            properties,
            options
          );
        case "MultiPolygon":
          if (multiFeatureIndex < 0)
            multiFeatureIndex = coords.length + multiFeatureIndex;
          if (geometryIndex < 0)
            geometryIndex = coords[multiFeatureIndex].length + geometryIndex;
          if (segmentIndex < 0)
            segmentIndex =
              coords[multiFeatureIndex][geometryIndex].length - segmentIndex - 1;
          return helpers$1.lineString(
            [
              coords[multiFeatureIndex][geometryIndex][segmentIndex],
              coords[multiFeatureIndex][geometryIndex][segmentIndex + 1],
            ],
            properties,
            options
          );
      }
      throw new Error("geojson is invalid");
    }

    /**
     * Finds a particular Point from a GeoJSON using `@turf/meta` indexes.
     *
     * Negative indexes are permitted.
     *
     * @param {FeatureCollection|Feature|Geometry} geojson Any GeoJSON Feature or Geometry
     * @param {Object} [options={}] Optional parameters
     * @param {number} [options.featureIndex=0] Feature Index
     * @param {number} [options.multiFeatureIndex=0] Multi-Feature Index
     * @param {number} [options.geometryIndex=0] Geometry Index
     * @param {number} [options.coordIndex=0] Coord Index
     * @param {Object} [options.properties={}] Translate Properties to output Point
     * @param {BBox} [options.bbox={}] Translate BBox to output Point
     * @param {number|string} [options.id={}] Translate Id to output Point
     * @returns {Feature<Point>} 2-vertex GeoJSON Feature Point
     * @example
     * var multiLine = turf.multiLineString([
     *     [[10, 10], [50, 30], [30, 40]],
     *     [[-10, -10], [-50, -30], [-30, -40]]
     * ]);
     *
     * // First Segment (defaults are 0)
     * turf.findPoint(multiLine);
     * // => Feature<Point<[10, 10]>>
     *
     * // First Segment of the 2nd Multi-Feature
     * turf.findPoint(multiLine, {multiFeatureIndex: 1});
     * // => Feature<Point<[-10, -10]>>
     *
     * // Last Segment of last Multi-Feature
     * turf.findPoint(multiLine, {multiFeatureIndex: -1, coordIndex: -1});
     * // => Feature<Point<[-30, -40]>>
     */
    function findPoint$1(geojson, options) {
      // Optional Parameters
      options = options || {};
      if (!helpers$1.isObject(options)) throw new Error("options is invalid");
      var featureIndex = options.featureIndex || 0;
      var multiFeatureIndex = options.multiFeatureIndex || 0;
      var geometryIndex = options.geometryIndex || 0;
      var coordIndex = options.coordIndex || 0;

      // Find FeatureIndex
      var properties = options.properties;
      var geometry;

      switch (geojson.type) {
        case "FeatureCollection":
          if (featureIndex < 0)
            featureIndex = geojson.features.length + featureIndex;
          properties = properties || geojson.features[featureIndex].properties;
          geometry = geojson.features[featureIndex].geometry;
          break;
        case "Feature":
          properties = properties || geojson.properties;
          geometry = geojson.geometry;
          break;
        case "Point":
        case "MultiPoint":
          return null;
        case "LineString":
        case "Polygon":
        case "MultiLineString":
        case "MultiPolygon":
          geometry = geojson;
          break;
        default:
          throw new Error("geojson is invalid");
      }

      // Find Coord Index
      if (geometry === null) return null;
      var coords = geometry.coordinates;
      switch (geometry.type) {
        case "Point":
          return helpers$1.point(coords, properties, options);
        case "MultiPoint":
          if (multiFeatureIndex < 0)
            multiFeatureIndex = coords.length + multiFeatureIndex;
          return helpers$1.point(coords[multiFeatureIndex], properties, options);
        case "LineString":
          if (coordIndex < 0) coordIndex = coords.length + coordIndex;
          return helpers$1.point(coords[coordIndex], properties, options);
        case "Polygon":
          if (geometryIndex < 0) geometryIndex = coords.length + geometryIndex;
          if (coordIndex < 0)
            coordIndex = coords[geometryIndex].length + coordIndex;
          return helpers$1.point(coords[geometryIndex][coordIndex], properties, options);
        case "MultiLineString":
          if (multiFeatureIndex < 0)
            multiFeatureIndex = coords.length + multiFeatureIndex;
          if (coordIndex < 0)
            coordIndex = coords[multiFeatureIndex].length + coordIndex;
          return helpers$1.point(coords[multiFeatureIndex][coordIndex], properties, options);
        case "MultiPolygon":
          if (multiFeatureIndex < 0)
            multiFeatureIndex = coords.length + multiFeatureIndex;
          if (geometryIndex < 0)
            geometryIndex = coords[multiFeatureIndex].length + geometryIndex;
          if (coordIndex < 0)
            coordIndex =
              coords[multiFeatureIndex][geometryIndex].length - coordIndex;
          return helpers$1.point(
            coords[multiFeatureIndex][geometryIndex][coordIndex],
            properties,
            options
          );
      }
      throw new Error("geojson is invalid");
    }

    js.coordAll = coordAll$1;
    js.coordEach = coordEach$1;
    js.coordReduce = coordReduce$1;
    js.featureEach = featureEach$2;
    js.featureReduce = featureReduce$1;
    js.findPoint = findPoint$1;
    js.findSegment = findSegment$1;
    js.flattenEach = flattenEach$1;
    js.flattenReduce = flattenReduce$1;
    js.geomEach = geomEach$1;
    js.geomReduce = geomReduce$1;
    js.lineEach = lineEach$1;
    js.lineReduce = lineReduce$1;
    js.propEach = propEach$1;
    js.propReduce = propReduce$1;
    js.segmentEach = segmentEach$1;
    js.segmentReduce = segmentReduce$1;

    var cjs$2 = {};

    var cjs$1 = {};

    var cjs = {};

    Object.defineProperty(cjs, "__esModule", {value: true});// index.ts
    var earthRadius = 63710088e-1;
    var factors = {
      centimeters: earthRadius * 100,
      centimetres: earthRadius * 100,
      degrees: 360 / (2 * Math.PI),
      feet: earthRadius * 3.28084,
      inches: earthRadius * 39.37,
      kilometers: earthRadius / 1e3,
      kilometres: earthRadius / 1e3,
      meters: earthRadius,
      metres: earthRadius,
      miles: earthRadius / 1609.344,
      millimeters: earthRadius * 1e3,
      millimetres: earthRadius * 1e3,
      nauticalmiles: earthRadius / 1852,
      radians: 1,
      yards: earthRadius * 1.0936
    };
    var areaFactors = {
      acres: 247105e-9,
      centimeters: 1e4,
      centimetres: 1e4,
      feet: 10.763910417,
      hectares: 1e-4,
      inches: 1550.003100006,
      kilometers: 1e-6,
      kilometres: 1e-6,
      meters: 1,
      metres: 1,
      miles: 386e-9,
      nauticalmiles: 29155334959812285e-23,
      millimeters: 1e6,
      millimetres: 1e6,
      yards: 1.195990046
    };
    function feature(geom, properties, options = {}) {
      const feat = { type: "Feature" };
      if (options.id === 0 || options.id) {
        feat.id = options.id;
      }
      if (options.bbox) {
        feat.bbox = options.bbox;
      }
      feat.properties = properties || {};
      feat.geometry = geom;
      return feat;
    }
    function geometry(type, coordinates, _options = {}) {
      switch (type) {
        case "Point":
          return point(coordinates).geometry;
        case "LineString":
          return lineString(coordinates).geometry;
        case "Polygon":
          return polygon(coordinates).geometry;
        case "MultiPoint":
          return multiPoint(coordinates).geometry;
        case "MultiLineString":
          return multiLineString(coordinates).geometry;
        case "MultiPolygon":
          return multiPolygon(coordinates).geometry;
        default:
          throw new Error(type + " is invalid");
      }
    }
    function point(coordinates, properties, options = {}) {
      if (!coordinates) {
        throw new Error("coordinates is required");
      }
      if (!Array.isArray(coordinates)) {
        throw new Error("coordinates must be an Array");
      }
      if (coordinates.length < 2) {
        throw new Error("coordinates must be at least 2 numbers long");
      }
      if (!isNumber(coordinates[0]) || !isNumber(coordinates[1])) {
        throw new Error("coordinates must contain numbers");
      }
      const geom = {
        type: "Point",
        coordinates
      };
      return feature(geom, properties, options);
    }
    function points(coordinates, properties, options = {}) {
      return featureCollection$1(
        coordinates.map((coords) => {
          return point(coords, properties);
        }),
        options
      );
    }
    function polygon(coordinates, properties, options = {}) {
      for (const ring of coordinates) {
        if (ring.length < 4) {
          throw new Error(
            "Each LinearRing of a Polygon must have 4 or more Positions."
          );
        }
        if (ring[ring.length - 1].length !== ring[0].length) {
          throw new Error("First and last Position are not equivalent.");
        }
        for (let j = 0; j < ring[ring.length - 1].length; j++) {
          if (ring[ring.length - 1][j] !== ring[0][j]) {
            throw new Error("First and last Position are not equivalent.");
          }
        }
      }
      const geom = {
        type: "Polygon",
        coordinates
      };
      return feature(geom, properties, options);
    }
    function polygons(coordinates, properties, options = {}) {
      return featureCollection$1(
        coordinates.map((coords) => {
          return polygon(coords, properties);
        }),
        options
      );
    }
    function lineString(coordinates, properties, options = {}) {
      if (coordinates.length < 2) {
        throw new Error("coordinates must be an array of two or more positions");
      }
      const geom = {
        type: "LineString",
        coordinates
      };
      return feature(geom, properties, options);
    }
    function lineStrings(coordinates, properties, options = {}) {
      return featureCollection$1(
        coordinates.map((coords) => {
          return lineString(coords, properties);
        }),
        options
      );
    }
    function featureCollection$1(features, options = {}) {
      const fc = { type: "FeatureCollection" };
      if (options.id) {
        fc.id = options.id;
      }
      if (options.bbox) {
        fc.bbox = options.bbox;
      }
      fc.features = features;
      return fc;
    }
    function multiLineString(coordinates, properties, options = {}) {
      const geom = {
        type: "MultiLineString",
        coordinates
      };
      return feature(geom, properties, options);
    }
    function multiPoint(coordinates, properties, options = {}) {
      const geom = {
        type: "MultiPoint",
        coordinates
      };
      return feature(geom, properties, options);
    }
    function multiPolygon(coordinates, properties, options = {}) {
      const geom = {
        type: "MultiPolygon",
        coordinates
      };
      return feature(geom, properties, options);
    }
    function geometryCollection(geometries, properties, options = {}) {
      const geom = {
        type: "GeometryCollection",
        geometries
      };
      return feature(geom, properties, options);
    }
    function round(num, precision = 0) {
      if (precision && !(precision >= 0)) {
        throw new Error("precision must be a positive number");
      }
      const multiplier = Math.pow(10, precision || 0);
      return Math.round(num * multiplier) / multiplier;
    }
    function radiansToLength(radians, units = "kilometers") {
      const factor = factors[units];
      if (!factor) {
        throw new Error(units + " units is invalid");
      }
      return radians * factor;
    }
    function lengthToRadians(distance, units = "kilometers") {
      const factor = factors[units];
      if (!factor) {
        throw new Error(units + " units is invalid");
      }
      return distance / factor;
    }
    function lengthToDegrees(distance, units) {
      return radiansToDegrees(lengthToRadians(distance, units));
    }
    function bearingToAzimuth(bearing) {
      let angle = bearing % 360;
      if (angle < 0) {
        angle += 360;
      }
      return angle;
    }
    function azimuthToBearing(angle) {
      angle = angle % 360;
      if (angle > 0)
        return angle > 180 ? angle - 360 : angle;
      return angle < -180 ? angle + 360 : angle;
    }
    function radiansToDegrees(radians) {
      const degrees = radians % (2 * Math.PI);
      return degrees * 180 / Math.PI;
    }
    function degreesToRadians(degrees) {
      const radians = degrees % 360;
      return radians * Math.PI / 180;
    }
    function convertLength(length, originalUnit = "kilometers", finalUnit = "kilometers") {
      if (!(length >= 0)) {
        throw new Error("length must be a positive number");
      }
      return radiansToLength(lengthToRadians(length, originalUnit), finalUnit);
    }
    function convertArea(area, originalUnit = "meters", finalUnit = "kilometers") {
      if (!(area >= 0)) {
        throw new Error("area must be a positive number");
      }
      const startFactor = areaFactors[originalUnit];
      if (!startFactor) {
        throw new Error("invalid original units");
      }
      const finalFactor = areaFactors[finalUnit];
      if (!finalFactor) {
        throw new Error("invalid final units");
      }
      return area / startFactor * finalFactor;
    }
    function isNumber(num) {
      return !isNaN(num) && num !== null && !Array.isArray(num);
    }
    function isObject(input) {
      return input !== null && typeof input === "object" && !Array.isArray(input);
    }
    function validateBBox(bbox) {
      if (!bbox) {
        throw new Error("bbox is required");
      }
      if (!Array.isArray(bbox)) {
        throw new Error("bbox must be an Array");
      }
      if (bbox.length !== 4 && bbox.length !== 6) {
        throw new Error("bbox must be an Array of 4 or 6 numbers");
      }
      bbox.forEach((num) => {
        if (!isNumber(num)) {
          throw new Error("bbox must only contain numbers");
        }
      });
    }
    function validateId(id) {
      if (!id) {
        throw new Error("id is required");
      }
      if (["string", "number"].indexOf(typeof id) === -1) {
        throw new Error("id must be a number or a string");
      }
    }































    cjs.areaFactors = areaFactors; cjs.azimuthToBearing = azimuthToBearing; cjs.bearingToAzimuth = bearingToAzimuth; cjs.convertArea = convertArea; cjs.convertLength = convertLength; cjs.degreesToRadians = degreesToRadians; cjs.earthRadius = earthRadius; cjs.factors = factors; cjs.feature = feature; cjs.featureCollection = featureCollection$1; cjs.geometry = geometry; cjs.geometryCollection = geometryCollection; cjs.isNumber = isNumber; cjs.isObject = isObject; cjs.lengthToDegrees = lengthToDegrees; cjs.lengthToRadians = lengthToRadians; cjs.lineString = lineString; cjs.lineStrings = lineStrings; cjs.multiLineString = multiLineString; cjs.multiPoint = multiPoint; cjs.multiPolygon = multiPolygon; cjs.point = point; cjs.points = points; cjs.polygon = polygon; cjs.polygons = polygons; cjs.radiansToDegrees = radiansToDegrees; cjs.radiansToLength = radiansToLength; cjs.round = round; cjs.validateBBox = validateBBox; cjs.validateId = validateId;

    Object.defineProperty(cjs$1, "__esModule", {value: true});// index.js
    var _helpers = cjs;
    function coordEach(geojson, callback, excludeWrapCoord) {
      if (geojson === null)
        return;
      var j, k, l, geometry, stopG, coords, geometryMaybeCollection, wrapShrink = 0, coordIndex = 0, isGeometryCollection, type = geojson.type, isFeatureCollection = type === "FeatureCollection", isFeature = type === "Feature", stop = isFeatureCollection ? geojson.features.length : 1;
      for (var featureIndex = 0; featureIndex < stop; featureIndex++) {
        geometryMaybeCollection = isFeatureCollection ? geojson.features[featureIndex].geometry : isFeature ? geojson.geometry : geojson;
        isGeometryCollection = geometryMaybeCollection ? geometryMaybeCollection.type === "GeometryCollection" : false;
        stopG = isGeometryCollection ? geometryMaybeCollection.geometries.length : 1;
        for (var geomIndex = 0; geomIndex < stopG; geomIndex++) {
          var multiFeatureIndex = 0;
          var geometryIndex = 0;
          geometry = isGeometryCollection ? geometryMaybeCollection.geometries[geomIndex] : geometryMaybeCollection;
          if (geometry === null)
            continue;
          coords = geometry.coordinates;
          var geomType = geometry.type;
          wrapShrink = excludeWrapCoord && (geomType === "Polygon" || geomType === "MultiPolygon") ? 1 : 0;
          switch (geomType) {
            case null:
              break;
            case "Point":
              if (callback(
                coords,
                coordIndex,
                featureIndex,
                multiFeatureIndex,
                geometryIndex
              ) === false)
                return false;
              coordIndex++;
              multiFeatureIndex++;
              break;
            case "LineString":
            case "MultiPoint":
              for (j = 0; j < coords.length; j++) {
                if (callback(
                  coords[j],
                  coordIndex,
                  featureIndex,
                  multiFeatureIndex,
                  geometryIndex
                ) === false)
                  return false;
                coordIndex++;
                if (geomType === "MultiPoint")
                  multiFeatureIndex++;
              }
              if (geomType === "LineString")
                multiFeatureIndex++;
              break;
            case "Polygon":
            case "MultiLineString":
              for (j = 0; j < coords.length; j++) {
                for (k = 0; k < coords[j].length - wrapShrink; k++) {
                  if (callback(
                    coords[j][k],
                    coordIndex,
                    featureIndex,
                    multiFeatureIndex,
                    geometryIndex
                  ) === false)
                    return false;
                  coordIndex++;
                }
                if (geomType === "MultiLineString")
                  multiFeatureIndex++;
                if (geomType === "Polygon")
                  geometryIndex++;
              }
              if (geomType === "Polygon")
                multiFeatureIndex++;
              break;
            case "MultiPolygon":
              for (j = 0; j < coords.length; j++) {
                geometryIndex = 0;
                for (k = 0; k < coords[j].length; k++) {
                  for (l = 0; l < coords[j][k].length - wrapShrink; l++) {
                    if (callback(
                      coords[j][k][l],
                      coordIndex,
                      featureIndex,
                      multiFeatureIndex,
                      geometryIndex
                    ) === false)
                      return false;
                    coordIndex++;
                  }
                  geometryIndex++;
                }
                multiFeatureIndex++;
              }
              break;
            case "GeometryCollection":
              for (j = 0; j < geometry.geometries.length; j++)
                if (coordEach(geometry.geometries[j], callback, excludeWrapCoord) === false)
                  return false;
              break;
            default:
              throw new Error("Unknown Geometry Type");
          }
        }
      }
    }
    function coordReduce(geojson, callback, initialValue, excludeWrapCoord) {
      var previousValue = initialValue;
      coordEach(
        geojson,
        function(currentCoord, coordIndex, featureIndex, multiFeatureIndex, geometryIndex) {
          if (coordIndex === 0 && initialValue === void 0)
            previousValue = currentCoord;
          else
            previousValue = callback(
              previousValue,
              currentCoord,
              coordIndex,
              featureIndex,
              multiFeatureIndex,
              geometryIndex
            );
        },
        excludeWrapCoord
      );
      return previousValue;
    }
    function propEach(geojson, callback) {
      var i;
      switch (geojson.type) {
        case "FeatureCollection":
          for (i = 0; i < geojson.features.length; i++) {
            if (callback(geojson.features[i].properties, i) === false)
              break;
          }
          break;
        case "Feature":
          callback(geojson.properties, 0);
          break;
      }
    }
    function propReduce(geojson, callback, initialValue) {
      var previousValue = initialValue;
      propEach(geojson, function(currentProperties, featureIndex) {
        if (featureIndex === 0 && initialValue === void 0)
          previousValue = currentProperties;
        else
          previousValue = callback(previousValue, currentProperties, featureIndex);
      });
      return previousValue;
    }
    function featureEach$1(geojson, callback) {
      if (geojson.type === "Feature") {
        callback(geojson, 0);
      } else if (geojson.type === "FeatureCollection") {
        for (var i = 0; i < geojson.features.length; i++) {
          if (callback(geojson.features[i], i) === false)
            break;
        }
      }
    }
    function featureReduce(geojson, callback, initialValue) {
      var previousValue = initialValue;
      featureEach$1(geojson, function(currentFeature, featureIndex) {
        if (featureIndex === 0 && initialValue === void 0)
          previousValue = currentFeature;
        else
          previousValue = callback(previousValue, currentFeature, featureIndex);
      });
      return previousValue;
    }
    function coordAll(geojson) {
      var coords = [];
      coordEach(geojson, function(coord) {
        coords.push(coord);
      });
      return coords;
    }
    function geomEach(geojson, callback) {
      var i, j, g, geometry, stopG, geometryMaybeCollection, isGeometryCollection, featureProperties, featureBBox, featureId, featureIndex = 0, isFeatureCollection = geojson.type === "FeatureCollection", isFeature = geojson.type === "Feature", stop = isFeatureCollection ? geojson.features.length : 1;
      for (i = 0; i < stop; i++) {
        geometryMaybeCollection = isFeatureCollection ? geojson.features[i].geometry : isFeature ? geojson.geometry : geojson;
        featureProperties = isFeatureCollection ? geojson.features[i].properties : isFeature ? geojson.properties : {};
        featureBBox = isFeatureCollection ? geojson.features[i].bbox : isFeature ? geojson.bbox : void 0;
        featureId = isFeatureCollection ? geojson.features[i].id : isFeature ? geojson.id : void 0;
        isGeometryCollection = geometryMaybeCollection ? geometryMaybeCollection.type === "GeometryCollection" : false;
        stopG = isGeometryCollection ? geometryMaybeCollection.geometries.length : 1;
        for (g = 0; g < stopG; g++) {
          geometry = isGeometryCollection ? geometryMaybeCollection.geometries[g] : geometryMaybeCollection;
          if (geometry === null) {
            if (callback(
              null,
              featureIndex,
              featureProperties,
              featureBBox,
              featureId
            ) === false)
              return false;
            continue;
          }
          switch (geometry.type) {
            case "Point":
            case "LineString":
            case "MultiPoint":
            case "Polygon":
            case "MultiLineString":
            case "MultiPolygon": {
              if (callback(
                geometry,
                featureIndex,
                featureProperties,
                featureBBox,
                featureId
              ) === false)
                return false;
              break;
            }
            case "GeometryCollection": {
              for (j = 0; j < geometry.geometries.length; j++) {
                if (callback(
                  geometry.geometries[j],
                  featureIndex,
                  featureProperties,
                  featureBBox,
                  featureId
                ) === false)
                  return false;
              }
              break;
            }
            default:
              throw new Error("Unknown Geometry Type");
          }
        }
        featureIndex++;
      }
    }
    function geomReduce(geojson, callback, initialValue) {
      var previousValue = initialValue;
      geomEach(
        geojson,
        function(currentGeometry, featureIndex, featureProperties, featureBBox, featureId) {
          if (featureIndex === 0 && initialValue === void 0)
            previousValue = currentGeometry;
          else
            previousValue = callback(
              previousValue,
              currentGeometry,
              featureIndex,
              featureProperties,
              featureBBox,
              featureId
            );
        }
      );
      return previousValue;
    }
    function flattenEach(geojson, callback) {
      geomEach(geojson, function(geometry, featureIndex, properties, bbox, id) {
        var type = geometry === null ? null : geometry.type;
        switch (type) {
          case null:
          case "Point":
          case "LineString":
          case "Polygon":
            if (callback(
              _helpers.feature.call(void 0, geometry, properties, { bbox, id }),
              featureIndex,
              0
            ) === false)
              return false;
            return;
        }
        var geomType;
        switch (type) {
          case "MultiPoint":
            geomType = "Point";
            break;
          case "MultiLineString":
            geomType = "LineString";
            break;
          case "MultiPolygon":
            geomType = "Polygon";
            break;
        }
        for (var multiFeatureIndex = 0; multiFeatureIndex < geometry.coordinates.length; multiFeatureIndex++) {
          var coordinate = geometry.coordinates[multiFeatureIndex];
          var geom = {
            type: geomType,
            coordinates: coordinate
          };
          if (callback(_helpers.feature.call(void 0, geom, properties), featureIndex, multiFeatureIndex) === false)
            return false;
        }
      });
    }
    function flattenReduce(geojson, callback, initialValue) {
      var previousValue = initialValue;
      flattenEach(
        geojson,
        function(currentFeature, featureIndex, multiFeatureIndex) {
          if (featureIndex === 0 && multiFeatureIndex === 0 && initialValue === void 0)
            previousValue = currentFeature;
          else
            previousValue = callback(
              previousValue,
              currentFeature,
              featureIndex,
              multiFeatureIndex
            );
        }
      );
      return previousValue;
    }
    function segmentEach(geojson, callback) {
      flattenEach(geojson, function(feature2, featureIndex, multiFeatureIndex) {
        var segmentIndex = 0;
        if (!feature2.geometry)
          return;
        var type = feature2.geometry.type;
        if (type === "Point" || type === "MultiPoint")
          return;
        var previousCoords;
        var previousFeatureIndex = 0;
        var previousMultiIndex = 0;
        var prevGeomIndex = 0;
        if (coordEach(
          feature2,
          function(currentCoord, coordIndex, featureIndexCoord, multiPartIndexCoord, geometryIndex) {
            if (previousCoords === void 0 || featureIndex > previousFeatureIndex || multiPartIndexCoord > previousMultiIndex || geometryIndex > prevGeomIndex) {
              previousCoords = currentCoord;
              previousFeatureIndex = featureIndex;
              previousMultiIndex = multiPartIndexCoord;
              prevGeomIndex = geometryIndex;
              segmentIndex = 0;
              return;
            }
            var currentSegment = _helpers.lineString.call(void 0, 
              [previousCoords, currentCoord],
              feature2.properties
            );
            if (callback(
              currentSegment,
              featureIndex,
              multiFeatureIndex,
              geometryIndex,
              segmentIndex
            ) === false)
              return false;
            segmentIndex++;
            previousCoords = currentCoord;
          }
        ) === false)
          return false;
      });
    }
    function segmentReduce(geojson, callback, initialValue) {
      var previousValue = initialValue;
      var started = false;
      segmentEach(
        geojson,
        function(currentSegment, featureIndex, multiFeatureIndex, geometryIndex, segmentIndex) {
          if (started === false && initialValue === void 0)
            previousValue = currentSegment;
          else
            previousValue = callback(
              previousValue,
              currentSegment,
              featureIndex,
              multiFeatureIndex,
              geometryIndex,
              segmentIndex
            );
          started = true;
        }
      );
      return previousValue;
    }
    function lineEach(geojson, callback) {
      if (!geojson)
        throw new Error("geojson is required");
      flattenEach(geojson, function(feature2, featureIndex, multiFeatureIndex) {
        if (feature2.geometry === null)
          return;
        var type = feature2.geometry.type;
        var coords = feature2.geometry.coordinates;
        switch (type) {
          case "LineString":
            if (callback(feature2, featureIndex, multiFeatureIndex, 0, 0) === false)
              return false;
            break;
          case "Polygon":
            for (var geometryIndex = 0; geometryIndex < coords.length; geometryIndex++) {
              if (callback(
                _helpers.lineString.call(void 0, coords[geometryIndex], feature2.properties),
                featureIndex,
                multiFeatureIndex,
                geometryIndex
              ) === false)
                return false;
            }
            break;
        }
      });
    }
    function lineReduce(geojson, callback, initialValue) {
      var previousValue = initialValue;
      lineEach(
        geojson,
        function(currentLine, featureIndex, multiFeatureIndex, geometryIndex) {
          if (featureIndex === 0 && initialValue === void 0)
            previousValue = currentLine;
          else
            previousValue = callback(
              previousValue,
              currentLine,
              featureIndex,
              multiFeatureIndex,
              geometryIndex
            );
        }
      );
      return previousValue;
    }
    function findSegment(geojson, options) {
      options = options || {};
      if (!_helpers.isObject.call(void 0, options))
        throw new Error("options is invalid");
      var featureIndex = options.featureIndex || 0;
      var multiFeatureIndex = options.multiFeatureIndex || 0;
      var geometryIndex = options.geometryIndex || 0;
      var segmentIndex = options.segmentIndex || 0;
      var properties = options.properties;
      var geometry;
      switch (geojson.type) {
        case "FeatureCollection":
          if (featureIndex < 0)
            featureIndex = geojson.features.length + featureIndex;
          properties = properties || geojson.features[featureIndex].properties;
          geometry = geojson.features[featureIndex].geometry;
          break;
        case "Feature":
          properties = properties || geojson.properties;
          geometry = geojson.geometry;
          break;
        case "Point":
        case "MultiPoint":
          return null;
        case "LineString":
        case "Polygon":
        case "MultiLineString":
        case "MultiPolygon":
          geometry = geojson;
          break;
        default:
          throw new Error("geojson is invalid");
      }
      if (geometry === null)
        return null;
      var coords = geometry.coordinates;
      switch (geometry.type) {
        case "Point":
        case "MultiPoint":
          return null;
        case "LineString":
          if (segmentIndex < 0)
            segmentIndex = coords.length + segmentIndex - 1;
          return _helpers.lineString.call(void 0, 
            [coords[segmentIndex], coords[segmentIndex + 1]],
            properties,
            options
          );
        case "Polygon":
          if (geometryIndex < 0)
            geometryIndex = coords.length + geometryIndex;
          if (segmentIndex < 0)
            segmentIndex = coords[geometryIndex].length + segmentIndex - 1;
          return _helpers.lineString.call(void 0, 
            [
              coords[geometryIndex][segmentIndex],
              coords[geometryIndex][segmentIndex + 1]
            ],
            properties,
            options
          );
        case "MultiLineString":
          if (multiFeatureIndex < 0)
            multiFeatureIndex = coords.length + multiFeatureIndex;
          if (segmentIndex < 0)
            segmentIndex = coords[multiFeatureIndex].length + segmentIndex - 1;
          return _helpers.lineString.call(void 0, 
            [
              coords[multiFeatureIndex][segmentIndex],
              coords[multiFeatureIndex][segmentIndex + 1]
            ],
            properties,
            options
          );
        case "MultiPolygon":
          if (multiFeatureIndex < 0)
            multiFeatureIndex = coords.length + multiFeatureIndex;
          if (geometryIndex < 0)
            geometryIndex = coords[multiFeatureIndex].length + geometryIndex;
          if (segmentIndex < 0)
            segmentIndex = coords[multiFeatureIndex][geometryIndex].length - segmentIndex - 1;
          return _helpers.lineString.call(void 0, 
            [
              coords[multiFeatureIndex][geometryIndex][segmentIndex],
              coords[multiFeatureIndex][geometryIndex][segmentIndex + 1]
            ],
            properties,
            options
          );
      }
      throw new Error("geojson is invalid");
    }
    function findPoint(geojson, options) {
      options = options || {};
      if (!_helpers.isObject.call(void 0, options))
        throw new Error("options is invalid");
      var featureIndex = options.featureIndex || 0;
      var multiFeatureIndex = options.multiFeatureIndex || 0;
      var geometryIndex = options.geometryIndex || 0;
      var coordIndex = options.coordIndex || 0;
      var properties = options.properties;
      var geometry;
      switch (geojson.type) {
        case "FeatureCollection":
          if (featureIndex < 0)
            featureIndex = geojson.features.length + featureIndex;
          properties = properties || geojson.features[featureIndex].properties;
          geometry = geojson.features[featureIndex].geometry;
          break;
        case "Feature":
          properties = properties || geojson.properties;
          geometry = geojson.geometry;
          break;
        case "Point":
        case "MultiPoint":
          return null;
        case "LineString":
        case "Polygon":
        case "MultiLineString":
        case "MultiPolygon":
          geometry = geojson;
          break;
        default:
          throw new Error("geojson is invalid");
      }
      if (geometry === null)
        return null;
      var coords = geometry.coordinates;
      switch (geometry.type) {
        case "Point":
          return _helpers.point.call(void 0, coords, properties, options);
        case "MultiPoint":
          if (multiFeatureIndex < 0)
            multiFeatureIndex = coords.length + multiFeatureIndex;
          return _helpers.point.call(void 0, coords[multiFeatureIndex], properties, options);
        case "LineString":
          if (coordIndex < 0)
            coordIndex = coords.length + coordIndex;
          return _helpers.point.call(void 0, coords[coordIndex], properties, options);
        case "Polygon":
          if (geometryIndex < 0)
            geometryIndex = coords.length + geometryIndex;
          if (coordIndex < 0)
            coordIndex = coords[geometryIndex].length + coordIndex;
          return _helpers.point.call(void 0, coords[geometryIndex][coordIndex], properties, options);
        case "MultiLineString":
          if (multiFeatureIndex < 0)
            multiFeatureIndex = coords.length + multiFeatureIndex;
          if (coordIndex < 0)
            coordIndex = coords[multiFeatureIndex].length + coordIndex;
          return _helpers.point.call(void 0, coords[multiFeatureIndex][coordIndex], properties, options);
        case "MultiPolygon":
          if (multiFeatureIndex < 0)
            multiFeatureIndex = coords.length + multiFeatureIndex;
          if (geometryIndex < 0)
            geometryIndex = coords[multiFeatureIndex].length + geometryIndex;
          if (coordIndex < 0)
            coordIndex = coords[multiFeatureIndex][geometryIndex].length - coordIndex;
          return _helpers.point.call(void 0, 
            coords[multiFeatureIndex][geometryIndex][coordIndex],
            properties,
            options
          );
      }
      throw new Error("geojson is invalid");
    }


















    cjs$1.coordAll = coordAll; cjs$1.coordEach = coordEach; cjs$1.coordReduce = coordReduce; cjs$1.featureEach = featureEach$1; cjs$1.featureReduce = featureReduce; cjs$1.findPoint = findPoint; cjs$1.findSegment = findSegment; cjs$1.flattenEach = flattenEach; cjs$1.flattenReduce = flattenReduce; cjs$1.geomEach = geomEach; cjs$1.geomReduce = geomReduce; cjs$1.lineEach = lineEach; cjs$1.lineReduce = lineReduce; cjs$1.propEach = propEach; cjs$1.propReduce = propReduce; cjs$1.segmentEach = segmentEach; cjs$1.segmentReduce = segmentReduce;

    Object.defineProperty(cjs$2, "__esModule", {value: true});// index.ts
    var _meta = cjs$1;
    function bbox(geojson, options = {}) {
      if (geojson.bbox != null && true !== options.recompute) {
        return geojson.bbox;
      }
      const result = [Infinity, Infinity, -Infinity, -Infinity];
      _meta.coordEach.call(void 0, geojson, (coord) => {
        if (result[0] > coord[0]) {
          result[0] = coord[0];
        }
        if (result[1] > coord[1]) {
          result[1] = coord[1];
        }
        if (result[2] < coord[0]) {
          result[2] = coord[0];
        }
        if (result[3] < coord[1]) {
          result[3] = coord[1];
        }
      });
      return result;
    }
    var turf_bbox_default = bbox;



    cjs$2.bbox = bbox; cjs$2.default = turf_bbox_default;

    var rbush = require$$0;
    var helpers = js$1;
    var meta = js;
    var turfBBox = cjs$2.default;
    var featureEach = meta.featureEach;
    meta.coordEach;
    helpers.polygon;
    var featureCollection = helpers.featureCollection;

    /**
     * GeoJSON implementation of [RBush](https://github.com/mourner/rbush#rbush) spatial index.
     *
     * @name rbush
     * @param {number} [maxEntries=9] defines the maximum number of entries in a tree node. 9 (used by default) is a
     * reasonable choice for most applications. Higher value means faster insertion and slower search, and vice versa.
     * @returns {RBush} GeoJSON RBush
     * @example
     * var geojsonRbush = require('geojson-rbush').default;
     * var tree = geojsonRbush();
     */
    function geojsonRbush(maxEntries) {
        var tree = new rbush(maxEntries);
        /**
         * [insert](https://github.com/mourner/rbush#data-format)
         *
         * @param {Feature} feature insert single GeoJSON Feature
         * @returns {RBush} GeoJSON RBush
         * @example
         * var poly = turf.polygon([[[-78, 41], [-67, 41], [-67, 48], [-78, 48], [-78, 41]]]);
         * tree.insert(poly)
         */
        tree.insert = function (feature) {
            if (feature.type !== 'Feature') throw new Error('invalid feature');
            feature.bbox = feature.bbox ? feature.bbox : turfBBox(feature);
            return rbush.prototype.insert.call(this, feature);
        };

        /**
         * [load](https://github.com/mourner/rbush#bulk-inserting-data)
         *
         * @param {FeatureCollection|Array<Feature>} features load entire GeoJSON FeatureCollection
         * @returns {RBush} GeoJSON RBush
         * @example
         * var polys = turf.polygons([
         *     [[[-78, 41], [-67, 41], [-67, 48], [-78, 48], [-78, 41]]],
         *     [[[-93, 32], [-83, 32], [-83, 39], [-93, 39], [-93, 32]]]
         * ]);
         * tree.load(polys);
         */
        tree.load = function (features) {
            var load = [];
            // Load an Array of Features
            if (Array.isArray(features)) {
                features.forEach(function (feature) {
                    if (feature.type !== 'Feature') throw new Error('invalid features');
                    feature.bbox = feature.bbox ? feature.bbox : turfBBox(feature);
                    load.push(feature);
                });
            } else {
                // Load a FeatureCollection
                featureEach(features, function (feature) {
                    if (feature.type !== 'Feature') throw new Error('invalid features');
                    feature.bbox = feature.bbox ? feature.bbox : turfBBox(feature);
                    load.push(feature);
                });
            }
            return rbush.prototype.load.call(this, load);
        };

        /**
         * [remove](https://github.com/mourner/rbush#removing-data)
         *
         * @param {Feature} feature remove single GeoJSON Feature
         * @param {Function} equals Pass a custom equals function to compare by value for removal.
         * @returns {RBush} GeoJSON RBush
         * @example
         * var poly = turf.polygon([[[-78, 41], [-67, 41], [-67, 48], [-78, 48], [-78, 41]]]);
         *
         * tree.remove(poly);
         */
        tree.remove = function (feature, equals) {
            if (feature.type !== 'Feature') throw new Error('invalid feature');
            feature.bbox = feature.bbox ? feature.bbox : turfBBox(feature);
            return rbush.prototype.remove.call(this, feature, equals);
        };

        /**
         * [clear](https://github.com/mourner/rbush#removing-data)
         *
         * @returns {RBush} GeoJSON Rbush
         * @example
         * tree.clear()
         */
        tree.clear = function () {
            return rbush.prototype.clear.call(this);
        };

        /**
         * [search](https://github.com/mourner/rbush#search)
         *
         * @param {BBox|FeatureCollection|Feature} geojson search with GeoJSON
         * @returns {FeatureCollection} all features that intersects with the given GeoJSON.
         * @example
         * var poly = turf.polygon([[[-78, 41], [-67, 41], [-67, 48], [-78, 48], [-78, 41]]]);
         *
         * tree.search(poly);
         */
        tree.search = function (geojson) {
            var features = rbush.prototype.search.call(this, this.toBBox(geojson));
            return featureCollection(features);
        };

        /**
         * [collides](https://github.com/mourner/rbush#collisions)
         *
         * @param {BBox|FeatureCollection|Feature} geojson collides with GeoJSON
         * @returns {boolean} true if there are any items intersecting the given GeoJSON, otherwise false.
         * @example
         * var poly = turf.polygon([[[-78, 41], [-67, 41], [-67, 48], [-78, 48], [-78, 41]]]);
         *
         * tree.collides(poly);
         */
        tree.collides = function (geojson) {
            return rbush.prototype.collides.call(this, this.toBBox(geojson));
        };

        /**
         * [all](https://github.com/mourner/rbush#search)
         *
         * @returns {FeatureCollection} all the features in RBush
         * @example
         * tree.all()
         */
        tree.all = function () {
            var features = rbush.prototype.all.call(this);
            return featureCollection(features);
        };

        /**
         * [toJSON](https://github.com/mourner/rbush#export-and-import)
         *
         * @returns {any} export data as JSON object
         * @example
         * var exported = tree.toJSON()
         */
        tree.toJSON = function () {
            return rbush.prototype.toJSON.call(this);
        };

        /**
         * [fromJSON](https://github.com/mourner/rbush#export-and-import)
         *
         * @param {any} json import previously exported data
         * @returns {RBush} GeoJSON RBush
         * @example
         * var exported = {
         *   "children": [
         *     {
         *       "type": "Feature",
         *       "geometry": {
         *         "type": "Point",
         *         "coordinates": [110, 50]
         *       },
         *       "properties": {},
         *       "bbox": [110, 50, 110, 50]
         *     }
         *   ],
         *   "height": 1,
         *   "leaf": true,
         *   "minX": 110,
         *   "minY": 50,
         *   "maxX": 110,
         *   "maxY": 50
         * }
         * tree.fromJSON(exported)
         */
        tree.fromJSON = function (json) {
            return rbush.prototype.fromJSON.call(this, json);
        };

        /**
         * Converts GeoJSON to {minX, minY, maxX, maxY} schema
         *
         * @private
         * @param {BBox|FeatureCollection|Feature} geojson feature(s) to retrieve BBox from
         * @returns {Object} converted to {minX, minY, maxX, maxY}
         */
        tree.toBBox = function (geojson) {
            var bbox;
            if (geojson.bbox) bbox = geojson.bbox;
            else if (Array.isArray(geojson) && geojson.length === 4) bbox = geojson;
            else if (Array.isArray(geojson) && geojson.length === 6) bbox = [geojson[0], geojson[1], geojson[3], geojson[4]];
            else if (geojson.type === 'Feature') bbox = turfBBox(geojson);
            else if (geojson.type === 'FeatureCollection') bbox = turfBBox(geojson);
            else throw new Error('invalid geojson')

            return {
                minX: bbox[0],
                minY: bbox[1],
                maxX: bbox[2],
                maxY: bbox[3]
            };
        };
        return tree;
    }

    geojsonRbush$1.exports = geojsonRbush;
    geojsonRbush$1.exports.default = geojsonRbush;

    var geojsonRbushExports = geojsonRbush$1.exports;
    var rbush$1 = /*@__PURE__*/getDefaultExportFromCjs(geojsonRbushExports);

    /**
     * Takes any LineString or Polygon GeoJSON and returns the intersecting point(s).
     *
     * @name lineIntersect
     * @param {GeoJSON} line1 any LineString or Polygon
     * @param {GeoJSON} line2 any LineString or Polygon
     * @returns {FeatureCollection<Point>} point(s) that intersect both
     * @example
     * var line1 = turf.lineString([[126, -11], [129, -21]]);
     * var line2 = turf.lineString([[123, -18], [131, -14]]);
     * var intersects = turf.lineIntersect(line1, line2);
     *
     * //addToMap
     * var addToMap = [line1, line2, intersects]
     */
    function lineIntersect(line1, line2) {
        var unique = {};
        var results = [];
        // First, normalize geometries to features
        // Then, handle simple 2-vertex segments
        if (line1.type === "LineString") {
            line1 = feature$2(line1);
        }
        if (line2.type === "LineString") {
            line2 = feature$2(line2);
        }
        if (line1.type === "Feature" &&
            line2.type === "Feature" &&
            line1.geometry !== null &&
            line2.geometry !== null &&
            line1.geometry.type === "LineString" &&
            line2.geometry.type === "LineString" &&
            line1.geometry.coordinates.length === 2 &&
            line2.geometry.coordinates.length === 2) {
            var intersect = intersects(line1, line2);
            if (intersect) {
                results.push(intersect);
            }
            return featureCollection$2(results);
        }
        // Handles complex GeoJSON Geometries
        var tree = rbush$1();
        tree.load(lineSegment(line2));
        featureEach$3(lineSegment(line1), function (segment) {
            featureEach$3(tree.search(segment), function (match) {
                var intersect = intersects(segment, match);
                if (intersect) {
                    // prevent duplicate points https://github.com/Turfjs/turf/issues/688
                    var key = getCoords(intersect).join(",");
                    if (!unique[key]) {
                        unique[key] = true;
                        results.push(intersect);
                    }
                }
            });
        });
        return featureCollection$2(results);
    }
    /**
     * Find a point that intersects LineStrings with two coordinates each
     *
     * @private
     * @param {Feature<LineString>} line1 GeoJSON LineString (Must only contain 2 coordinates)
     * @param {Feature<LineString>} line2 GeoJSON LineString (Must only contain 2 coordinates)
     * @returns {Feature<Point>} intersecting GeoJSON Point
     */
    function intersects(line1, line2) {
        var coords1 = getCoords(line1);
        var coords2 = getCoords(line2);
        if (coords1.length !== 2) {
            throw new Error("<intersects> line1 must only contain 2 coordinates");
        }
        if (coords2.length !== 2) {
            throw new Error("<intersects> line2 must only contain 2 coordinates");
        }
        var x1 = coords1[0][0];
        var y1 = coords1[0][1];
        var x2 = coords1[1][0];
        var y2 = coords1[1][1];
        var x3 = coords2[0][0];
        var y3 = coords2[0][1];
        var x4 = coords2[1][0];
        var y4 = coords2[1][1];
        var denom = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1);
        var numeA = (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3);
        var numeB = (x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3);
        if (denom === 0) {
            if (numeA === 0 && numeB === 0) {
                return null;
            }
            return null;
        }
        var uA = numeA / denom;
        var uB = numeB / denom;
        if (uA >= 0 && uA <= 1 && uB >= 0 && uB <= 1) {
            var x = x1 + uA * (x2 - x1);
            var y = y1 + uA * (y2 - y1);
            return point$2([x, y]);
        }
        return null;
    }

    /**
     * Takes a {@link Point} and a {@link LineString} and calculates the closest Point on the (Multi)LineString.
     *
     * @name nearestPointOnLine
     * @param {Geometry|Feature<LineString|MultiLineString>} lines lines to snap to
     * @param {Geometry|Feature<Point>|number[]} pt point to snap from
     * @param {Object} [options={}] Optional parameters
     * @param {string} [options.units='kilometers'] can be degrees, radians, miles, or kilometers
     * @returns {Feature<Point>} closest point on the `line` to `point`. The properties object will contain three values: `index`: closest point was found on nth line part, `dist`: distance between pt and the closest point, `location`: distance along the line between start and the closest point.
     * @example
     * var line = turf.lineString([
     *     [-77.031669, 38.878605],
     *     [-77.029609, 38.881946],
     *     [-77.020339, 38.884084],
     *     [-77.025661, 38.885821],
     *     [-77.021884, 38.889563],
     *     [-77.019824, 38.892368]
     * ]);
     * var pt = turf.point([-77.037076, 38.884017]);
     *
     * var snapped = turf.nearestPointOnLine(line, pt, {units: 'miles'});
     *
     * //addToMap
     * var addToMap = [line, pt, snapped];
     * snapped.properties['marker-color'] = '#00f';
     */
    function nearestPointOnLine(lines, pt, options) {
        if (options === void 0) { options = {}; }
        var closestPt = point$2([Infinity, Infinity], {
            dist: Infinity,
        });
        var length = 0.0;
        flattenEach$2(lines, function (line) {
            var coords = getCoords(line);
            for (var i = 0; i < coords.length - 1; i++) {
                //start
                var start = point$2(coords[i]);
                start.properties.dist = distance(pt, start, options);
                //stop
                var stop_1 = point$2(coords[i + 1]);
                stop_1.properties.dist = distance(pt, stop_1, options);
                // sectionLength
                var sectionLength = distance(start, stop_1, options);
                //perpendicular
                var heightDistance = Math.max(start.properties.dist, stop_1.properties.dist);
                var direction = bearing(start, stop_1);
                var perpendicularPt1 = destination(pt, heightDistance, direction + 90, options);
                var perpendicularPt2 = destination(pt, heightDistance, direction - 90, options);
                var intersect = lineIntersect(lineString$1([
                    perpendicularPt1.geometry.coordinates,
                    perpendicularPt2.geometry.coordinates,
                ]), lineString$1([start.geometry.coordinates, stop_1.geometry.coordinates]));
                var intersectPt = null;
                if (intersect.features.length > 0) {
                    intersectPt = intersect.features[0];
                    intersectPt.properties.dist = distance(pt, intersectPt, options);
                    intersectPt.properties.location =
                        length + distance(start, intersectPt, options);
                }
                if (start.properties.dist < closestPt.properties.dist) {
                    closestPt = start;
                    closestPt.properties.index = i;
                    closestPt.properties.location = length;
                }
                if (stop_1.properties.dist < closestPt.properties.dist) {
                    closestPt = stop_1;
                    closestPt.properties.index = i + 1;
                    closestPt.properties.location = length + sectionLength;
                }
                if (intersectPt &&
                    intersectPt.properties.dist < closestPt.properties.dist) {
                    closestPt = intersectPt;
                    closestPt.properties.index = i;
                }
                // update length
                length += sectionLength;
            }
        });
        return closestPt;
    }

    /**
     * Takes two {@link Point|points} and returns a point midway between them.
     * The midpoint is calculated geodesically, meaning the curvature of the earth is taken into account.
     *
     * @name midpoint
     * @param {Coord} point1 first point
     * @param {Coord} point2 second point
     * @returns {Feature<Point>} a point midway between `pt1` and `pt2`
     * @example
     * var point1 = turf.point([144.834823, -37.771257]);
     * var point2 = turf.point([145.14244, -37.830937]);
     *
     * var midpoint = turf.midpoint(point1, point2);
     *
     * //addToMap
     * var addToMap = [point1, point2, midpoint];
     * midpoint.properties['marker-color'] = '#f00';
     */
    function midpoint(point1, point2) {
      var dist = distance(point1, point2);
      var heading = bearing(point1, point2);
      var midpoint = destination(point1, dist / 2, heading);

      return midpoint;
    }

    var MapboxSnap = /** @class */ (function () {
        function MapboxSnap(options) {
            var _a;
            this.status = (_a = options.status) !== null && _a !== void 0 ? _a : false;
            this.map = options.map;
            this.drawing = options.drawing;
            this.options = options.options;
            this.onSnapped = function (fc) {
                if (options.onSnapped !== undefined) {
                    options.onSnapped(fc);
                }
            };
            this.features = {};
            this.snapStatus = false;
            this.snapCoords = [];
            this.radiusInMeters = 0;
            this.addRadiusCircleLayer();
            this.addEvents();
        }
        MapboxSnap.prototype.changeSnappedPoints = function () {
            var drawings = this.drawing.getAll();
            var arr = [];
            for (var i = 0; i < drawings.features.length; i++) {
                var feature = drawings.features[i];
                var id = feature['id'];
                if (this.features[id]) {
                    var snapP = this.features[id].snapPoints;
                    if (this.features['unknow'] !== undefined) {
                        var snapEx = this.features['unknow'].snapPoints;
                        snapP = __assign(__assign({}, snapP), snapEx);
                    }
                    var newFeature = this.doSnap(feature, snapP);
                    arr.push(newFeature);
                }
                else {
                    if (this.features['unknow'] !== undefined) {
                        var snapEx = this.features['unknow'].snapPoints;
                        var newFeature = this.doSnap(feature, snapEx);
                        arr.push(newFeature);
                    }
                    else {
                        arr.push(feature);
                    }
                }
            }
            var fc = { type: 'FeatureCollection', features: arr };
            this.drawing.set(fc);
            if (this.onSnapped) {
                this.onSnapped(fc);
            }
        };
        MapboxSnap.prototype.isPointSnapped = function (p1, p2) {
            var dist = distance(point$2(p1), point$2(p2), { units: 'meters' });
            if (dist < this.radiusInMeters) {
                return true;
            }
            else {
                return false;
            }
        };
        MapboxSnap.prototype.doSnap = function (feature, snaps) {
            switch (feature.geometry.type) {
                case 'Point': {
                    var ownCoords1 = feature.geometry.coordinates;
                    for (var i in snaps) {
                        if (this.isPointSnapped(ownCoords1, snaps[i])) {
                            feature.geometry.coordinates = snaps[i];
                        }
                    }
                    break;
                }
                case 'Polygon': {
                    var ownCoords3 = feature.geometry.coordinates;
                    var newCoords = [];
                    for (var r = 0; r < ownCoords3.length; r++) {
                        var ring = ownCoords3[r];
                        var arr = [];
                        for (var j = 0; j < ring.length; j++) {
                            var coord3 = ring[j];
                            var isOk = false;
                            for (var i in snaps) {
                                if (this.isPointSnapped(coord3, snaps[i])) {
                                    isOk = true;
                                    arr.push(snaps[i]);
                                    break;
                                }
                            }
                            if (isOk == false) {
                                arr.push(coord3);
                            }
                        }
                        newCoords.push(arr);
                    }
                    feature.geometry.coordinates = newCoords;
                    break;
                }
                case 'LineString': {
                    var ownCoords2 = feature.geometry.coordinates;
                    var arr = [];
                    for (var j = 0; j < ownCoords2.length; j++) {
                        var coord = ownCoords2[j];
                        var isOk = false;
                        for (var i in snaps) {
                            if (this.isPointSnapped(coord, snaps[i])) {
                                isOk = true;
                                arr.push(snaps[i]);
                                break;
                            }
                        }
                        if (isOk == false) {
                            arr.push(coord);
                        }
                    }
                    feature.geometry.coordinates = arr;
                }
            }
            return feature;
        };
        MapboxSnap.prototype.getMe = function () {
            return this;
        };
        MapboxSnap.prototype.setStatus = function (s) {
            this.status = s;
        };
        MapboxSnap.prototype.snapToClosestPoint = function (e) {
            if (this.status) {
                var point = e.point;
                var lngLat = this.map.unproject(point);
                var pointAtRadius = [point.x + this.options.radius, point.y];
                var lngLatAtRadius = this.map.unproject(pointAtRadius);
                var radiusInMeters = distance(point$2([lngLat.lng, lngLat.lat]), point$2([lngLatAtRadius.lng, lngLatAtRadius.lat]), { units: 'meters' });
                this.radiusInMeters = radiusInMeters;
                var circle$1 = false;
                var snappedPoint = this.getCloseFeatures(e, radiusInMeters);
                if (snappedPoint) {
                    this.snapStatus = true;
                    this.snapCoords = snappedPoint.coords;
                    circle$1 = circle(snappedPoint.coords, radiusInMeters, { steps: 64, units: 'meters', properties: { color: snappedPoint.color } });
                }
                else {
                    this.snapStatus = false;
                    this.snapCoords = [];
                }
                var circleGeoJSON = featureCollection$2(circle$1 == false ? [] : [circle$1]);
                var source = this.map.getSource('snap-helper-circle');
                if (source) {
                    source.setData(circleGeoJSON);
                }
            }
        };
        MapboxSnap.prototype.addEvents = function () {
            var _this = this;
            this.map.on('mousemove', function (e) {
                _this.snapToClosestPoint(e);
            });
            this.map.on('draw.delete', function (e) {
                setTimeout(function () {
                    _this.changeSnappedPoints();
                }, 100);
            });
            this.map.on('draw.update', function (e) {
                setTimeout(function () {
                    _this.changeSnappedPoints();
                }, 100);
            });
            this.map.on('draw.create', function (e) {
                setTimeout(function () {
                    _this.changeSnappedPoints();
                }, 100);
            });
            this.map.on('draw.selectionchange', function (e) {
                if (e.features.length > 0) {
                    _this.status = true;
                }
                else {
                    setTimeout(function () {
                        _this.changeSnappedPoints();
                    }, 100);
                    _this.status = false;
                }
            });
            this.map.on('draw.modechange', function (e) {
                _this.status = true;
                if (e.mode == 'simple_select') {
                    _this.status = false;
                }
            });
            this.map.on('draw.render', function (e) {
                var source = _this.map.getSource('mapbox-gl-draw-hot');
                if (source) {
                    var data = source._data;
                    if (_this.snapStatus) {
                        var coord = [_this.snapCoords[0], _this.snapCoords[1]];
                        if (data.features.length > 0) {
                            data.features[0].geometry.coordinates = coord;
                        }
                    }
                }
            });
            this.map.on('mouseup', function () {
                _this.drawingSnapCheck();
            });
            this.map.on('click', function () {
                _this.drawingSnapCheck();
            });
        };
        MapboxSnap.prototype.drawingSnapCheck = function () {
            if (this.snapStatus) {
                var source = this.map.getSource('mapbox-gl-draw-hot');
                var coord = [this.snapCoords[0], this.snapCoords[1]];
                var lng = coord[0].toFixed(6);
                var lat = coord[1].toFixed(6);
                var points = {};
                points["".concat(lng, "_").concat(lat)] = coord;
                if (source) {
                    var data = source._data;
                    if (data.features.length > 0) {
                        var f = data.features.find(function (a) { return a.properties.meta == 'feature'; });
                        if (f) {
                            var id = f.properties['id'];
                            if (!this.features[id]) {
                                this.features[id] = {
                                    id: id,
                                    snapPoints: points
                                };
                            }
                            else {
                                this.features[id].snapPoints["".concat(lng, "_").concat(lat)] = coord;
                            }
                        }
                    }
                    else {
                        if (!this.features['unknow']) {
                            this.features['unknow'] = {
                                id: id,
                                snapPoints: points
                            };
                        }
                        else {
                            this.features['unknow'].snapPoints["".concat(lng, "_").concat(lat)] = coord;
                        }
                    }
                }
            }
        };
        MapboxSnap.prototype.searchInVertex = function (feature, mouse, radius) {
            var allCords = coordAll$2(feature);
            var closest = [];
            allCords.map(function (coords) {
                var dist = distance(point$2(coords), point$2([mouse.lng, mouse.lat]), { units: 'meters' });
                if (dist < radius) {
                    closest.push({ coords: coords, dist: dist, color: '#8bc34a' });
                }
            });
            if (closest.length > 0) {
                closest.sort(function (a, b) { return a.dist - b.dist; });
                return closest[0];
            }
        };
        MapboxSnap.prototype.getLines = function (feature, mouse, radius) {
            var lines = [];
            switch (feature.geometry.type) {
                case 'LineString': {
                    lines.push(feature);
                    break;
                }
                case 'MultiLineString': {
                    feature.geometry.coodinates.map(function (c) {
                        lines.push(lineString$1(c));
                    });
                    break;
                }
                case 'Polygon': {
                    var line = polygonToLine(feature.geometry);
                    lines.push(line);
                    break;
                }
                case 'MultiPolygon': {
                    var mlines = polygonToLine(feature.geometry);
                    mlines.coodinates.map(function (c) {
                        lines.push(lineString$1(c));
                    });
                    break;
                }
            }
            return lines;
        };
        MapboxSnap.prototype.searchInMidPoint = function (feature, mouse, radius) {
            var lines = this.getLines(feature, mouse, radius);
            var segments = [];
            lines.map(function (line) {
                segments = segments.concat(lineSegment(line).features);
            });
            var closest = [];
            segments.map(function (seg) {
                var midPoint = midpoint(seg.geometry.coordinates[0], seg.geometry.coordinates[1]);
                var dist = distance(midPoint, point$2([mouse.lng, mouse.lat]), { units: 'meters' });
                if (dist < radius) {
                    closest.push({ coords: midPoint.geometry.coordinates, dist: dist, color: '#03a9f4' });
                }
            });
            if (closest.length > 0) {
                closest.sort(function (a, b) { return a.dist - b.dist; });
                return closest[0];
            }
        };
        MapboxSnap.prototype.searchInEdge = function (feature, mouse, radius) {
            var lines = this.getLines(feature, mouse, radius);
            var closest = [];
            for (var i = 0; i < lines.length; i++) {
                var p = nearestPointOnLine(lines[i], point$2([mouse.lng, mouse.lat]), { units: 'meters' });
                if (p.properties['dist'] !== undefined) {
                    if (p.properties.dist < radius) {
                        closest.push({ coords: p.geometry.coordinates, dist: p.properties.dist, color: '#ff9800' });
                    }
                }
            }
            if (closest.length > 0) {
                closest.sort(function (a, b) { return a.dist - b.dist; });
                return closest[0];
            }
        };
        MapboxSnap.prototype.getCloseFeatures = function (e, radiusInMeters) {
            var features = this.map.queryRenderedFeatures(e.point, {
                layers: this.options.layers
            });
            if (features.length > 0) {
                var snappedPoint;
                var isSnapped = false;
                for (var i = 0; i < features.length; i++) {
                    var mostClose = features[i];
                    var rules = this.options.rules;
                    var mouseCoord = e.lngLat;
                    isSnapped = false;
                    if (rules.indexOf('vertex') !== -1 && snappedPoint == undefined) {
                        snappedPoint = this.searchInVertex(mostClose, mouseCoord, radiusInMeters);
                        if (snappedPoint) {
                            isSnapped = true;
                            break;
                        }
                    }
                    if (rules.indexOf('midpoint') !== -1 && snappedPoint == undefined) {
                        snappedPoint = this.searchInMidPoint(mostClose, mouseCoord, radiusInMeters);
                        if (snappedPoint) {
                            isSnapped = true;
                            break;
                        }
                    }
                    if (rules.indexOf('edge') !== -1 && snappedPoint == undefined) {
                        snappedPoint = this.searchInEdge(mostClose, mouseCoord, radiusInMeters);
                        if (snappedPoint) {
                            isSnapped = true;
                            break;
                        }
                    }
                }
                if (isSnapped) {
                    return snappedPoint;
                }
                else {
                    return false;
                }
            }
            else {
                return false;
            }
        };
        MapboxSnap.prototype.addRadiusCircleLayer = function () {
            this.map.addSource('snap-helper-circle', { type: 'geojson', data: { type: 'FeatureCollection', features: [] } });
            this.map.addLayer({
                id: "snap-helper-circle",
                type: 'fill',
                source: "snap-helper-circle",
                paint: {
                    'fill-color': ['get', 'color'],
                    'fill-opacity': 0.6
                }
            });
        };
        return MapboxSnap;
    }());

    return MapboxSnap;

}));
