class MapboxSnap {
  constructor(options) {
    this.status = options.status ?? false;
    this.map = options.map;
    this.drawing = options.drawing;
    this.options = options.options;
    this.onSnapped = (fc) => {
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

  changeSnappedPoints() {
    var drawings = this.drawing.getAll();
    var arr = [];
    for (var i = 0; i < drawings.features.length; i++) {
      var feature = drawings.features[i];
      var id = feature['id'];
      if (this.features[id]) {
        var snapP = this.features[id].snapPoints;
        if (this.features['unknow'] !== undefined) {
          var snapEx = this.features['unknow'].snapPoints;
          snapP = { ...snapP, ...snapEx };
        }
        var newFeature = this.doSnap(feature, snapP);
        arr.push(newFeature);
      } else {
        if (this.features['unknow'] !== undefined) {
          var snapEx = this.features['unknow'].snapPoints;
          var newFeature = this.doSnap(feature, snapEx);
          arr.push(newFeature);
        } else {
          arr.push(feature);
        }
      }
    }
    var fc = { type: 'FeatureCollection', features: arr };
    this.drawing.set(fc);
    if (this.onSnapped) {
      this.onSnapped(fc);
    }
  }

  isPointSnapped(p1, p2) {
    var dist = turf.distance(turf.point(p1), turf.point(p2), { units: 'meters' });
    return dist < this.radiusInMeters;
  }

  doSnap(feature, snaps) {
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
            if (!isOk) {
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
          if (!isOk) {
            arr.push(coord);
          }
        }
        feature.geometry.coordinates = arr;
      }
    }
    return feature;
  }

  getMe() {
    return this;
  }

  setStatus(s) {
    this.status = s;
  }

  snapToClosestPoint(e) {
    if (this.status) {
      const point = e.point;
      const lngLat = this.map.unproject(point);
      const pointAtRadius = [point.x + this.options.radius, point.y];
      const lngLatAtRadius = this.map.unproject(pointAtRadius);
      const radiusInMeters = turf.distance(
          turf.point([lngLat.lng, lngLat.lat]),
          turf.point([lngLatAtRadius.lng, lngLatAtRadius.lat]),
        { units: 'meters' }
      );
      this.radiusInMeters = radiusInMeters;
      var circle = false;
      const snappedPoint = this.getCloseFeatures(e, radiusInMeters);
      if (snappedPoint) {
        this.snapStatus = true;
        this.snapCoords = snappedPoint.coords;
        circle = turf.circle(snappedPoint.coords, radiusInMeters, { steps: 64, units: 'meters', properties: { color: snappedPoint.color } });
      } else {
        this.snapStatus = false;
        this.snapCoords = [];
      }
      const circleGeoJSON = turf.featureCollection(circle ? [circle] : []);
      let source = this.map.getSource('snap-helper-circle');
      if (source) {
        source.setData(circleGeoJSON);
      }
    }
  }

  addEvents() {
    this.map.on('mousemove', (e) => {
      this.snapToClosestPoint(e);
    });

    this.map.on('draw.selectionchange', (e) => {
      if (e.features.length > 0) {
        this.status = true;
      } else {
        setTimeout(() => {
          this.changeSnappedPoints();
        }, 100);
        this.status = false;
      }
    });

    this.map.on('draw.modechange', (e) => {
      this.status = true;
      if (e.mode === 'simple_select') {
        this.status = false;
      }
    });

    this.map.on('draw.render', (e) => {
      let source = this.map.getSource('mapbox-gl-draw-hot');
      if (source) {
        var data = source._data;
        if (this.snapStatus) {
          var coord = [this.snapCoords[0], this.snapCoords[1]];
          if (data.features.length > 0) {
            data.features[0].geometry.coordinates = coord;
          }
        }
      }
    });

    this.map.on('mouseup', () => {
      this.drawingSnapCheck();
    });

    this.map.on('click', () => {
      this.drawingSnapCheck();
    });
  }

  drawingSnapCheck() {
    if (this.snapStatus) {
      let source = this.map.getSource('mapbox-gl-draw-hot');
      var coord = [this.snapCoords[0], this.snapCoords[1]];
      var lng = coord[0].toFixed(6);
      var lat = coord[1].toFixed(6);
      var points = {};
      points[`${lng}_${lat}`] = coord;
      if (source) {
        var data = source._data;
        if (data.features.length > 0) {
          var f = data.features.find((a) => a.properties.meta === 'feature');
          if (f) {
            var id = f.properties['id'];
            if (!this.features[id]) {
              this.features[id] = {
                id: id,
                snapPoints: points,
              };
            } else {
              this.features[id].snapPoints[`${lng}_${lat}`] = coord;
            }
          }
        } else {
          if (!this.features['unknow']) {
            this.features['unknow'] = {
              id: id,
              snapPoints: points,
            };
          } else {
            this.features['unknow'].snapPoints[`${lng}_${lat}`] = coord;
          }
        }
      }
    }
  }

  searchInVertex(feature, mouse, radius) {
    const allCords = turf.coordAll(feature);
    const closest = [];
    allCords.map((coords) => {
      var dist = turf.distance(turf.point(coords), turf.point([mouse.lng, mouse.lat]), { units: 'meters' });
      if (dist < radius) {
        closest.push({ coords: coords, dist: dist, color: '#8bc34a' });
      }
    });
    if (closest.length > 0) {
      closest.sort((a, b) => a.dist - b.dist);
      return closest[0];
    }
  }

  getLines(feature, mouse, radius) {
    var lines = [];
    switch (feature.geometry.type) {
      case 'LineString': {
        lines.push(feature);
        break;
      }
      case 'MultiLineString': {
        feature.geometry.coodinates.map((c) => {
          lines.push(turf.lineString(c));
        });
        break;
      }
      case 'Polygon': {
        var line = turf.polygonToLine(feature.geometry);
        lines.push(line);
        break;
      }
      case 'MultiPolygon': {
        const mlines = turf.polygonToLine(feature.geometry);
        mlines.coodinates.map((c) => {
          lines.push(turf.lineString(c));
        });
        break;
      }
    }
    return lines;
  }

  searchInMidPoint(feature, mouse, radius) {
    var lines = this.getLines(feature, mouse, radius);
    var segments = [];
    lines.map((line) => {
      segments = segments.concat(turf.lineSegment(line).features);
    });
    var closest = [];
    segments.map((seg) => {
      var midPoint = turf.midpoint(seg.geometry.coordinates[0], seg.geometry.coordinates[1]);
      var dist = turf.distance(midPoint, turf.point([mouse.lng, mouse.lat]), { units: 'meters' });
      if (dist < radius) {
        closest.push({ coords: midPoint.geometry.coordinates, dist: dist, color: '#03a9f4' });
      }
    });
    if (closest.length > 0) {
      closest.sort((a, b) => a.dist - b.dist);
      return closest[0];
    }
  }

  searchInEdge(feature, mouse, radius) {
    var lines = this.getLines(feature, mouse, radius);
    var closest = [];
    for (var i = 0; i < lines.length; i++) {
      var p = turf.nearestPointOnLine(lines[i], turf.point([mouse.lng, mouse.lat]), { units: 'meters' });
      if (p.properties['dist'] !== undefined) {
        if (p.properties.dist < radius) {
          closest.push({ coords: p.geometry.coordinates, dist: p.properties.dist, color: '#ff9800' });
        }
      }
    }
    if (closest.length > 0) {
      closest.sort((a, b) => a.dist - b.dist);
      return closest[0];
    }
  }

  getCloseFeatures(e, radiusInMeters) {
    const features = this.map.queryRenderedFeatures(e.point, {
      layers: this.options.layers,
    });
    if (features.length > 0) {
      var snappedPoint;
      for(var i=0; i<features.length; i++){
        var mostClose = features[i];
        var rules = this.options.rules;
        var mouseCoord = e.lngLat;
        var isSnapped = false;
  
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
      }else{
        return false;
      }
    } else {
      return false;
    }
  }

  addRadiusCircleLayer() {
    this.map.addSource('snap-helper-circle', { type: 'geojson', data: { type: 'FeatureCollection', features: [] } });
    this.map.addLayer({
      id: `snap-helper-circle`,
      type: 'fill',
      source: `snap-helper-circle`,
      paint: {
        'fill-color': ['get', 'color'],
        'fill-opacity': 0.6,
      },
    });
  }
}
