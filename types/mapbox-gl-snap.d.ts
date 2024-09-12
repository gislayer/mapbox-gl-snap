declare module 'mapbox-gl-snap' {
  import {
    Map
  } from 'mapbox-gl';
  import MapboxDraw from '@mapbox/mapbox-gl-draw';
  import {
    FeatureCollection,
    Geometry
  } from 'geojson';

  export interface MapboxSnapOptions {
    layers: string[];
    rules: ('vertex' | 'edge' | 'middle')[];
    radius: number;
  }

  export interface MapboxSnapProps {
    map: Map;
    drawing: MapboxDraw;
    options: MapboxSnapOptions;
    status ? : boolean;
    onSnapped ? : (features: FeatureCollection < Geometry > ) => void;
  }

  class MapboxSnap {
    constructor(options: MapboxSnapProps);
    setStatus(status: boolean): void;
    getMe(): MapboxSnap;
  }

  export default MapboxSnap;
}