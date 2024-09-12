/// <reference types="mapbox__mapbox-gl-draw" />
import { MapMouseEvent, LngLat } from 'mapbox-gl';
import MapboxDraw from "@mapbox/mapbox-gl-draw";
import { Feature, FeatureCollection, Geometry } from 'geojson';
interface SpatialResult {
    coords: number[];
    dist: number;
    color: string;
}
export interface MapboxSnapOptions {
    layers: string[];
    rules: ('vertex' | 'edge' | 'middle')[];
    radius: number;
}
export interface MapboxSnapProps {
    map: mapboxgl.Map;
    drawing: MapboxDraw;
    options: MapboxSnapOptions;
    status?: boolean;
    onSnapped?: Function;
}
interface SnapSettings {
    id: string;
    snapPoints: Record<string, number[]>;
}
declare class MapboxSnap {
    map: mapboxgl.Map;
    drawing: MapboxDraw;
    options: MapboxSnapOptions;
    snapStatus: boolean;
    snapCoords: number[];
    status: boolean;
    features: Record<string, SnapSettings>;
    radiusInMeters: number;
    onSnapped: (features: FeatureCollection<Geometry>) => void;
    constructor(options: MapboxSnapProps);
    changeSnappedPoints(): void;
    isPointSnapped(p1: number[], p2: number[]): boolean;
    doSnap(feature: Feature, snaps: Record<string, number[]>): Feature<Geometry, import("geojson").GeoJsonProperties>;
    getMe(): this;
    setStatus(s: boolean): void;
    snapToClosestPoint(e: MapMouseEvent): void;
    addEvents(): void;
    drawingSnapCheck(): void;
    searchInVertex(feature: any, mouse: LngLat, radius: number): SpatialResult | undefined;
    getLines(feature: any, mouse: LngLat, radius: number): any[];
    searchInMidPoind(feature: any, mouse: LngLat, radius: number): SpatialResult | undefined;
    searchInEdge(feature: any, mouse: LngLat, radius: number): SpatialResult | undefined;
    getCloseFeatures(e: MapMouseEvent, radiusInMeters: number): any;
    addRadiusCircleLayer(): void;
}
export default MapboxSnap;
