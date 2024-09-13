const commonjs = require('@rollup/plugin-commonjs');
const resolve = require('@rollup/plugin-node-resolve');
const typescript = require('@rollup/plugin-typescript');

exports.default = {
  input: 'src/MapboxSnap.ts', // Giriş dosyanız
  output: {
    file: 'dist/purejs/mapbox-gl-snap.js', // Çıktı dosyası
    format: 'umd',
    name: 'MapboxSnap'
  },
  plugins: [resolve(), commonjs(), typescript()]
};
