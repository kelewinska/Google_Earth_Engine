/*
#+
# :AUTHOR: Katarzyna Ewa Lewinska [lewinska@wisc.edu]
# :DATE: 5 February 2019
#
# :Description: Script allows on comparison between spectral unmixing executed with unmix and matrixSolve function.
#               Due to limitations the latter approach shade endmember is introduced as [1,1,1,1,1,1] (an endmember cannot be [0,0,0,0,0,0]).
#               The code comes with the solveMartix solution being commented out.
#               In the second part of a code, a comparison among results of the unmix approach run for 4 combinations of
#               sumToOne and nonNegative constrains is shown. In addition to endmembers' abundances, rmse for each solution is calculated
#               through a 'reversed' approach. (A reversed image is constructed based on endmembers and their obtained abundances.
#               Band-wise differences between the original image and the reversed image are used to calculate rmse for each model).
#               The script allows on visualization on all results in a quadruple display, where each considered model
#               is presented in a different map window.
#               In the current example the code runs for a single Landsat 8 OLI data, converting it to the ETM+ space. The code can be
#               adjusted to work for TM and ETM+ scenes, as well as long time series of Landsat data. Implemented cloud masking is optional.
#               A scene to be used is identified as a first image in the defined search window overlapping with a used-defined geometry point.
#               Endmembers definitions are based on [10.1016/j.rse.2005.07.013].
#
# :Input:       point - geometry point to identify a corresponding Landsat tile
#               start - a string YYYY-MM-DD defining the first day of the scene search window
#               end - a string YYYY-MM-DD defining the last day of the scene search window
#               soil, gv, npv, shade - endmembers introduced as separate lists
#
# :Output:      Dynamic quadruple linked display. Values can be inspected in the most-left map view.
#
# :Updates:      2019-02-05:  quadruple linked display;  RMSE
#               2021-10-24:   'retired' as Collection 1 solution
#
# :2Do:
#
# :Disclaimer:  The author of this code accepts no responsibility for errors or omissions in this work
#               and shall not be liable for any damage caused by these.
#-
*/


// ### INPUTS ### \\

var point = ee.Geometry.Point([46.267220603552914, 41.37497420893961]);

var start = '2018-07-20';
var end = '2018-08-01';

var soil = ee.Array([900, 1000, 1600, 3000, 2900, 4800]);
var gv = ee.Array([500, 900, 400, 6100, 3000, 1000]);
var npv = ee.Array([1400, 1700, 2200, 3000, 5500, 3000]);
var shade = ee.Array([1, 1, 1, 1, 1, 1]);

var vis = {"opacity":1,"bands":["RMSE"],"max":1500.0,"gamma":1};


// ### ENVIRONMENT ### \\

// ETM+ cloud mask
var cloudMaskL8 = function(in_image) {
  var qa = in_image.select('pixel_qa');
  var cloud = qa.bitwiseAnd(1 << 5)              // mask clouds of low, medium and high confidence level
                .or(qa.bitwiseAnd(1 << 6).and(qa.bitwiseAnd(1 << 7))) // this statement is equal to (1 << 5)
                .or(qa.bitwiseAnd(1 << 4))       // mask snow and ice
                .or(qa.bitwiseAnd(1 << 3))       // mask shadows
                .or(qa.bitwiseAnd(1 << 8).and(qa.bitwiseAnd(1 << 9)));     // cirrus confidence - if medium of high then mask    .or(qa.bitwiseAnd(1 << 9))
  var mask2 = in_image.mask().reduce(ee.Reducer.min());  // Remove edge pixels that don't occur in all bands
  return in_image.updateMask(cloud.not()).updateMask(mask2);
};

// OLI to ETM function
  // after doi:10.1016/j.rse.2015.12.024
var OLI2ETMf = function(image) {
  var blue = image.expression('0.0183 + 0.8850 * B2', {'B2':image.select('B2')});
  var green = image.expression('0.0123 + 0.9317 * B3', {'B3':image.select('B3')});
  var red = image.expression('0.0123 + 0.9372 * B4', {'B4':image.select('B4')});
  var nir = image.expression('0.0448 + 0.8339 * B5', {'B5':image.select('B5')});
  var swir1 = image.expression('0.0306 + 0.8639 * B6', {'B6':image.select('B6')});
  var swir2 = image.expression('0.0116 + 0.9165 * B7', {'B7':image.select('B7')});

  var temp = ee.Image(image).addBands(blue).addBands(green).addBands(red).addBands(nir).addBands(swir1).addBands(swir2);

  return temp.select(['constant', 'constant_1', 'constant_2', 'constant_3', 'constant_4', 'constant_5' ],
                   ['B1','B2','B3','B4','B5','B7']).int16();
};


// ### CODE ### \\

var l8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR');

var ts_l8 = ee.ImageCollection(l8.filterDate(start, end).filterBounds(point));
    // ts_l8 = ts_l8.map(cloudMaskL8);        // cloud mask
    ts_l8 = ts_l8.map(OLI2ETMf);              // OLI -> ETM+ normalization

var image = ts_l8.first();
print(image);
var bands = image.bandNames();




// // unmixing through matrixSolve()

// var end_arr = ee.Array.cat([soil, gv, npv, shade],1);
// print('end_arr', end_arr, end_arr.length());                //[6x5]

// // Make an Array Image with a 2-D Array per pixel, 6x1.
// var arrayImage = image.toArray().toArray(1);                //[6x1]
// print('arrayImage', arrayImage);

// //  [4x1]              [6x4]                 [6x1]  -> [6x4]x[4x1]=[6x1]
// var unmixed = ee.Image(end_arr).matrixSolve(arrayImage);

// // var unmixedImage = unmixed.arrayProject([0]).arrayFlatten([['soil', 'gv', 'npv', 'shade']]); // 0 image axis

// print('unmixedImage', unmixedImage);
// Map.addLayer(unmixedImage, {}, 'matrixSolve');





// unmix with unmix()
var endmembers = ee.Dictionary({
    'soil': [2000, 3000, 3400, 5800, 6000, 5800],
    'gv': [500, 900, 400, 6100, 3000, 1000],
    'npv': [1400, 1700, 2200, 3000, 5500, 3000],
    'shade': [1, 1, 1, 1, 1, 1],
      });
var u_soil = endmembers.get('soil');
var u_gv = endmembers.get('gv');
var u_npv =endmembers.get('npv');
var u_shade = endmembers.get('shade');



// sumToOne = F and nonNegative=F
var unmixFF = image.unmix([u_soil, u_gv, u_npv, u_shade], false, false).rename(['soil','gv','npv','shade']);

  // sumToOne = T and nonNegative=T
var unmixTT = image.unmix([u_soil, u_gv, u_npv, u_shade], true, true).rename(['soil','gv','npv','shade']);

  // sumToOne = T and nonNegative=false
var unmixTF = image.unmix([u_soil, u_gv, u_npv, u_shade], true, false).rename(['soil','gv','npv','shade']);

  // sumToOne = F and nonNegative=T
var unmixFT = image.unmix([u_soil, u_gv, u_npv, u_shade], false, true).rename(['soil','gv','npv','shade']);



// ### ERROR TERMS ### \\

var end_arr = ee.Array.cat([soil, gv, npv, shade],1);

// reverse the image from spectra and endmembers
var unmixArFF = unmixFF.toArray().toArray(1);  //[5x1]
var revImFF = ee.Image(end_arr).matrixMultiply(unmixArFF);    //  [6x5] x [5x1] -> [6x1]
var revImFF = revImFF.arrayProject([0]).arrayFlatten([bands]);

var unmixArTT = unmixTT.toArray().toArray(1);  //[5x1]
var revImTT = ee.Image(end_arr).matrixMultiply(unmixArTT);    //  [6x5] x [5x1] -> [6x1]
var revImTT = revImTT.arrayProject([0]).arrayFlatten([bands]);

var unmixArTF = unmixTF.toArray().toArray(1);  //[5x1]
var revImTF = ee.Image(end_arr).matrixMultiply(unmixArTF);    //  [6x5] x [5x1] -> [6x1]
var revImTF = revImTF.arrayProject([0]).arrayFlatten([bands]);

var unmixArFT = unmixFT.toArray().toArray(1);  //[5x1]
var revImFT = ee.Image(end_arr).matrixMultiply(unmixArFT);    //  [6x5] x [5x1] -> [6x1]
var revImFT = revImFT.arrayProject([0]).arrayFlatten([bands]);

var difFF = image.subtract(revImFF);
var difTT = image.subtract(revImTT);
var difTF = image.subtract(revImTF);
var difFT = image.subtract(revImFT);


var difFFsq = ((difFF.select('B1').pow(ee.Image(2)))
         .add(difFF.select('B2').pow(ee.Image(2)))
         .add(difFF.select('B3').pow(ee.Image(2)))
         .add(difFF.select('B4').pow(ee.Image(2)))
         .add(difFF.select('B5').pow(ee.Image(2)))
         .add(difFF.select('B7').pow(ee.Image(2))).divide(ee.Image(6))).sqrt().rename('RMSE');
var difTTsq = ((difTT.select('B1').pow(ee.Image(2)))
         .add(difTT.select('B2').pow(ee.Image(2)))
         .add(difTT.select('B3').pow(ee.Image(2)))
         .add(difTT.select('B4').pow(ee.Image(2)))
         .add(difTT.select('B5').pow(ee.Image(2)))
         .add(difTT.select('B7').pow(ee.Image(2))).divide(ee.Image(6))).sqrt().rename('RMSE');
var difTFsq = ((difTF.select('B1').pow(ee.Image(2)))
         .add(difTF.select('B2').pow(ee.Image(2)))
         .add(difTF.select('B3').pow(ee.Image(2)))
         .add(difTF.select('B4').pow(ee.Image(2)))
         .add(difTF.select('B5').pow(ee.Image(2)))
         .add(difTF.select('B7').pow(ee.Image(2))).divide(ee.Image(6))).sqrt().rename('RMSE');
var difFTsq = ((difFT.select('B1').pow(ee.Image(2)))
         .add(difFT.select('B2').pow(ee.Image(2)))
         .add(difFT.select('B3').pow(ee.Image(2)))
         .add(difFT.select('B4').pow(ee.Image(2)))
         .add(difFT.select('B5').pow(ee.Image(2)))
         .add(difFT.select('B7').pow(ee.Image(2))).divide(ee.Image(6))).sqrt().rename('RMSE');




/// ### MAPS ### \\\

  var long1= point.getInfo()['coordinates'][0];
  var lat1=point.getInfo()['coordinates'][1];

  function initMap(map) {
  map.setCenter(long1,lat1, 9);
  }

// Initialize Map
  initMap(Map);

  function createMap(title) {
  var map = ui.Map();
  ui.Label(title, {position:'bottom-center'});
  map.add(title);
  return map;
  }

  function getMapSize() {
  var scale = Map.getScale();
  var bounds = ee.Geometry(Map.getBounds(true)).coordinates().get(0).getInfo();

  var ll = bounds[0];
  var ur = bounds[2];
  var width = (ur[0] - ll[0]) / scale;
  var height = (ur[1] - ll[1]) / scale;

  return { w: Math.floor(width), h: Math.floor(height) };
  }

  var height = getMapSize().h;

var map0 = ui.Map();
var map1 = ui.Map();
var map2 = ui.Map();
var map3 = ui.Map();

    map0.addLayer(difFTsq, vis, 'FT Error', false);
    map0.addLayer(difTFsq, vis, 'TF Error', false);
    map0.addLayer(difTTsq, vis, 'TT Error', false);
    map0.addLayer(difFFsq, vis, 'FF Error', false);
    map0.addLayer(unmixFT, {}, 'FT', false);
    map0.addLayer(unmixTF, {}, 'TF', false);
    map0.addLayer(unmixTT, {}, 'TT', false);
    map0.addLayer(unmixFF, {}, 'FF');

    map1.addLayer(difTTsq, vis, 'TT Error');
    map1.addLayer(unmixTT, {}, 'TT');
    map2.addLayer(difTFsq, vis, 'TF Error');
    map2.addLayer(unmixTF, {}, 'TF');
    map3.addLayer(difFTsq, vis, 'FT Error');
    map3.addLayer(unmixFT, {}, 'FT');

    map0.add(ui.Label('SumToOne= F  nonNegative= F', {position:'bottom-center'}));
    map1.add(ui.Label('SumToOne= T  nonNegative= T', {position:'bottom-center'}));
    map2.add(ui.Label('SumToOne= T  nonNegative= F', {position:'bottom-center'}));
    map3.add(ui.Label('SumToOne= F  nonNegative= T', {position:'bottom-center'}));

var maps = [map0, map1, map2, map3];

// Create a panel with vertical flow layout.
  var panel = ui.Panel({
  layout: ui.Panel.Layout.flow('horizontal'),
  style: {width: '100vw', height: height + '300px'}
  });

  var linker = ui.Map.Linker(maps);

  maps.map(function(map) {
  initMap(map);
  // map.setOptions('HYBRID');
  panel.add(map);
  });

  ui.root.clear();
  ui.root.add(maps[0]);
  ui.root.add(maps[1]);
  ui.root.add(maps[2]);
  ui.root.add(maps[3]);

  var linker = ui.Map.Linker(maps);

// this is the end of the code \\
