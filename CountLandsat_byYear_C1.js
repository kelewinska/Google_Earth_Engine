/*
#+
# :AUTHOR: Katarzyna Ewa Lewinska [lewinska@wisc.edu]
# :DATE: 14 October 2019
#
# :Description: Scripts provides a per-pixel overview of cloud- and snow-free Landsat observations (L4-8) for a given time period.
#               Data are sumarized on yearly bases for all months (Jan-Dec), or for a subset of several consecutive months
#               (e.g. Apr-Sep). The scripts assumes Jan to be the first, and Dec to be the last month.
#               All layers are by default OFF and need to be toggled ON to be seen.
#
# :Input:       startY - [integer YYYY] a starting year for the inquiry
#               endY - [integer YYYY] an ending year for the inquiry
#               startM - [integer M or MM] a starting month for the inquiry; min=1, max=12
#               endM - [integerM or MM] an endingmonth for the inquiry; min=1, max=12
#               AOI - area of interest defined either as ee.GeometryPolygon or hand drawn geometry feature
#
# :Output:      Visual Layers (one year - one layer; separately scalled and collored)
#               multilayer image exported to Drive (one year - one band)
#
# :Updates:     2019-10-15:   No duplicated based of He Yin's code
#               2021-10-24:   'retired' as Collection 1 solution
#
# :2Do:
#
# :Disclaimer:  The author of this code accepts no responsibility for errors or omissions in this work
#               and shall not be liable for any damage caused by these.
#-
*/



// ### INPUTS ### \\
var startY = 1984
var endY = 1990

var startM = 1
var endM = 12


// AOI definition
var AOI =  ee.Geometry.Polygon(
        [[[-89.73775572597816, 43.444981295302384], [-89.73775572597816, 42.80350046836383],
          [-88.55123228847816, 42.80350046836383],  [-88.55123228847816, 43.444981295302384]]]);


// ### ENVIRONMENT ### \\


// ### DATA ### \\
var l8 = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR').filter(ee.Filter.calendarRange(startY, endY,'year'))
                                                     .filter(ee.Filter.calendarRange(startM, endM, 'month'))
                                                     .filterBounds(AOI);
var l7 = ee.ImageCollection('LANDSAT/LE07/C01/T1_SR').filter(ee.Filter.calendarRange(startY, endY,'year'))
                                                     .filter(ee.Filter.calendarRange(startM, endM, 'month'))
                                                     .filterBounds(AOI);
var l5 = ee.ImageCollection('LANDSAT/LT05/C01/T1_SR').filter(ee.Filter.calendarRange(startY, endY,'year'))
                                                     .filter(ee.Filter.calendarRange(startM, endM, 'month'))
                                                     .filterBounds(AOI);
var l4 = ee.ImageCollection('LANDSAT/LT04/C01/T1_SR').filter(ee.Filter.calendarRange(startY, endY,'year'))
                                                     .filter(ee.Filter.calendarRange(startM, endM, 'month'))
                                                     .filterBounds(AOI);


// ### FUNCTIONALITY ### \\

//cloud mask for L4-7
var cloudMaskL457 = function(in_image) {
  var qa = in_image.select('pixel_qa');
  var cloud = qa.bitwiseAnd(1 << 5)
          .and(qa.bitwiseAnd(1 << 7))
          .or(qa.bitwiseAnd(1 << 3))
          .or(qa.bitwiseAnd(1 << 4));
  var mask2 = in_image.mask().reduce(ee.Reducer.min());  // Remove edge pixels that don't occur in all bands
  var mask3 = in_image.select('B1').lt(20000);
  return in_image.updateMask(cloud.not()).updateMask(mask2).updateMask(mask3);
  };

// OLI mask
var cloudMaskL8 = function(in_image) {
  var qa = in_image.select('pixel_qa');
  var cloud = qa.bitwiseAnd(1 << 5)
                .and(qa.bitwiseAnd(1 << 6).or(qa.bitwiseAnd(1 << 7)))
                .or(qa.bitwiseAnd(1 << 4))
                .or(qa.bitwiseAnd(1 << 3))
                .or(qa.bitwiseAnd(1 << 8).and(qa.bitwiseAnd(1 << 9)));
                // .or(qa.bitwiseAnd(1 << 6)) //.and(qa.bitwiseAnd(1 << 7))) // in fact, while taking all three cloud confidence levels, this statement is equal to (1 << 5)
                // .or(qa.bitwiseAnd(1 << 4))       // mask snow and ice
                // .or(qa.bitwiseAnd(1 << 3))       // mask shadows
                // .or(qa.bitwiseAnd(1 << 8).and(qa.bitwiseAnd(1 << 9)));     // cirrus confidence - if medium of high then mask    .or(qa.bitwiseAnd(1 << 9))
  var mask2 = in_image.mask().reduce(ee.Reducer.min());  // Remove edge pixels that don't occur in all bands
  return in_image.updateMask(cloud.not()).updateMask(mask2);
};

var DateBand_F = function(image){
  var date = ee.Number.parse(image.date().format('YYYYMMdd'));
  var DateBand = ee.Image.constant(date).uint32().rename('date')
  DateBand = DateBand.updateMask(image.select('B4').mask())
      image = image.addBands(DateBand.updateMask(image.select('B4')));
  return image.select('date').selfMask()
};

var years_F = function(y){
  var y_count = col_B4.filter(ee.Filter.calendarRange(y,y,'year')).count().rename(ee.String(y).slice(0,4));
  return ee.Image(y_count);
};

var years_noDup_F = function(y){
  var y_count_dis = col_B4.filter(ee.Filter.calendarRange(y,y,'year')).reduce(ee.Reducer.countDistinct()).rename('no-duplicates');
  return ee.Image(y_count_dis);
};

// ### CODE ### \\

// mask Landsat collections and select NIR band
var l4M_B4 = l4.map(cloudMaskL457).select("B4");
var l5M_B4 = l5.map(cloudMaskL457).select("B4");
var l7M_B4 = l7.map(cloudMaskL457).select("B4");
var l8M_B4 = l8.map(cloudMaskL8).select("B5").map(function(im){ return im.select("B5").rename('B4')});


// merge all collections
var col_B4 = l4M_B4.merge(l5M_B4).merge(l7M_B4).merge(l8M_B4);

var years = ee.List.sequence(startY, endY);
var yearsN = ee.List.sequence(startY-1, endY-1);


var yearly_count = ee.ImageCollection(years.map(years_noDup_F));

var years2bands = yearly_count.toBands().rename(years.map(function(y){
                                                return ee.String(y).slice(0,4);
                                                }));

    years2bands = years2bands.clip(AOI).int();

Map.centerObject(AOI, 7);


// Add Visual Layers
for (var no = startY; no <= endY; no++){
  var vis = {min:0, max: 28};
  Map.addLayer(years2bands.select(no-startY), vis, 'year'+no, 0);
}

Map.setOptions('HYBRID');


// /// ### Export 2 Asset ### \\\
Export.image.toDrive({
  image:years2bands,
  description: 'CountLandsatByYear'+startY+'_'+endY,
  scale: 30,
  region: AOI,
  maxPixels:1e13,
});


/*
*/

// this is the end of the code \\
