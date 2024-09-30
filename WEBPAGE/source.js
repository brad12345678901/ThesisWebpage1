var map = L.map('map', {
    zoomControl:true, maxZoom:28, minZoom:1
}).fitBounds([[14.66996750972913,121.0946753134224],[14.690671168315529,121.1325971661867]]);
var hash = new L.Hash(map);
map.attributionControl.setPrefix('<a href="https://github.com/tomchadwin/qgis2web" target="_blank">qgis2web</a> &middot; <a href="https://leafletjs.com" title="A JS library for interactive maps">Leaflet</a> &middot; <a href="https://qgis.org">QGIS</a>');
var autolinker = new Autolinker({truncate: {length: 30, location: 'smart'}});

function removeEmptyRowsFromPopupContent(content, feature) {
    var tempDiv = document.createElement('div');
    tempDiv.innerHTML = content;
    var rows = tempDiv.querySelectorAll('tr');
    for (var i = 0; i < rows.length; i++) {
        var td = rows[i].querySelector('td.visible-with-data');
        var key = td ? td.id : '';
        if (td && td.classList.contains('visible-with-data') && feature.properties[key] == null) {
            rows[i].parentNode.removeChild(rows[i]);
        }
    }
    return tempDiv.innerHTML;
}

document.querySelector(".leaflet-popup-pane").addEventListener("load", function(event) {
    var tagName = event.target.tagName,
        popup = map._popup;
    // Also check if flag is already set.
    if (tagName === "IMG" && popup && !popup._updated) {
        popup._updated = true; // Set flag to prevent looping.
        popup.update();
    }
}, true);

var bounds_group = new L.featureGroup([]);

function setBounds() {
}

map.createPane('pane_OpenStreetMap_0');
map.getPane('pane_OpenStreetMap_0').style.zIndex = 400;

var layer_OpenStreetMap_0 = L.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png', {
    pane: 'pane_OpenStreetMap_0',
    opacity: 1.0,
    attribution: '',
    minZoom: 1,
    maxZoom: 28,
    minNativeZoom: 0,
    maxNativeZoom: 19
});

layer_OpenStreetMap_0;
map.addLayer(layer_OpenStreetMap_0);
function pop_RoadLayerV3_1(feature, layer) {
    var popupContent = '<table>\
            <tr>\
                <td colspan="2"><strong>name</strong><br />' + (feature.properties['name'] !== null ? autolinker.link(feature.properties['name'].toLocaleString()) : '') + '</td>\
            </tr>\
            <tr>\
                <td colspan="2"><strong>highway</strong><br />' + (feature.properties['highway'] !== null ? autolinker.link(feature.properties['highway'].toLocaleString()) : '') + '</td>\
            </tr>\
            <tr>\
                <th scope="row">access</th>\
                <td class="visible-with-data" id="access">' + (feature.properties['access'] !== null ? autolinker.link(feature.properties['access'].toLocaleString()) : '') + '</td>\
            </tr>\
            <tr>\
                <th scope="row">oneway</th>\
                <td class="visible-with-data" id="oneway">' + (feature.properties['oneway'] !== null ? autolinker.link(feature.properties['oneway'].toLocaleString()) : '') + '</td>\
            </tr>\
        </table>';
    layer.bindPopup(popupContent, {maxHeight: 400});
    var popup = layer.getPopup();
    var content = popup.getContent();
    var updatedContent = removeEmptyRowsFromPopupContent(content, feature);
    popup.setContent(updatedContent);
}

function style_RoadLayerV3_1_0() {
    return {
        pane: 'pane_RoadLayerV3_1',
        opacity: 1,
        color: 'rgba(35,35,35,1)',
        dashArray: '',
        lineCap: 'round',
        lineJoin: 'bevel',
        weight: 3.0,
        fillOpacity: 0,
        interactive: false,
    }
}
function style_RoadLayerV3_1_1() {
    return {
        pane: 'pane_RoadLayerV3_1',
        opacity: 1,
        color: 'rgba(255,255,255,1)',
        dashArray: '',
        lineCap: 'square',
        lineJoin: 'bevel',
        weight: 1.0,
        fillOpacity: 0,
        interactive: false,
    }
}
map.createPane('pane_RoadLayerV3_1');
map.getPane('pane_RoadLayerV3_1').style.zIndex = 401;
map.getPane('pane_RoadLayerV3_1').style['mix-blend-mode'] = 'normal';

var layer_RoadLayerV3_1 = new L.geoJson.multiStyle(json_RoadLayerV3_1, {
    attribution: '',
    interactive: true,
    dataVar: 'json_RoadLayerV3_1',
    layerName: 'layer_RoadLayerV3_1',
    pane: 'pane_RoadLayerV3_1',
    onEachFeature: pop_RoadLayerV3_1,
    styles: [style_RoadLayerV3_1_0,style_RoadLayerV3_1_1,]
});

bounds_group.addLayer(layer_RoadLayerV3_1);
map.addLayer(layer_RoadLayerV3_1);
setBounds();

var roadLayer = null; // Variable to store the road layer
var roadlayerbounds = layer_RoadLayerV3_1.getBounds();
console.log(roadlayerbounds)
function toggleRoadLayer() {
    if (roadLayer) {
      map.removeLayer(roadLayer);
      map.fitBounds(roadlayerbounds);
      roadLayer = null; // Set to null when hidden
    } else {
      map.addLayer(layer_RoadLayerV3_1);
      roadLayer = layer_RoadLayerV3_1; // Store the reference for toggling
      roadlayerbounds = roadLayer.getBounds();
        map.fitBounds(roadlayerbounds);
    }
}
var roadLayerButtons = document.querySelector("#road-layer");
roadLayerButtons.addEventListener('click',toggleRoadLayer);
