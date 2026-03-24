// source: https://github.com/rstudio/leaflet/issues/496#issuecomment-650122985
window.LeafletWidget.methods.setStyle = function(category, layerId, style){
    var map = this;
    if (!layerId){
	return;
    } else if (!(typeof(layerId) === "object" && layerId.length)){ // in case a single layerid is given
	layerId = [layerId];
    }
    //convert columnstore to row store
    style = HTMLWidgets.dataframeToD3(style);
    layerId.forEach(function(d,i){
	var layer = map.layerManager.getLayer(category, d);
	if (layer){ // or should this raise an error?
	    layer.setStyle(style[i]);
	}
    });
};

// based on https://github.com/rstudio/leaflet/issues/496#issuecomment-651625559
window.LeafletWidget.methods.setLabel = function(category, layerId, label, options){
    var map = this;
    if (!layerId){
	return;
    } else if (!(typeof(layerId) === "object" && layerId.length)){ // in case a single layerid is given
	layerId = [layerId];
    }
    layerId.forEach(function(d,i){
	var layer = map.layerManager.getLayer(category, d);
	if (layer){ // or should this raise an error?
	    layer.unbindTooltip();
	    layer.bindTooltip(label[i], options)
	}
    });
};
