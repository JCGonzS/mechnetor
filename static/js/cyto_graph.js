"use strict";

var json_elements = (function () {
  var json_elements = null;
  $.ajax({
      'async': false,
      'global': false,
      'url': graph_json,
      'dataType': "json",
      'success': function (data) {
        json_elements = data;
      }
  });
  return json_elements;
})();

var json_style = (function () {
  var json_style = null;
  $.ajax({
      'async': false,
      'global': false,
      'url': style_json,
      'dataType': "json",
      'success': function (data) {
        json_style = data;
      }
  });
  return json_style;
})();

var cy;
document.addEventListener("DOMContentLoaded", function() {
  cy = cytoscape({
  	container: document.getElementById("cy"),
  	elements: json_elements,
  	layout: {
    	name: 'preset',
    	fit: true
  	},
    // motionBlur: true,
    style: json_style,
    selectionType: 'single',
    wheelSensitivity: 0.3
	});
});
