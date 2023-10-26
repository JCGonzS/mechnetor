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

// var json_style = (function () {
//   var json_style = null;
//   $.ajax({
//       'async': false,
//       'global': false,
//       'url': style_json,
//       'dataType': "json",
//       'success': function (data) {
//         json_style = data;
//       }
//   });
//   return json_style;
// })();

var json_style = [
	{
		"selector": "node[role='protein_main']",
		"style": {
			"shape": "roundrectangle",
			"padding": "40px",
			"background-opacity": "0",
			"border-width": "2px",
			"border-color": "#2C2C2C",
			"border-opacity": "1",
			"label": "data(label)",
			"font-size": "30px",
			"font-weight": "bold",
			"text-halign": "center",
			"text-valign": "bottom",
			"text-margin-y": "-8px",
			"color": "#2C2C2C",
			"text-outline-width": "0",
			"text-outline-color": "#FFF",
			"text-background-opacity": "1",
			"text-background-color": "#FFF",
			"text-background-shape": "roundrectangle",
			"events": "yes",
			"text-events": "yes",
			"display": "data(display)"
		}
	},
	{
		"selector": "node[role='protein_main'].hl, node[role='protein_main']:selected",
		"style": {
			"border-color": "#17A589",
			"color": "#17A589",
			"overlay-padding": "5px",
			"overlay-color": "#64c2af",
			"overlay-opacity": "0.2",
			"z-index": "1000"
		}
	},
	{
		"selector": "node[role='protein_main'].hl2",
		"style": {
			"border-color": "#DD8C3D",
			"color": "#DD8C3D",
      		"overlay-padding": "5px",
			"overlay-color": "#DD8C3D",
			"overlay-opacity": "0.2",
		}
	},
	{
		"selector": "node[role='protein_seq']",
		"style": {
			"shape": "rectangle",
			"width": "data(length)",
			"height": "6px",
			"background-opacity": "1",
			"background-color": "#AAB7B8",
			"events": "no",
			"z-index-compare": "manual",
			"z-index": "1"
		}
	},
	{
		"selector": "node[role='domain']",
		"style": {
			"shape": "roundrectangle",
			"width": "data(length)",
			"height": "30px",
			"background-gradient-stop-colors": "data(colorgrad)",
			"background-gradient-stop-positions": "0% 25% 75% 100%",
			"background-fill": "linear-gradient",
			"background-opacity": "0.8",
			"border-width": "1px",
			"border-color": "data(color)",
			"label": "data(label)",
			"font-size": "10px",
			"font-weight": "bold",
			"min-zoomed-font-size": "15px",
			"text-valign": "bottom",
			"text-margin-y": "-5px",
			"text-wrap": "wrap",
			"text-max-width": "data(length)",
			"color": "#2C2C2C",
			"text-outline-color": "#FFF",
			"text-outline-width": "0.3px",
			"events": "yes",
			"text-events": "yes"
		}
	},
	{
		"selector": "node[role='domain'].hl, node[role='domain']:selected",
		"style": {
			"background-color": "#5D6D7E",
			"border-style": "dashed",
			"border-width": "2px",
			"label": function(ele){ return ele.data("label")+"\n"+ele.data("start")+"-"+ele.data("end"); },
      		"text-wrap": "wrap",
			"font-size": "15px",
			"text-outline-width": "2px",
			"min-zoomed-font-size": "5px"
		}
	},
	{
		"selector": "node[role='domain'].hl2",
		"style": {
			"border-style": "dashed",
			"border-width": "1px",
			"text-outline-width": "2px",
			"min-zoomed-font-size": "5px",
			"z-index": "999"
		}
	},
	{
		"selector": "node[role='elm'], node[role='3dlm']",
		"style": {
			"shape": "roundrectangle",
			"width": "data(length)",
			"height": "15px",
			"background-opacity": "0",
			"border-width": "1px",
			"border-color": "data(color)",
			"label": "data(label)",
			"font-size": "8px",
			"min-zoomed-font-size": "17px",
			"text-valign": "bottom",
			"color": "#2C2C2C",
			"events": "yes",
			"text-events": "yes",
			"display": "none"
		}
	},
	{
		"selector": "node[role='elm'].hl, node[role='elm'].hl2, node[role='elm']:selected, "+
                "node[role='3dlm'].hl, node[role='3dlm'].hl2, node[role='3dlm']:selected",
		"style": {
			"background-opacity": "0.5",
			"background-color": "data(color)",
      		"label": function(ele){ return ele.data("label")+"\n"+ele.data("start")+"-"+ele.data("end"); },
      		"text-wrap": "wrap",
			"font-size": "12px",
			"font-weight": "bold",
      		"text-outline-width": "1px",
     		"text-outline-color": "#FFF",
			"text-valign": "top",
			"text-margin-y": "2px",
			"z-index": "1000"
		}
	},
	{
		"selector": "node[role ^='uni']",
		"style": {
			"shape": "polygon",
			"shape-polygon-points": "1, 0.5, 0, 1, -1, 0.5, -1, -0.5, 0, -1, 1, -0.5",
			"height": "20px",
			"width": "data(length)",
			"background-opacity": "0.5",
			"border-color": "#000",
			"font-size": "10px",
			"min-zoomed-font-size": "15px",
			"text-background-opacity": "1",
			"text-background-color": "#FFF",
			"text-background-shape": "rectangle",
			"text-background-padding": "10px",
			"text-valign": "top",
			"text-margin-y": "-10px",
			"text-wrap": "wrap",
			"text-max-width": "300px",
			// "text-justification": "left",
			"color": "#000",
			"display": "none",
			"events": "yes",
			"text-events": "yes"
		}
	},
	{
		"selector": "node[role ^='uni'].hl, node[role ^='uni'].hl2",
		"style": {
			"height": "25px",
			"background-opacity": "1",
			"border-style": "dashed",
			"label": "data(label)",
			"overlay-padding": "3px",
			"overlay-opacity": "0.3",
			"z-index": "1000"
		}
	},

	{
		"selector": "node[role='uni_region']",
		"style": {
			"background-color": "#F1948A",
			"border-width": "1px"
		}
	},
	{
		"selector": "node[role='uni_region'].hl, node[role='uni_region'].hl2",
		"style": {
			"background-color": "#EC7063",
			"overlay-color": "#EC7063"
		}
	},
	{
		"selector": "node[role='uni_var']",
		"style": {
			"background-color": "#58D68D",
			"border-width": "0.3px"
		}
	},
	{
		"selector": "node[role='uni_var'].hl, node[role='uni_var'].hl2",
		"style": {
			"background-color": "#28B463",
			"overlay-color": "#28B463"
		}
	},
	{
		"selector": "node[role='uni_mtg']",
		"style": {
			"background-color": "#C39BD3",
			"border-width": "0.3px"
		}
	},
	{
		"selector": "node[role='uni_mtg'].hl, node[role='uni_mtg'].hl2",
		"style": {
			"background-color": "#9B59B6",
			"overlay-color": "#9B59B6"
		}
	},
	{
		"selector": "node[role='uni_binding']",
		"style": {
			"background-color": "#7fb3d5",
			"border-width": "0.3px"
		}
	},
	{
		"selector": "node[role='uni_binding'].hl, node[role='uni_binding'].hl2",
		"style": {
			"background-color": "#5499c7",
			"overlay-color": "#5499c7"
		}
	},
	{
		"selector": "node[role='uni_dnabind']",
		"style": {
			"background-color": "#d57fba",
			"border-width": "0.3px"
		}
	},
	{
		"selector": "node[role='uni_dnabind'].hl, node[role='uni_dnabind'].hl2",
		"style": {
			"background-color": "#d844a9",
			"overlay-color": "#d844a9"
		}
	},
	{
		"selector": "node[role='uni_metal']",
		"style": {
			"background-color": "#66666b",
			"border-width": "0.3px"
		}
	},
	{
		"selector": "node[role='uni_metal'].hl, node[role='uni_metal'].hl2",
		"style": {
			"background-color": "#4c4c50",
			"overlay-color": "#4c4c50"
		}
	},
	{
		"selector": "node[role='uni_transmem']",
		"style": {
			"shape": "hexagon",
			"background-color": "#789c4c",
			"border-width": "0.3px"
		}
	},
	{
		"selector": "node[role='uni_transmem'].hl, node[role='uni_transmem'].hl2",
		"style": {
			"background-color": "#648c33",
			"overlay-color": "#648c33"
		}
	},
	{
		"selector": "node[role='uni_disulfid']",
		"style": {
			"width": "2px",
			"shape": "round-tag",
			"background-color": "#f7dc6f",
			"border-width": "0.3px"
		}
	},
	{
		"selector": "node[role='uni_disulfid'].hl, node[role='uni_disulfid'].hl2",
		"style": {
			"background-color": "#f4d03f",
			"overlay-color": " #f4d03f "
		}
	},
	{
		"selector": "node[role='iprets']",
		"style": {
			"shape": "ellipse",
			"width": "data(length)",
			"background-opacity": "0",
			"border-width": "1px",
			"border-color": "#e44445",
			"label": "data(label)",
			"font-size": "12px",
			"color": "#2C2C2C",
			"text-outline-color": "#FFF",
			"text-outline-width": "0.3px",
			"min-zoomed-font-size": "24px",
			"text-valign": "bottom",
			"text-wrap": "wrap",
			"events": "yes",
			"text-events": "yes",
			"display": "none"
		}
	},
	{
		"selector": "node[role='iprets'].hl",
		"style": {
			"height": "40px",
			"background-opacity": "0.2",
			"background-color": "#e44445",
			"border-style": "dashed",
			"border-width": "3px",
			"font-size": "15px",
			"font-weight": "bold",
			"text-outline-width": "2px",
			"min-zoomed-font-size": "15px"
		}
	},
	{
		"selector": "node[role='iprets'].hl2",
		"style": {
			"border-style": "dashed",
			"border-width": "3px",
			"font-size": "15px",
			"font-weight": "bold",
			"text-outline-width": "2px",
			"min-zoomed-font-size": "15px"
		}
	},
	{
		"selector": "node[role ^='mod']",
		"style": {
			"border-width": "0.6px",
			"border-color": "#000",
			"color": "#2C2C2C",
			"font-size": "8px",
			"min-zoomed-font-size": "14px",
			"text-wrap": "wrap",
			"text-max-width": "75px",
			"text-margin-y": "-2px",
			"text-background-opacity": "1",
			"text-background-shape": "rectangle",
			"text-background-padding": "0.5px",
			"text-background-color": "#FFF",
			"text-border-style": "solid",
			"text-border-opacity": "1",
			"text-border-width": "0.3px",
			"text-border-color": "#2C2C2C",
			"display": "none",
			"events": "yes",
			"text-events": "yes"
		}
	},
	{
		"selector": "node[role ^='mod'].hl",
		"style": {
			"overlay-padding": "2px",
			"overlay-opacity": "0.3",
			"z-index": "1000"
		}
	},
	{
		"selector": "node[role='mod_cosmic']",
		"style": {
			"shape": "polygon",
			"shape-polygon-points": "1, -1, 0.1, -1, 0.1, 1, -0.1, 1, -0.1, -1, -1, -1",
			"height": "data(height)",
			"width": "3px",
			"border-width": "0.8px",
			"border-color": "#28446F",
			"background-color": "#2389AF"
		}
	},
	{
		"selector": "node[role='mod_cosmic'].hl",
		"style": {
			"shape": "polygon",
			"shape-polygon-points": "-1, -1, 1, -1, 0, 1",
			"height": "function(){ return data(height)+3; }",
			"width": "4px",
			"border-width": "0.5px",
			"border-color": "#000",
			"background-color": "#2389AF",
			"label": "data(aa_mut)",
			"overlay-color": "#2389AF"
		}
	},
	{
		"selector": "node[role='mod_input']",
		"style": {
			"shape": "polygon",
			"shape-polygon-points": "-1, -1, 1, -1, 0, 0.4",
			"height": "14px",
			"width": "3px",
			"background-color": "#FF0000"
		}
	},
	{
		"selector": "node[role='mod_input'].hl",
		"style": {
			"height": "16px",
			"width": "4px",
			"background-color": "#971031",
			"label": "data(label)",
			"overlay-color": "#FF0000"
		}
	},
	{
		"selector": "node[role='mod_phos']",
		"style": {
			"shape": "polygon",
			"shape-polygon-points": "0.2, -0.6, 0, -0.6, 0, 1, 0, -0.6, -0.2, -0.6, -0.2, -1, 0.2, -1",
			"height": "10px",
			"width": "10px",
			"background-color": "#FFFF00"
		}
	},
	{
		"selector": "node[role='mod_phos'].hl",
		"style": {
			"height": "12px",
			"background-color": "#FFB533",
			"label": "data(label)",
			"overlay-color": "#FFFF00"
		}
	},
	{
		"selector": "node[role='mod_acet']",
		"style": {
			"shape": "polygon",
			"shape-polygon-points": "0.2, -0.6, 0, -0.6, 0, 1, 0, -0.6, -0.2, -0.6, -0.2, -1, 0.2, -1",
			"height": "10px",
			"width": "10px",
			"background-color": "#008000"
		}
	},
	{
		"selector": "node[role='mod_acet'].hl",
		"style": {
			"height": "12px",
			"background-color": "#10972B",
			"label": "data(label)",
			"overlay-color": "#008000"
		}
	},
	{
		"selector": "edge[role='protein_sequence']",
		"style": {
			"width": "6px",
			"line-color": "#5F6A6A",
			"events": "no",
			"z-index-compare": "manual",
			"z-index": "0"
		}
	},
	{
		"selector": "edge[role='prot_prot_interaction']",
		"style": {
			"width": "mapData(evidence_n, 1, 50, 4, 30)",
			"line-color": "#5F6A6A",
			"line-style": "solid",
			"line-cap": "round",
			"opacity": "1",
			"arrow-scale": "0.8",
			"events": "yes",
			"z-index-compare": "manual",
			"z-index": "2500"
		}
	},
	{
		"selector": "edge[role='prot_prot_interaction']:selected, "+
                "edge[role='prot_prot_interaction'].hl, edge[role='prot_prot_interaction'].hl2",
		"style": {
			"width": "mapData(evidence_n, 1, 50, 8, 35)",
			"line-color": "#17A589",
			"source-arrow-color": "#17A589",
			"target-arrow-color": "#17A589",
      		"overlay-padding": "mapData(evidence_n, 1, 50, 8, 34)",
			"overlay-color": "#17A589",
			"overlay-opacity": "0.2"
		}
	},
	{
		"selector": "edge[role='user_interaction']",
		"style": {
			"width": "14px",
			"line-color": "#460DD4",
			"line-style": "solid",
			"opacity": "1",
			"arrow-scale": "0.8",
			"events": "yes"
		}
	},
	{
		"selector": "edge[role^='uni_']",
		"style": {
			"width": "4px",
			"curve-style": "bezier",
			"line-color": "#000",
			"line-style": "solid",
			"target-arrow-shape": "triangle",
			"target-arrow-color": "#000",
			"opacity": "0.7",
			"display": "none",
			"events": "yes"
		}
	},
	{
		"selector": "edge[role^='uni_']:selected, edge[role^='uni_'].hl",
		"style": {
			"line-style": "dashed",
			"opacity": "1",
			"overlay-opacity": "0.2",
			"label": "data(label)",
			"font-size": "10px",
			"min-zoomed-font-size": "15px",
			"text-background-opacity": "1",
			"text-background-color": "#FFF",
			"text-background-shape": "rectangle",
			"text-background-padding": "0",
			"text-valign": "top",
			"text-wrap": "wrap",
			"text-max-width": "100px",
			"color": "#000",
			"z-index": "1000"
		}
	},
	{
		"selector": "edge[role='uni_region_interaction']",
		"style": {
			"overlay-color": "#EC7063"
		}
	},
	{
		"selector": "edge[role='uni_var_interaction']",
		"style": {
			"overlay-color": "#28B463"
		}
	},
	{
		"selector": "edge[role='uni_mtg_interaction']",
		"style": {
			"overlay-color": "#9B59B6"
		}
	},
	{
		"selector": "edge[role='DOM_interaction']",
		"style": {
			"width": "mapData(pdb_n, 1, 50, 3, 20)",
			"curve-style": "bezier",
			"line-color": "#33a1c2",
			"line-style": "solid",
			"opacity": "0.7",
			"display": "none",
			"events": "yes"
		}
	},
	{
		"selector": "edge[role='DOM_interaction']:selected, "+
                	"edge[role='DOM_interaction'].hl, edge[role='DOM_interaction'].hl2",
		"style": {
			"line-style": "dashed",
			"opacity": "1",
			"overlay-padding": "mapData(pdb_n, 1, 50, 5, 22)",
			"overlay-color": "#33a1c2",
			"overlay-opacity": "0.2"
		}
	},
	{
		"selector": "edge[role='iDOM_interaction']",
		"style": {
			"width": "mapData(lo, 2, 10, 3, 20)",
			"curve-style": "bezier",
			"line-color": "#D4AC0D",
			"line-style": "solid",
			"opacity": "0.7",
			"display": "none",
			"events": "yes"
		}
	},
	{
		"selector": "edge[role='iDOM_interaction']:selected, "+
                	"edge[role='iDOM_interaction'].hl, edge[role='iDOM_interaction'].hl2",
		"style": {
			"line-style": "dashed",
			"opacity": "1",
			"overlay-padding": "mapData(lo, 2, 10, 6, 23)",
			"overlay-color": "#D4AC0D",
			"overlay-opacity": "0.2"
		}
	},
	{
		"selector": "edge[role='ELM_interaction']",
		"style": {
			"width": "4px",
			"curve-style": "bezier",
			"line-color": "#AF7AC5",
			"line-style": "solid",
			"opacity": "0.7",
			"display": "none",
			"events": "yes"
		}
	},
	{
		"selector": "edge[role='ELM_interaction'].hl",
		"style": {
			"line-style": "dashed",
			"line-color": "#884EA0",
			"opacity": "1",
			"overlay-padding": "6px",
			"overlay-color": "#884EA0",
			"overlay-opacity": "0.2"
		}
	},
	{
		"selector": "edge[role='iELM_interaction']",
		"style": {
			"width": "mapData(lo, 2, 10, 6, 23)",
			"line-style": "solid",
			"line-color": "#E74C3C",
			"curve-style": "bezier",
			"opacity": "0.7",
			"display": "none",
			"events": "yes"
		}
	},
	{
		"selector": "edge[role='iELM_interaction'].hl",
		"style": {
			"width": "mapData(lo, 2, 10, 3, 20)",
			"line-style": "dashed",
			"line-color": "#B03A2E",
			"opacity": "1"
		}
	},
	{
		"selector": "edge[role='DMI_interaction']",
		"style": {
			"width": "4px",
			"curve-style": "bezier",
			"line-color": "#e34a92",
			"line-style": "solid",
			"opacity": "0.7",
			"display": "none",
			"events": "yes"
		}
	},
	{
		"selector": "edge[role='DMI_interaction'].hl",
		"style": {
			"line-style": "dashed",
			"line-color": "#e34a92",
			"opacity": "1",
			"overlay-padding": "6px",
			"overlay-color": "#e34a92",
			"overlay-opacity": "0.2"
		}
	},
	{
		"selector": "edge[role='INT_interaction']",
		"style": {
			"width": "4px",
			"line-style": "solid",
			"line-color": "#e44445",
			"curve-style": "bezier",
			"opacity": "0.7",
			"events": "no",
			"display": "none",
			"events": "yes"
		}
	},
	{
		"selector": "edge[role='INT_interaction']:selected, edge[role='INT_interaction'].hl",
		"style": {
			"line-style": "dashed",
			"opacity": "1",
			"overlay-padding": "6px",
			"overlay-color": "#e44445",
			"overlay-opacity": "0.2",
			"z-index": "1000"
		}
	}
]


var cy;
document.addEventListener("DOMContentLoaded", function() {
  cy = cytoscape({
  	container: document.getElementById("cy"),
  	elements: json_elements,
    style: json_style,
  	layout: {
    	name: 'preset',
    	fit: true
  	},
    // motionBlur: true,
    minZoom: 0.05,
    maxZoom : 3.0,
    selectionType: 'single',
    wheelSensitivity: 0.3,
    // pixelRatio: 1,
    hideEdgesOnViewport: true, // These 2 options could be used only if network is big enough
    textureOnViewport: true
	});

  cy.fit( cy.$("node:visible") );
});
