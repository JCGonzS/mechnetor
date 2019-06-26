"use strict";

// Show/Hide proteins without interactions
function showProt(label) {
	var node = cy.$("node[role='whole'][label='"+label+"'], node[role='user_seq'][label='"+label+"']");
	var display = node.style("display");
	if (display=="none"){
		node.style("display", "element");
		// cy.center(node);
	} else {
		node.style("display", "none");
	}
};

function removeThisEle(id) {
	var node = cy.$("node[id='"+id+"']");
	cy.remove(node);
};

function removeAllinProt(label, protein) {
	var node = cy.$("node[label='"+label+"'][protein='"+protein+"']");
	cy.remove(node);
};

function removeAll(label) {
	var node = cy.$("node[label='"+label+"']");
	cy.remove(node);
};

$(document).ready(function(){

	// Save original positions of all graph elements
	cy.nodes().forEach(function(n){
		var pos = n.position();
		n.data("orPos", {
			x: pos.x,
			y: pos.y
		});
	});

	// BUTTON: reset graph
	$("#reset").click(function(){
		// Remove all elements
		cy.elements().remove();
		// Add elements again
		cy.add( json_elements );
		// Set elements to original positions
		cy.nodes().forEach(function(n){
			var pos = n.data('orPos');
			n.position({
				x: pos.x,
				y: pos.y
		 	});
		});
		// Resets the zoom and pan.
		cy.reset();
		// Fit view to visible elements
		cy.fit( cy.$("node:visible") );
		// Check/Uncheck to reflect default status
		$("#toggle_box").prop("checked", true);
		$("#toggle_doms").prop("checked", true);
		$("#toggle_elms").prop("checked", false);
		$("#toggle_iprets").prop("checked", false);
		$("#toggle_pp_int").prop("checked", true);
		$("#toggle_dom_int").prop("checked", false);
		$("#toggle_idom_int").prop("checked", false);
		$("#toggle_elmdom_int").prop("checked", false);
		$("#toggle_prets_int").prop("checked", false);
		$("#toggle_phos").prop("checked", false);
		$("#toggle_acet").prop("checked", false);
		$("#toggle_input_mut").prop("checked", false);
		$("#toggle_cosmic_mut").prop("checked", false);
	});

	// BUTTON: Center the graph on the page
  $("#center").click(function(){
		// cy.center();
		cy.fit( cy.$("node:visible") );
  });

  // BUTTON: Lock/Unlock element positions
  $("#lock").click(function(){
		if ( cy.autolock()==false ) {
			cy.autolock(true);
      $("#lock i").toggleClass("fa-unlock fa-lock");
    } else if ( cy.autolock()==true ) {
      cy.autolock(false);
      $("#lock i").toggleClass("fa-lock fa-unlock");
    }
  });

	// BUTTON: Zoom-In / Zoom-Out
	$("#zoom_in").click(function(){
	 var z = cy.zoom() + 0.2
	 cy.zoom( z )
	});
	$("#zoom_out").click(function(){
	 var z = cy.zoom() - 0.2
	 cy.zoom( z )
	});

	/// TOGGLE TOOLS:
  // BUTTON: Toggle box surrounding the whole protein
  $("#toggle_box").click(function(){
		var eles = cy.$("node[role='whole'], node[role='user_seq']");
    var eles2 = cy.$("node[role='whole']:selected, node[role='user_seq']:selected");
    var checked = document.getElementById("toggle_box").checked;
    if (checked) {
      eles.style("border-opacity", "1");
    } else {
			eles.style("border-opacity", "0");
			eles2.style("border-opacity", "1");
    }
  });

  // BUTTON: Toggle ALL protein elements
	$("#toggle_all_ele").click(function(){
		var dom_nodes = cy.$("node[role='domain']");
		var dom_edges = dom_nodes.connectedEdges();
		var elm_nodes = cy.$("node[role='elm']");
		var elm_edges = elm_nodes.connectedEdges();
		var ip_nodes = cy.$("node[role='iprets']");
		var ip_edges = ip_nodes.connectedEdges();
		var checked = document.getElementById("toggle_all_ele").checked;
		if (checked) {
			dom_nodes.style("display", "element");
			elm_nodes.style("display", "element");
			ip_nodes.style("display", "element");
			$("#toggle_doms").prop("checked", true);
			$("#toggle_elms").prop("checked", true);
			$("#toggle_iprets").prop("checked", true);
		} else {
			dom_nodes.style("display", "none");
			dom_edges.style("display", "none");
			elm_nodes.style("display", "none");
			elm_edges.style("display", "none");
			ip_nodes.style("display", "none");
			ip_edges.style("display", "none");
			$("#toggle_doms").prop("checked", false);
			$("#toggle_elms").prop("checked", false);
			$("#toggle_iprets").prop("checked", false);
			$("#toggle_dom_int").prop("checked", false);
			$("#toggle_idom_int").prop("checked", false);
			$("#toggle_elmdom_int").prop("checked", false);
			$("#toggle_prets_int").prop("checked", false);
		}
	});

  // BUTTON: Toggle domains
  $("#toggle_doms").click(function(){
    var nodes = cy.$("node[role='domain']");
		var edges = nodes.connectedEdges();
    var checked = document.getElementById("toggle_doms").checked;
    if (checked) {
			nodes.style("display", "element");
    } else {
			nodes.style("display", "none");
			edges.style("display", "none");
			$("#toggle_dom_int").prop("checked", false);
			$("#toggle_idom_int").prop("checked", false);
			$("#toggle_elmdom_int").prop("checked", false);
    }
  });

  // BUTTON: Toggle elms
  $("#toggle_elms").click(function(){
    var nodes = cy.$("node[role='elm']");
		var edges = nodes.connectedEdges();
    var checked = document.getElementById("toggle_elms").checked;
    if (checked) {
      nodes.style("display", "element");
    } else {
			nodes.style("display", "none");
			edges.style("display", "none");
			$("#toggle_elmdom_int").prop("checked", false);
    }
  });

  // BUTTON: Toggle InterPreTs regions
  $("#toggle_iprets").click(function(){
    var nodes = cy.$("node[role='iprets']");
		var edges = nodes.connectedEdges();
    var checked = document.getElementById("toggle_iprets").checked;
    if (checked) {
      nodes.style("display", "element");
    } else {
			nodes.style("display", "none");
			edges.style("display", "none");
			$("#toggle_prets_int").prop("checked", false);
    }
  });

  // BUTTON: Toggle ALL interactions
	$("#toggle_all_int").click(function(){
		var prot_edges = cy.$("edge[role='prot_prot_interaction']");
		var user_edges = cy.$("edge[role='user_interaction']");
		var dom_nodes = cy.$("node[role='domain']");
		var dom_edges = dom_nodes.connectedEdges();
		var elm_nodes = cy.$("node[role='elm']");
		var elm_edges = elm_nodes.connectedEdges();
		var ip_nodes = cy.$("node[role='iprets']");
		var ip_edges = ip_nodes.connectedEdges();
		var checked = document.getElementById("toggle_all_int").checked;
		if (checked) {
			prot_edges.style("display", "element");
			user_edges.style("display", "element");
			dom_nodes.style("display", "element");
			dom_edges.style("display", "element");
			elm_nodes.style("display", "element");
			elm_edges.style("display", "element");
			ip_nodes.style("display", "element");
			ip_edges.style("display", "element");
			$("#toggle_doms").prop("checked", true);
			$("#toggle_elms").prop("checked", true);
			$("#toggle_iprets").prop("checked", true);
			$("#toggle_pp_int").prop("checked", true);
			$("#toggle_user_int").prop("checked", true);
			$("#toggle_dom_int").prop("checked", true);
			$("#toggle_idom_int").prop("checked", true);
			$("#toggle_elmdom_int").prop("checked", true);
			$("#toggle_prets_int").prop("checked", true);
		} else {
			prot_edges.style("display", "none");
			user_edges.style("display", "none");
			dom_edges.style("display", "none");
			elm_edges.style("display", "none");
			ip_edges.style("display", "none");
			$("#toggle_pp_int").prop("checked", false);
			$("#toggle_user_int").prop("checked", false);
			$("#toggle_dom_int").prop("checked", false);
			$("#toggle_idom_int").prop("checked", false);
			$("#toggle_elmdom_int").prop("checked", false);
			$("#toggle_prets_int").prop("checked", false);
		}
	});

	// BUTTON: Toggle user-input interactions
  $("#toggle_user_int").click(function(){
		var edges = cy.$("edge[role='user_interaction']");
    var checked = document.getElementById("toggle_user_int").checked;
		if (checked) {
      edges.style("display", "element");
    } else {
      edges.style("display", "none");
    }
  });

  // BUTTON: Toggle protein-protein interactions
  $("#toggle_pp_int").click(function(){
    var edges = cy.$("edge[role='prot_prot_interaction']");
    var checked = document.getElementById("toggle_pp_int").checked;
		if (checked) {
      edges.style("display", "element");
    } else {
      edges.style("display", "none");
    }
  });

  // BUTTON: Toggle domain-domain interactions (3did)
  $("#toggle_dom_int").click(function(){
    var edges = cy.$("edge[role='DOM_interaction']");
		var nodes = edges.connectedNodes();
    var checked = document.getElementById("toggle_dom_int").checked;
		if (checked) {
      edges.style("display", "element");
			nodes.style("display", "element");
			$("#toggle_doms").prop("checked", true);
    } else {
      edges.style("display", "none");
    }
  });

  // BUTTON: Toggle domain-domain interactions (statistical prediction)
  $("#toggle_idom_int").click(function(){
    var edges = cy.$("edge[role='iDOM_interaction']");
		var nodes = edges.connectedNodes();
    var checked = document.getElementById("toggle_idom_int").checked;
		if (checked) {
      edges.style("display", "element");
			nodes.style("display", "element");
			$("#toggle_doms").prop("checked", true);
    } else {
      edges.style("display", "none");
    }
  });

  // BUTTON: Toggle elm-domain interactions
  $("#toggle_elmdom_int").click(function(){
    var edges = cy.$("edge[role='ELM_interaction']");
		var nodes = edges.connectedNodes();
    var checked = document.getElementById("toggle_elmdom_int").checked;
		if (checked) {
      edges.style("display", "element");
			nodes.style("display", "element");
			$("#toggle_elms").prop("checked", true);
    } else {
      edges.style("display", "none");
    }
  });

  // BUTTON: Toggle InterPreTS interactions
  $("#toggle_prets_int").click(function(){
    var edges = cy.$("edge[role='INT_interaction']");
		var nodes = edges.connectedNodes();
    var checked = document.getElementById("toggle_prets_int").checked;
		if (checked) {
      edges.style("display", "element");
			nodes.style("display", "element");
			$("#toggle_iprets").prop("checked", true);
    } else {
      edges.style("display", "none");
    }
  });

	// BUTTON: Toggle ALL protein modifications
	$("#toggle_all_mod").click(function(){
		var phos = cy.$("node[role='phosphorylation']");
		var ace = cy.$("node[role='acetylation']");
		var input_mut = cy.$("node[role='input_mut']");
		var cosmic_mut = cy.$("node[role='cosmic_mut']");
		var checked = document.getElementById("toggle_all_mod").checked;
		if (checked) {
			phos.style("display", "element");
			ace.style("display", "element");
			input_mut.style("display", "element");
			cosmic_mut.style("display", "element");
			$("#toggle_phos").prop("checked", true);
			$("#toggle_acet").prop("checked", true);
			$("#toggle_input_mut").prop("checked", true);
			$("#toggle_cosmic_mut").prop("checked", true);
		} else {
			phos.style("display", "none");
			ace.style("display", "none");
			input_mut.style("display", "none");
			cosmic_mut.style("display", "none");
			$("#toggle_phos").prop("checked", false);
			$("#toggle_acet").prop("checked", false);
			$("#toggle_input_mut").prop("checked", false);
			$("#toggle_cosmic_mut").prop("checked", true);
		}
	});

	// BUTTON: Toggle mutations
	$("#toggle_input_mut").click(function(){
		var eles = cy.$("node[role='input_mut']");
		var checked = document.getElementById("toggle_input_mut").checked;
		if (checked) {
			eles.style("display", "element");
		} else {
			eles.style("display", "none");
		}
	});

	$("#toggle_cosmic_mut").click(function(){
		var eles = cy.$("node[role='cosmic_mut']");
		var checked = document.getElementById("toggle_cosmic_mut").checked;
		if (checked) {
			eles.style("display", "element");
		} else {
			eles.style("display", "none");
		}
	});

	// BUTTON: Toggle phosphorylations
  $("#toggle_phos").click(function(){
    var eles = cy.$("node[role='phosphorylation']");
    var checked = document.getElementById("toggle_phos").checked;
		if (checked) {
      eles.style("display", "element");
    } else {
      eles.style("display", "none");
    }
  });

  // BUTTON: Toggle acetylations
  $("#toggle_acet").click(function(){
    var eles = cy.$("node[role='acetylation']");
    var checked = document.getElementById("toggle_acet").checked;
		if (checked) {
      eles.style("display", "element");
    } else {
      eles.style("display", "none");
    }
  });

	/// DOWNLOAD BUTTONS
	// BUTTON: Download graph as JSON
	$("#dl_json").click(function(){
		var jsonBlob = new Blob([ JSON.stringify( cy.json() ) ], { type: "application/javascript;charset=utf-8" });
		saveAs( jsonBlob, "graph.json" );
	});
  // BUTTON: Get snapshot as PNG
  $("#dl_png").click(function(){
		var image = cy.png()
		var iframe = "<iframe src='"+image+`' frameborder='0'
			style='border:0; top:0px; left:0px; bottom:0px; right:0px; width:100%; height:100%;'
			allowfullscreen></iframe>`;
  	var win = window.open();
	  win.document.write(iframe);
  });
  // BUTTON: Get snapshot as JPG
  $("#dl_jpg").click(function(){
		var image = cy.jpg()
		var iframe = "<iframe src='"+image+`' frameborder='0'
			style='border:0; top:0px; left:0px; bottom:0px; right:0px; width:100%; height:100%;'
			allowfullscreen></iframe>`;
  	var win = window.open();
	  win.document.write(iframe);
  });
  $("#dl_svg").click(function(filename) {
		var svgContent = cy.svg({scale: 1, full: true});
		var blob = new Blob([svgContent], {type:"image/svg+xml;charset=utf-8"});
		saveAs(blob, "graph.svg");
	};

	// Create double-Tap event
	var tappedBefore;
	var tappedTimeout;
	cy.on("tap", function(event) {
	  var tappedNow = event.target;
	  if (tappedTimeout && tappedBefore) {
	    clearTimeout(tappedTimeout);
	  }
	  if(tappedBefore === tappedNow) {
	    tappedNow.trigger("doubleTap");
	    tappedBefore = null;
	  } else {
	    tappedTimeout = setTimeout(function(){ tappedBefore = null; }, 300);
	    tappedBefore = tappedNow;
	  }
	});

	// FEATURE: Double-Tap on protein node to zoom-in
	cy.on("doubleTap", "node[role='whole'], node[role='user_seq']", function(event) {
	  var node = event.target;
		cy.animate({
			fit: { eles: node,
						 padding: 100}
		})
	});

	// FEATURE: Some nodes & edges' sizes will vary with zoom levels
	cy.on("render zoom", function(event) {
		var node = cy.$("node[role='whole'], node[role='user_seq']");
		var edges = node.connectedEdges();
		var dim = 16/cy.zoom();
		var maxDim = Math.max(dim,30);
		node.style({"font-size": maxDim,
								"text-outline-width": maxDim/10,
								"border-width": maxDim/18,
								"text-margin-y": maxDim/-4
							});
	// if (cy.zoom() <= 0.5){
	// 	node.style("visibility", "hidden")
	// }
	});

	// MOUSEOVER/MOUSEOUT changes
	cy.autoungrabify(true);
  cy.on("tapdragover tapdragout","node[role='whole'], node[role='user_seq']", function(event) {
    var node = event.target;
		node.toggleClass("highlight");
		node.connectedEdges().toggleClass("highlight");
		node.neighborhood().toggleClass("highlight2");
		if (event.type=="tapdragover"){
			cy.autoungrabify(false);
		} else if (event.type=="tapdragout"){
			cy.autoungrabify(true);
		}
  });

	cy.on("tapdragover tapdragout",
				"node[role='domain'], node[role='elm'], node[role='iprets'], "+
				"node[role='phosphorylation'], node[role='acetylation'], "+
				"node[role='input_mut'], node[role='cosmic_mut']", function(event) {
		var node = event.target;
		node.toggleClass("highlight");
		node.connectedEdges().toggleClass("highlight");
		node.neighborhood().toggleClass("highlight2");
		// if (event.type=="tapdragover"){
		// 	node.ungrabify();
		// } else if (event.type=="tapdragout") {
		// 	node.grabify();
		// }
	});

	cy.on("mouseover mouseout", "edge[role='prot_prot_interaction']", function(event) {
		var edge = event.target;
		edge.toggleClass("highlight");
		edge.connectedNodes().toggleClass("highlight2");
	});

	cy.on("mouseover mouseout", "edge[role='DOM_interaction']", function(event) {
		var edge = event.target;
		edge.toggleClass("highlight");
		edge.connectedNodes().toggleClass("highlight2");
		edge.connectedNodes().parent().toggleClass("highlight2");
	});

	cy.on("mouseover mouseout", "edge[role='iDOM_interaction']", function(event) {
		var edge = event.target;
		edge.toggleClass("highlight");
		edge.connectedNodes().toggleClass("highlight2");
		edge.connectedNodes().parent().toggleClass("highlight2");
	});

	cy.on("mouseover mouseout", "edge[role='ELM_interaction']", function(event) {
		var edge = event.target;
		edge.toggleClass("highlight");
		edge.connectedNodes().toggleClass("highlight2");
		edge.connectedNodes().parent().toggleClass("highlight2");
	});

	cy.on("mouseover mouseout", "edge[role='INT_interaction']", function(event) {
		var edge = event.target;
		edge.toggleClass("highlight");
		edge.connectedNodes().toggleClass("highlight2");
		edge.connectedNodes().parent().toggleClass("highlight2");
	});

	// QTIPs on node/edge click
	cy.nodes().on("tap", function( e ){
		var node = e.target;
		var role = node.data("role");
		if( role == "whole" ){
			var gene = node.data("label");
			var des = node.data("des");
			var acc = node.data("protein");
			var length = node.data("length");
			node.qtip({
				content:
					"<span class='tip'>"+
					"<span class='tipProt'>Gene</span> | <b>"+gene+"</b><br>" +
					"<span class='tipProt'>Protein</span> | <b>"+des+"</b><br>" +
					"<span class='tipProt'>Accession</span> | " +
					"<a href='https://www.uniprot.org/uniprot/"+acc+"'>" +
					acc+" <i class='fas fa-external-link-alt fa-xs'></i>" +
					"</a><br>" +
					"<span class='tipProt'>Length</span> | "+
					"<a href='https://www.uniprot.org/uniprot/"+acc+".fasta'>" +
					length+" AA <i class='fas fa-external-link-alt fa-xs'></i>"+
					"</span>",
				position: {
					my: "top center",
					at: "bottom center"
				},
				style: {
					classes: "qtip-bootstrap",
					tip: {
						width: 20,
						height: 10
					}
				},
				show: { event: "directtap" }
			});

			this.trigger("directtap");
		}
	});

	cy.nodes().on("tap", function( e ){
		var node = e.target;
		var role = node.data("role");
		if( role == "user_seq" ){
			var label = node.data("label");
			var length = node.data("length");
			var link = node.data("link");
			var blast = node.data("blast");
			var blast_hits = "<span class='tipProt'>Blast hits:</span><br>\n";
			blast_hits += "<div style='margin: 0 auto;'>";
			blast.forEach(function(hit) {
				var coor = hit.split("|")[0]
				var acc = hit.split("|")[1]
				var gene = hit.split("|")[2]
				var info = hit.split("|")[3]
				blast_hits+="\t<span style='float:left;width:80px;'><b>"+coor+"</b></span>"
				blast_hits+=" | <a href='https://www.uniprot.org/uniprot/"+acc+"'>"+acc+"</a>"
				blast_hits+=" <b>"+gene+"</b> "+info+"<br>\n"
			});
			blast_hits += "</div>"
			node.qtip({
				content:
					"<span class='tip'>User-input sequence<br>\n"+
					"\t<span class='tipProt'>Label</span> | <b>"+label+"</b><br>\n" +
					"\t<span class='tipProt'>Length</span> | "+
					"<a href='"+link+"' target='_blank'_>" +
					length+" AA <i class='fas fa-external-link-alt fa-xs'></i></a><br>"+
					blast_hits+
					"</span>",
				position: {
					my: "top center",
					at: "bottom center"
				},
				style: {
					classes: "qtip-bootstrap",
					tip: {
						width: 20,
						height: 10
					}
				},
				show: { event: "directtap" }
			});

			this.trigger("directtap");
		}
	});

	cy.on("click", "node[role='domain']", function(event) {
		var node = event.target;
		var name = node.data("label");
		var acc = node.data("acc");
		var des = node.data("des");
		var start = node.data("start");
		var end = node.data("end");
		var prot = node.data("protein");
		node.qtip({
			content:
				"<span class='tip' style='color: #074987;'>" +
				"<b>Source: <a href='https://pfam.xfam.org'>Pfam</a></b>" +
				"</span><br>" +
				"<span class='tip'>" +
				"<span class='tipPfam'>Family</span> | " +
				"<b><i>"+name+"</i></b> (<a href='https://pfam.xfam.org/family/"+acc+"'>"+acc+" <i class='fas fa-external-link-alt fa-xs'></i></a>)<br>" +
				"<span class='tipPfam'>Description</span> | <b>"+des+"</b><br>" +
				"<span class='tipPfam'>Start - End</span> | " +
				"<b>"+start+"</b> - <b>"+end+"</b>" +
				" (<a href='https://pfam.xfam.org/protein/"+prot+"'>"+prot+" <i class='fas fa-external-link-alt fa-xs'></i></a>)"+
				"</span>",
			position: {
				my: "top center",
				at: "bottom center"
			},
			style: {
				classes: "qtip-bootstrap",
				tip: {
					width: 20,
					height: 10
				}
			}
		});
	});

	cy.on("click","node[role='elm']", function(event) {
		var node = event.target;
		var name = node.data("label");
		var id = node.data("id");
		var acc = node.data("acc");
		var i = 0;
		var des = "";
		node.data("des").split("").forEach(function(char) {
			i += 1;
			des += char;
			if (i == 60) {
				i = 0;
				des += "<br>";
			}
		});
		var regex = node.data("regex");
		var start = node.data("start");
		var end = node.data("end");
		var seq = node.data("seq");
		var prot = node.data("protein");
		node.qtip({
			content:
				"<span class='tip' style='color: #7f7c7b;'>\n" +
					"<b>Source: <a href='https://elm.eu.org'>ELM</a></b>\n" +
				"</span><br>\n" +
				"<span class='tip'>\n"+
					"<span class='tipELM'>Identifier</span> | "+
					"<a href='http://elm.eu.org/elms/"+name+"'>"+name+" <i class='fas fa-external-link-alt fa-xs'></i></a><br>\n"+
					"<span class='tipELM'>Accession</span> | "+ acc + "<br>\n"+
					"<span class='tipELM'>Description</span> | <b>"+des+"</b><br>\n"+
					"<span class='tipELM'>Start - End</span> | " +
					"<b>"+start+"</b> - <b>"+end+"</b>" +
					" (<a href='http://elm.eu.org/instances/"+name+"/"+prot+"/'>"+prot+" <i class='fas fa-external-link-alt fa-xs'></i></a>)<br>\n"+
					"<span class=tipELM>Subsequence</span> | <i>"+seq+"</i><br>\n" +
				"</span>\n"+
				"<div class='row'>\n"+
					"<div class='col-sm-4 text-center'>\n"+
						"<button class='txt-btn' onclick=removeThisEle(\""+id+"\") >\n"+
							"Remove\n"+
						"</button>\n"+
					"</div>\n"+
					"<div class='col-sm-4 text-center'>\n"+
						"<button class='txt-btn' onclick=removeAllinProt(\""+name+"\",\""+prot+"\") >\n"+
							"Remove all in this protein\n"+
						"</button>\n"+
					"</div>\n"+
					"<div class='col-sm-4 text-center'>\n"+
						"<button class='txt-btn' onclick=removeAll(\""+name+"\") >\n"+
							"Remove all in network\n"+
						"</button>\n"+
					"</div>\n"+
				"</div>",
			position: {
				my: "top center",
				at: "bottom center"
			},
			style: {
				classes: "qtip-bootstrap",
				tip: {
					width: 20,
					height: 10
				}
			}
		});
	});

	cy.on("click","node[role='iprets']", function(event) {
		var node = event.target;
		var pdb = node.data("pdb").split("|")[1];
		var chain = node.data("pdb").split("|")[2];
		var prot = node.data("protein");
		var start = node.data("start");
		var end = node.data("end");
		var pdb_start = node.data("pdb_start");
		var pdb_end = node.data("pdb_end");
		var ev = node.data("eval");
		var pcid = node.data("pcid");

		node.qtip({
			content:
				"<span class='tip' style='color: #7b241c;'>" +
				"<b>Predicted with <a href='http://www.russelllab.org/cgi-bin/tools/interprets.pl/interprets.pl'>InterPreTS</a></b>" +
				"</span><br>" +
				"<span class='tip'>"+
					"<span class='tipInP'>PDB Template ID</span> | " +
					"<a href='https://www.rcsb.org/structure/"+pdb+"'><b>"+pdb+"</b> <i class='fas fa-external-link-alt fa-xs'></i></a> chain "+chain+"<br>"+
					"<span class='tipInP'>Alignment</span> |  "+
					"Protein/<b>"+start+"-"+end+"</b>; Template/"+pdb_start+"-"+pdb_end+"<br>"+
					"<span class='tipInP'>Alignment score</span> | BLAST e-val=<i>"+ev+"</i>, <i>"+pcid+"%</i> id<br>"+
				"</span>",
			position: {
				my: "top center",
				at: "bottom center"
			},
			style: {
				classes: "qtip-bootstrap",
				tip: {
					width: 20,
					height: 10
				}
			}
		});
	});

	cy.on("click","node[role='cosmic_mut']", function(event) {
		var node = event.target;
		var cos_id = node.data("cos_id").split(";");
		var mut_aas = node.data("aa_mut").split(";");
		var mut_cds = node.data("cds").split(";");
		var counts = node.data("count").split(";");
		var html = "<table style='width:100%'>\n";
		html += "\t<tr>\n";
		html += "\t\t<th>COSMIC ID</th>\n";
		html += "\t\t<th>AA mutation</th>\n";
		html += "\t\t<th>CDS mutation</th>\n";
		html += "\t\t<th>Samples</th>\n";
		html += "\t</tr>\n";
		for (var i = 0; i < mut_aas.length; i++){
			html += "\t<tr>\n";
			html += "\t\t<td><a href='https://cancer.sanger.ac.uk/cosmic/mutation/overview?id="+cos_id[i].split("COSM")[1]+"'><span>"+cos_id[i]+"</span></a></td>\n";
			html += "\t\t<td>"+mut_aas[i]+"</td>\n";
			html += "\t\t<td>"+mut_cds[i]+"</td>\n";
			html += "\t\t<td>"+counts[i]+"</td>\n";
			html += "\t</tr>\n";
		};
		html += "</table>\n";
		// console.log(html);
		var gene = node.parent().data("label");
		var gene_link = "https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln="+gene

		node.qtip({
			content:
				"<span class='tip' style='color: #28446f;'>" +
					"<b><a href='https://cancer.sanger.ac.uk/cosmic'>COSMIC <i class='fas fa-external-link-alt fa-xs'></i></a></b><br>" +
					"<a href="+gene_link+">View gene "+gene+" <i class='fas fa-external-link-alt fa-xs'></i></a>"+
				"</span><br>\n"+html,
			position: {
				my: "top center",
				at: "bottom center"
			},
			style: {
				classes: "qtip-bootstrap",
				tip: {
					width: 20,
					height: 10
				}
			}
		});

	});

	cy.on("click","edge[role='prot_prot_interaction']", function(event) {
			var edge = event.target;
			var nodes = edge.connectedNodes();
			var genes = []
			var bioids = []
			nodes.forEach(function(node) {
				genes.push(node.data("label"));
				bioids.push(node.data("biogrid_id"));
			});
			var ds = edge.data("ds");
			var low = edge.data("low");
			var low_len = low.length;
			var high = edge.data("high");
			var high_len = high.length;
			var low_links = [];
			var low_link = ""
			low.slice(0,3).forEach(function(link) {
				low_links.push("<a href='https://thebiogrid.org/interaction/"+link+"'>"+link+" <i class='fas fa-external-link-alt fa-xs'></i></a>");
			});
			if (low_links.length > 0){
				low_link = " eg. ".concat(low_links.join(", "));
				if (low_len > 3){
					low_link = low_link.concat(" ...");
				}
			}
			var high_links = [];
			var high_link = ""
			high.slice(0,3).forEach(function(link) {
				high_links.push("<a href='https://thebiogrid.org/interaction/"+link+"'>"+link+" <i class='fas fa-external-link-alt fa-xs'></i></a>");
			});
			if (high_links.length > 0){
				high_link = " eg. ".concat(high_links.join(", "));
				if (high_len > 3){
					high_link = high_link.concat(" ...");
				}
			}
			edge.qtip({
				content:
					"<span class='tip' style='color: #33a1c2;'>" +
					 	"Source: <a href='https://thebiogrid.org/'><b>BioGRID</b></a></b>" +
					"</span><br>" +
					"<span class='tip'>"+
						"<span class='tipBioG'>Interacting Proteins</span> | "+
						"<a href='https://thebiogrid.org/"+bioids[0]+"'><b>"+genes[0]+"</b> <i class='fas fa-external-link-alt fa-xs'></i></a> - "+
						"<a href='https://thebiogrid.org/"+bioids[1]+"'><b>"+genes[1]+"</b> <i class='fas fa-external-link-alt fa-xs'></i></a>"+
					  "<br>"+
						"<span class='tipBioG'>LT evidence</span> | "+
						"<b>"+low_len+"</b><span style='font-size:12px;'>"+low_link+"</span><br>"+
						"<span class='tipBioG'>HT evidence</span> | "+
						"<b>"+high_len+"</b><span style='font-size:12px;'>"+high_link+"</span>"+
					"</span><br>",
				position: {
					my: "top center",
					at: "bottom center"
				},
				style: {
					classes: "qtip-bootstrap",
					tip: {
						width: 20,
						height: 10
			    }
			  }
			});
	});

	cy.on("click","edge[role='DOM_interaction']", function(event) {
			var edge = event.target;
			var nodes = edge.connectedNodes();
			var doms = [];
			var pfams = [];
			var prots = [];
			nodes.forEach(function(node) {
				pfams.push(node.data("label"));
				doms.push("<span style='color: "+node.data("color")+";'><b>"+node.data("label")+"</b></span>");
				prots.push(node.parent().data("label"));
			});
			var ds = edge.data("ds");

			edge.qtip({
				content:
					"<span class='tip' style='color: #16a085;'>" +
						"Source: <a href='https://3did.irbbarcelona.org/'><b>3did</b></a></b>" +
					"</span><br>" +
					"<span class='tip'>"+
						"<span class='tip3did'>Interacting Domains</span> | "+
				 		"<a href='https://3did.irbbarcelona.org/dispatch.php?type=interaction&type1=domain&type2=domain&value1="+pfams[0]+"&value2="+pfams[1]+"'>"+
							doms.join(" - ")+" <i class='fas fa-external-link-alt fa-xs'></i>"+
						"</a>"+
					"</span><br>"+
					"<span class='tip'>"+
						"<span class='tip3did'>Interacting Proteins</span> | "+prots.join(" - ")+
					"</span>",
				position: {
					my: "top center",
					at: "bottom center"
				},
				style: {
					classes: "qtip-bootstrap",
					tip: {
						width: 20,
						height: 10
			    }
			  }
			});
	});

	cy.on("click","edge[role='iDOM_interaction']", function(event) {
			var edge = event.target;
			var nodes = edge.connectedNodes();
			var doms = [];
			var pfams = [];
			var prots = [];
			nodes.forEach(function(node) {
				pfams.push(node.data("label"));
				doms.push("<span style='color: "+node.data("color")+";'><b>"+node.data("label")+"</b></span>");
				prots.push(node.parent().data("label"));
			});
			var lo = edge.data("lo");
			var ds = edge.data("ds");

			edge.qtip({
				content:
					"<span class='tip' style='color: #d4ac0d;'>" +
						"<b>Inferred Domain-Domain ("+ds+")</b>" +
					"</span><br>" +
					"<span class='tip'>"+
						"<span class='tipIdom'>Interacting Domains</span> | "+doms.join(" - ")+
					"</span><br>"+
					"<span class='tip'>"+
						"<span class='tipIdom'>Interacting Proteins</span> | "+prots.join(" - ")+
					"</span><br>" +
					"<span class='tip'>"+
						"<span class='tipIdom'>Association Score</span> | "+lo+
					"</span>",

				position: {
					my: "top center",
					at: "bottom center"
				},
				style: {
					classes: "qtip-bootstrap",
					tip: {
						width: 20,
						height: 10
			    }
			  }
			});
	});

	cy.on("click","edge[role='ELM_interaction']", function(event) {
			var edge = event.target;
			var nodes = edge.connectedNodes();
			var doms = [];
			var pfams = [];
			var prots = [];
			nodes.forEach(function(node) {
				pfams.push(node.data("label"));
				doms.push("<span style='color: "+node.data("color")+";'><b>"+node.data("label")+"</b></span>");
				prots.push(node.parent().data("label"));
			});
			var ds = edge.data("ds");

			edge.qtip({
				content:
					"<span class='tip' style='color: #b95db9;'>" +
						"Source: <a href='http://elm.eu.org/downloads.html#interactions'><b>ELM</b> (interactions)</a></b>" +
					"</span><br>" +
					"<span class='tip'>"+
						"<span class='tipELMint'>Interacting Elements</span> | "+
							doms.join(" - ")+
						"</a><br>"+
						"<span class='tipELMint'>Interacting Proteins</span> | "+prots.join(" - ")+
					"</span>",
				position: {
           my: "top center",
           at: "bottom center"
       	},
				style: {
					classes: "qtip-bootstrap",
					tip: {
						width: 20,
						height: 10
			    }
			  }
			});
	});

	cy.on("click","edge[role='INT_interaction']", function(event) {
		var edge = event.target;
		var nodes = edge.connectedNodes();
		var pdb = edge.data("pdb").split("|")[1];
		var z = edge.data("z-score");
		var pval = edge.data("p-value");
		var chains = [];
		var prots = [];
		nodes.forEach(function(node) {
			chains.push(node.data("pdb").split("|")[2]);
			prots.push(node.parent().data("label"));
		});

		edge.qtip({
			content:
				"<span class='tip' style='color: #7b241c;'>" +
				"<b>Predicted with <a href='http://www.russelllab.org/cgi-bin/tools/interprets.pl/interprets.pl'>InterPreTS</a></b>" +
				"</span><br>" +
				"<span class='tip'>"+
					"<span class='tipInP'>PDB Template ID</span> |  "+
					"<a href='https://www.rcsb.org/structure/"+pdb+"'><b>"+pdb+"</b> <i class='fas fa-external-link-alt fa-xs'></i></a><br>"+
					"<span class='tipInP'>PDB Template Chains</span> | "+chains.join(" - ")+"<br>"+
					"<span class='tipInP'>Interacting Proteins</span> | "+prots.join(" - ")+"<br>"+
					"<span class='tipInP'>Z-score</span> | "+z+"<br>"+
					"<span class='tipInP'>P-value</span> | "+pval+
				"</span>",
			 position: {
					my: "top center",
					at: "bottom center"
			 },
			 style: {
				 classes: "qtip-bootstrap",
					tip: {
						width: 20,
						height: 10
					}
				}
		});
	});
});
