"use strict";

$(document).ready(function(){

	// BUTTON: reset graphic
	// 1- Save original positions:
	cy.nodes().forEach(function(n){
		var positions = n.position();

		n.data("orgPos", {
			x: positions.x,
			y: positions.y
			});
	});
	// 2 - Reset on click
	$("#reset").click(function(){
	 cy.nodes().forEach(function(n){
	 var p = n.data('orgPos');
	 n.position({ x: p.x, y: p.y });
	 });
	 cy.reset();
	 cy.fit();
	});

	// BUTTON: Center the graph on the page
  $("#center").click(function(){
    // cy.center();
		cy.fit();
  });

  // BUTTON: Enables/Disables the option of draggin the graph elements
  var lock_action = 0;
  $("#lock").click(function(){
    if ( lock_action == 0 ) {
      cy.autolock(true);
      $("#lock i").toggleClass("fa-unlock fa-lock");
      lock_action = 1;
    } else {
      cy.autolock(false);
      $("#lock i").toggleClass("fa-lock fa-unlock");
      lock_action = 0;
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
    var eles = cy.$('node[role="whole"]');
    var eles2 = cy.$('node[role="whole"]:selected');
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
		var dom_nodes = cy.$('node[role="domain"]');
		var dom_edges = dom_nodes.connectedEdges();
		var elm_nodes = cy.$('node[role="elm"]');
		var elm_edges = elm_nodes.connectedEdges();
		var ip_nodes = cy.$('node[role="iprets"]');
		var ip_edges = ip_nodes.connectedEdges();
		var checked = document.getElementById("toggle_all_ele").checked;
		if (checked) {
			dom_nodes.style("visibility", "visible");
			elm_nodes.style("visibility", "visible");
			ip_nodes.style("visibility", "visible");
			$("#toggle_doms").prop("checked", true);
			$("#toggle_elms").prop("checked", true);
			$("#toggle_iprets").prop("checked", true);
		} else {
			dom_nodes.style("visibility", "hidden");
			dom_edges.style("visibility", "hidden");
			elm_nodes.style("visibility", "hidden");
			elm_edges.style("visibility", "hidden");
			ip_nodes.style("visibility", "hidden");
			ip_edges.style("visibility", "hidden");
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
    var nodes = cy.$('node[role="domain"]');
		var edges = nodes.connectedEdges();
    var checked = document.getElementById("toggle_doms").checked;
    if (checked) {
			nodes.style("visibility", "visible");
    } else {
			nodes.style("visibility", "hidden");
			edges.style("visibility", "hidden");
			$("#toggle_dom_int").prop("checked", false);
			$("#toggle_idom_int").prop("checked", false);
			$("#toggle_elmdom_int").prop("checked", false);
    }
  });

  // BUTTON: Toggle elms
  $("#toggle_elms").click(function(){
    var nodes = cy.$('node[role="elm"]');
		var edges = nodes.connectedEdges();
    var checked = document.getElementById("toggle_elms").checked;
    if (checked) {
      nodes.style("visibility", "visible");
    } else {
			nodes.style("visibility", "hidden");
			edges.style("visibility", "hidden");
			$("#toggle_elmdom_int").prop("checked", false);
    }
  });

  // BUTTON: Toggle InterPreTs regions
  $("#toggle_iprets").click(function(){
    var nodes = cy.$('node[role="iprets"]');
		var edges = nodes.connectedEdges();
    var checked = document.getElementById("toggle_iprets").checked;
    if (checked) {
      nodes.style("visibility", "visible");
    } else {
			nodes.style("visibility", "hidden");
			edges.style("visibility", "hidden");
			$("#toggle_prets_int").prop("checked", false);
    }
  });

  // BUTTON: Toggle ALL interactions
	$("#toggle_all_int").click(function(){
		var prot_edges = cy.$('edge[role="prot_prot_interaction"]');
		var user_edges = cy.$('edge[role="user_interaction"]');
		var dom_nodes = cy.$('node[role="domain"]');
		var dom_edges = dom_nodes.connectedEdges();
		var elm_nodes = cy.$('node[role="elm"]');
		var elm_edges = elm_nodes.connectedEdges();
		var ip_nodes = cy.$('node[role="iprets"]');
		var ip_edges = ip_nodes.connectedEdges();
		var checked = document.getElementById("toggle_all_int").checked;
		if (checked) {
			prot_edges.style("visibility", "visible");
			user_edges.style("visibility", "visible");
			dom_nodes.style("visibility", "visible");
			dom_edges.style("visibility", "visible");
			elm_nodes.style("visibility", "visible");
			elm_edges.style("visibility", "visible");
			ip_nodes.style("visibility", "visible");
			ip_edges.style("visibility", "visible");
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
			prot_edges.style("visibility", "hidden");
			user_edges.style("visibility", "hidden");
			dom_edges.style("visibility", "hidden");
			elm_edges.style("visibility", "hidden");
			ip_edges.style("visibility", "hidden");
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
    var edges = cy.$('edge[role="user_interaction"]');
    var checked = document.getElementById("toggle_user_int").checked;
		if (checked) {
      edges.style("visibility", "visible");
    } else {
      edges.style("visibility", "hidden");
    }
  });

  // BUTTON: Toggle protein-protein interactions
  $("#toggle_pp_int").click(function(){
    var edges = cy.$('edge[role="prot_prot_interaction"]');
    var checked = document.getElementById("toggle_pp_int").checked;
		if (checked) {
      edges.style("visibility", "visible");
    } else {
      edges.style("visibility", "hidden");
    }
  });

  // BUTTON: Toggle domain-domain interactions (3did)
  $("#toggle_dom_int").click(function(){
    var edges = cy.$('edge[role="DOM_interaction"]');
		var nodes = edges.connectedNodes();
    var checked = document.getElementById("toggle_dom_int").checked;
		if (checked) {
      edges.style("visibility", "visible");
			nodes.style("visibility", "visible");
			$("#toggle_doms").prop("checked", true);
    } else {
      edges.style("visibility", "hidden");
    }
  });

  // BUTTON: Toggle domain-domain interactions (statistical prediction)
  $("#toggle_idom_int").click(function(){
    var edges = cy.$('edge[role="iDOM_interaction"]');
		var nodes = edges.connectedNodes();
    var checked = document.getElementById("toggle_idom_int").checked;
		if (checked) {
      edges.style("visibility", "visible");
			nodes.style("visibility", "visible");
			$("#toggle_doms").prop("checked", true);
    } else {
      edges.style("visibility", "hidden");
    }
  });

  // BUTTON: Toggle elm-domain interactions
  $("#toggle_elmdom_int").click(function(){
    var edges = cy.$('edge[role="ELM_interaction"]');
		var nodes = edges.connectedNodes();
    var checked = document.getElementById("toggle_elmdom_int").checked;
		if (checked) {
      edges.style("visibility", "visible");
			nodes.style("visibility", "visible");
			$("#toggle_elms").prop("checked", true);
    } else {
      edges.style("visibility", "hidden");
    }
  });

  // BUTTON: Toggle InterPreTS interactions
  $("#toggle_prets_int").click(function(){
    var edges = cy.$('edge[role="INT_interaction"]');
		var nodes = edges.connectedNodes();
    var checked = document.getElementById("toggle_prets_int").checked;
		if (checked) {
      edges.style("visibility", "visible");
			nodes.style("visibility", "visible");
			$("#toggle_iprets").prop("checked", true);
    } else {
      edges.style("visibility", "hidden");
    }
  });

	// BUTTON: Toggle ALL protein modifications
	$("#toggle_all_mod").click(function(){
		var phos = cy.$('node[role="phosphorylation"]');
		var ace = cy.$('node[role="acetylation"]');
		var input_mut = cy.$('node[role="input_mut"]');
		var cosmic_mut = cy.$('node[role="cosmic_mut"]');
		var checked = document.getElementById("toggle_all_mod").checked;
		if (checked) {
			phos.style("visibility", "visible");
			ace.style("visibility", "visible");
			input_mut.style("visibility", "visible");
			cosmic_mut.style("visibility", "visible");
			$("#toggle_phos").prop("checked", true);
			$("#toggle_acet").prop("checked", true);
			$("#toggle_input_mut").prop("checked", true);
			$("#toggle_cosmic_mut").prop("checked", true);
		} else {
			phos.style("visibility", "hidden");
			ace.style("visibility", "hidden");
			input_mut.style("visibility", "hidden");
			cosmic_mut.style("visibility", "hidden");
			$("#toggle_phos").prop("checked", false);
			$("#toggle_acet").prop("checked", false);
			$("#toggle_input_mut").prop("checked", false);
			$("#toggle_cosmic_mut").prop("checked", true);
		}
	});

	// BUTTON: Toggle mutations
	$("#toggle_input_mut").click(function(){
		var eles = cy.$('node[role="input_mut"]');
		var checked = document.getElementById("toggle_input_mut").checked;
		if (checked) {
			eles.style("visibility", "visible");
		} else {
			eles.style("visibility", "hidden");
		}
	});

	$("#toggle_cosmic_mut").click(function(){
		var eles = cy.$('node[role="cosmic_mut"]');
		var checked = document.getElementById("toggle_cosmic_mut").checked;
		if (checked) {
			eles.style("visibility", "visible");
		} else {
			eles.style("visibility", "hidden");
		}
	});

	// BUTTON: Toggle phosphorylations
  $("#toggle_phos").click(function(){
    var eles = cy.$('node[role="phosphorylation"]');
    var checked = document.getElementById("toggle_phos").checked;
		if (checked) {
      eles.style("visibility", "visible");
    } else {
      eles.style("visibility", "hidden");
    }
  });

  // BUTTON: Toggle acetylations
  $("#toggle_acet").click(function(){
    var eles = cy.$('node[role="acetylation"]');
    var checked = document.getElementById("toggle_acet").checked;
    if (checked) {
      eles.style("visibility", "visible");
    } else {
      eles.style("visibility", "hidden");
    }
  });

	/// DOWNLOAD BUTTONS
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

  // BUTTON: Download graph as JSON
  $("#dl_json").click(function(){
	  var jsonBlob = new Blob([ JSON.stringify( cy.json() ) ], { type: 'application/javascript;charset=utf-8' });
	  saveAs( jsonBlob, 'graph.json' );
  });

	// Create double-Tap event
	var tappedBefore;
	var tappedTimeout;
	cy.on('tap', function(event) {
	  var tappedNow = event.target;
	  if (tappedTimeout && tappedBefore) {
	    clearTimeout(tappedTimeout);
	  }
	  if(tappedBefore === tappedNow) {
	    tappedNow.trigger('doubleTap');
	    tappedBefore = null;
	  } else {
	    tappedTimeout = setTimeout(function(){ tappedBefore = null; }, 300);
	    tappedBefore = tappedNow;
	  }
	});

	// FEATURE: Double-Tap on protein node to zoom-in
	cy.on('doubleTap', 'node[role=\"whole\"]', function(event) {
	  var node = event.target;
		cy.animate({
			fit: { eles: node,
						 padding: 100}
		})
	});

	// FEATURE: Some nodes & edges' sizes will vary with zoom levels
	cy.on('render zoom', function(event) {
		var node = cy.$('node[role="whole"]');
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

	cy.autoungrabify(true);

	// MOUSEOVER/MOUSEOUT changes
  cy.on('tapdragover tapdragout','node[role=\"whole\"]', function(event) {
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

	cy.on('tapdragover tapdragout',
				'node[role=\"domain\"], node[role=\"elm\"], node[role=\"iprets\"], '+
				'node[role=\"phosphorylation\"], node[role=\"acetylation\"], '+
				'node[role=\"input_mut\"], node[role=\"cosmic_mut\"]', function(event) {
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

	cy.on('mouseover mouseout','edge[role=\"prot_prot_interaction\"]', function(event) {
		var edge = event.target;
		edge.toggleClass("highlight");
		edge.connectedNodes().toggleClass("highlight2");
	});

	cy.on('mouseover mouseout','edge[role=\"DOM_interaction\"]', function(event) {
		var edge = event.target;
		edge.toggleClass("highlight");
		edge.connectedNodes().toggleClass("highlight2");
		edge.connectedNodes().parent().toggleClass("highlight2");
	});

	cy.on('mouseover mouseout','edge[role=\"iDOM_interaction\"]', function(event) {
		var edge = event.target;
		edge.toggleClass("highlight");
		edge.connectedNodes().toggleClass("highlight2");
		edge.connectedNodes().parent().toggleClass("highlight2");
	});

	cy.on('mouseover mouseout','edge[role=\"ELM_interaction\"]', function(event) {
		var edge = event.target;
		edge.toggleClass("highlight");
		edge.connectedNodes().toggleClass("highlight2");
		edge.connectedNodes().parent().toggleClass("highlight2");
	});

	cy.on('mouseover mouseout','edge[role=\"INT_interaction\"]', function(event) {
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
					"<a href='https://www.uniprot.org/uniprot/"+acc+"''>" +
					acc+" <i class='fas fa-external-link-alt fa-xs'></i>" +
					"</a><br>" +
					"<span class='tip'><span class='tipProt'>Length</span> | "+length+
					"</span>",
				position: {
					my: 'top center',
					at: 'bottom center'
				},
				style: {
					classes: 'qtip-bootstrap',
					tip: {
						width: 20,
						height: 10
					}
				},
				show: { event: 'directtap' }
			});

			this.trigger('directtap');
		}
	});

	cy.on('click','node[role=\"domain\"]', function(event) {
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
				my: 'top center',
				at: 'bottom center'
			},
			style: {
				classes: 'qtip-bootstrap',
				tip: {
					width: 20,
					height: 10
				}
			}
		});
	});

	cy.on('click','node[role=\"elm\"]', function(event) {
		var node = event.target;
		var name = node.data("label");
		var acc = node.data("acc");
		var des = node.data("des");
		var regex = node.data("regex");
		var start = node.data("start");
		var end = node.data("end");
		var seq = node.data("seq");
		var prot = node.data("protein");
		node.qtip({
			content:
				"<span class='tip' style='color: #7f7c7b;'>" +
					"<b>Source: <a href='https://elm.eu.org'>ELM</a></b>" +
				"</span><br>" +
				"<span class='tip'>"+
					"<span class='tipELM'>Identifier</span> | "+
					"<a href='http://elm.eu.org/elms/"+name+"'>"+name+" <i class='fas fa-external-link-alt fa-xs'></i></a><br>"+
					"<span class='tipELM'>Accession</span> | "+ acc + "<br>"+
					"<span class='tipELM'>Description</span> | <b>"+des+"</b><br>"+
					"<span class='tipELM'>Start - End</span> | " +
					"<b>"+start+"</b> - <b>"+end+"</b>" +
					" (<a href='http://elm.eu.org/instances/"+name+"/"+prot+"/'>"+prot+" <i class='fas fa-external-link-alt fa-xs'></i></a>)<br>"+
					"<span class=tipELM>Subsequence</span> | <i>"+seq+"</i>" +
				"</span>",
			position: {
				my: 'top center',
				at: 'bottom center'
			},
			style: {
				classes: 'qtip-bootstrap',
				tip: {
					width: 20,
					height: 10
				}
			}
		});
	});

	cy.on('click','node[role=\"iprets\"]', function(event) {
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
					"<a href='https://www.rcsb.org/structure/1JM7'><b>"+pdb+"</b> <i class='fas fa-external-link-alt fa-xs'></i></a> chain "+chain+"<br>"+
					"<span class='tipInP'>Alignment</span> |  "+
					"Protein/<b>"+start+"-"+end+"</b>; Template/"+pdb_start+"-"+pdb_end+"<br>"+
					"<span class='tipInP'>Alignment score</span> | BLAST e-val=<i>"+ev+"</i>, <i>"+pcid+"%</i> id<br>"+
				"</span>",
			position: {
				my: 'top center',
				at: 'bottom center'
			},
			style: {
				classes: 'qtip-bootstrap',
				tip: {
					width: 20,
					height: 10
				}
			}
		});
	});

	cy.on('click','node[role=\"cosmic_mut\"]', function(event) {
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
				my: 'top center',
				at: 'bottom center'
			},
			style: {
				classes: 'qtip-bootstrap',
				tip: {
					width: 20,
					height: 10
				}
			}
		});

	});

	cy.on('click','edge[role=\"prot_prot_interaction\"]', function(event) {
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
           my: 'top center',
           at: 'bottom center'
       	},
				style: {
					classes: 'qtip-bootstrap',
					tip: {
						width: 20,
						height: 10
			    }
			  }
			});
	});

	cy.on('click','edge[role=\"DOM_interaction\"]', function(event) {
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
           my: 'top center',
           at: 'bottom center'
       	},
				style: {
					classes: 'qtip-bootstrap',
					tip: {
						width: 20,
						height: 10
			    }
			  }
			});
	});

	cy.on('click','edge[role=\"iDOM_interaction\"]', function(event) {
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
					"<span class='tip' style='color: #d4ac0d;'>" +
						"<b>statistical prediction</b>" +
					"</span><br>" +
					"<span class='tip'>"+
						"<span class='tipIdom'>Interacting Domains</span> | "+doms.join(" - ")+
					"</span><br>"+
					"<span class='tip'>"+
						"<span class='tipIdom'>Interacting Proteins</span> | "+prots.join(" - ")+
					"</span>",
				position: {
           my: 'top center',
           at: 'bottom center'
       	},
				style: {
					classes: 'qtip-bootstrap',
					tip: {
						width: 20,
						height: 10
			    }
			  }
			});
	});

	cy.on('click','edge[role=\"ELM_interaction\"]', function(event) {
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
           my: 'top center',
           at: 'bottom center'
       	},
				style: {
					classes: 'qtip-bootstrap',
					tip: {
						width: 20,
						height: 10
			    }
			  }
			});
	});

	cy.on('click','edge[role=\"INT_interaction\"]', function(event) {
		var edge = event.target;
		var nodes = edge.connectedNodes();
		var pdb = edge.data("pdb").split("|")[1];
		var z = edge.data("z-score");
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
				"<a href='https://www.rcsb.org/structure/1JM7'><b>"+pdb+"</b> <i class='fas fa-external-link-alt fa-xs'></i></a><br>"+
				"<span class='tipInP'>PDB Template Chains</span> | "+chains.join(" - ")+"<br>"+
				"<span class='tipInP'>Interacting Proteins</span> | "+prots.join(" - ")+"<br>"+
				"<span class='tipInP'>Z-score</span> | "+z+
			"</span>",
			position: {
				 my: 'top center',
				 at: 'bottom center'
			},
			style: {
				classes: 'qtip-bootstrap',
				tip: {
					width: 20,
					height: 10
				}
			}
		});
	});
});
