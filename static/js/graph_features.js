"use strict";

// function makeEchart3() {
// 	var chart3 = echarts.init(document.getElementById('chart3'));
// 	var dom_sets = new Map();
// 	var dom_ints = cy.$("edge[role='DOM_interaction']");
// 	dom_ints.forEach(function(dom_int){
// 		if (dom_int.visible()) {
// 			var source = dom_int.source();
// 			var target = dom_int.target();
// 			if (source.data("label") in dom_sets) {
// 				dom_sets[source.data("label")].add(target.data("parent"));
// 			} else {
// 				dom_sets[source.data("label")] = new Set();
// 				dom_sets[source.data("label")].add(target.data("parent"));
// 			}
// 			if (target.data("label") in dom_sets) {
// 				dom_sets[target.data("label")].add(source.data("parent"));
// 			} else {
// 				dom_sets[target.data("label")] = new Set();
// 				dom_sets[target.data("label")].add(source.data("parent"));
// 			}
//
// 			// var doms = dom_int.connectedNodes();
// 			// doms.forEach(function(dom){
// 			// 	var dom_name = dom.data("label");
// 			// 	if (dom_name in dom_count){
// 			// 		dom_count[dom_name]+=1;
// 			// 	} else {
// 			// 		dom_count[dom_name]=1;
// 			// 	}
// 			// 	console.log(dom_name, dom_count[dom_name]);
// 			// });
// 		}
// 	});
// 	var dom_count = new Map();
// 	Object.keys(dom_sets).forEach(function(dom){
// 		console.log(dom+"dom: "+dom_sets[dom].size);
// 		dom_count[dom] = dom_sets[dom].size;
// 	});
// 	dom_sets.clear();
// 	// var dom_count_sorted = new Map([...Object.entries(dom_count)].sort((a, b) => a.value - b.value));
// 	// // dom_count.clear();
// 	// console.log("ads:"+dom_count_sorted);
// 	// dom_count_sorted.forEach(function(dom){
// 	// 	console.log(dom, dom_count_sorted.dom);
// 	// });
// 	var option = {
// 			title: {
// 					text: dom_ints.length+" domains"
// 			},
// 			grid: {
// 				left: '15%',
// 			},
// 			tooltip: {
// 				trigger: 'axis',
// 				axisPointer: {
// 					type: 'shadow'
// 				}
// 			},
// 			xAxis: {
// 				minInterval: 1
// 			},
// 			yAxis: {
// 					type: 'category',
// 					data: Object.keys(dom_count),
// 					axisLabel: {fontWeight: 'bold'}
// 			},
// 			series: [{
// 					name: 'Interactors',
// 					type: 'bar',
// 					color: ['#16A085'],
// 					data: Object.values(dom_count)
// 			}]
// 	};
// 	chart3.setOption(option);
// };

$(document).ready(function(){

	function toggle_nodes(element, node_description) {
		var nodes = cy.$(node_description);
		if (element.checked) {
			nodes.style("display", "element");
		} else {
			nodes.style("display", "none");
		}
	};
	
	function toggle_nodes_and_connectedEdges(element, node_description) {
		var nodes = cy.$(node_description);
		var edges = nodes.connectedEdges();
		if (element.checked) {
			nodes.style("display", "element");
			edges.style("display", "element");
		} else {
			nodes.style("display", "none");
			edges.style("display", "none");
		}
	};
	
	// Show/Hide proteins without interactions
	function showProt(label) {
		var node = cy.$("node[role='protein_main'][label='"+label+"']");
		var display = node.style("display");
		if (display=="none"){
			node.style("display", "element");
			// cy.center(node);
		} else {
			node.style("display", "none");
		}
	};
	
	function toggle_int_above_pval(role) {
		var pval = document.getElementById("pval_cutoff").value;
		cy.$("edge[role='"+role+"'][p_val<="+pval+"]").style("display", "element");
	}

	function toggle_cosmic_muts_by_count() {
		var min_count = document.getElementById("sample_min").value;
		cy.$("node[role='mod_cosmic'][tot_count>="+min_count+"]").style("display", "element");
		cy.$("node[role='mod_cosmic'][tot_count<"+min_count+"]").style("display", "none");
	}

	// Save original positions of all graph elements
	cy.nodes().forEach(function(n){
		var pos = n.position();
		n.data("orPos", {
			x: pos.x,
			y: pos.y
		});
	});

	// Set initial zoom-dependent size for protein labels
	var node = cy.$("node[role='protein_main']");
	var dim = 16/cy.zoom();
	var maxDim = Math.max(dim,30);
	node.style({"font-size": maxDim,
							"text-outline-width": maxDim/10,
							"border-width": maxDim/18,
							"text-margin-y": maxDim/-2.5
	});

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
	cy.on("doubleTap", "node[role='protein_main']", function(event) {
		var node = event.target;
		  cy.animate({
			fit: { eles: node, padding: 100}
		})
	});
  
	// FEATURE: Some nodes & edges' sizes will vary with zoom levels
	cy.on("render zoom", function(event) {
		var node = cy.$("node[role='protein_main']");
		var dim = 16/cy.zoom();
		var maxDim = Math.max(dim,30);
		node.style({"font-size": maxDim,
					"text-outline-width": maxDim/10,
					"border-width": maxDim/18,
					"text-margin-y": maxDim/-2.5
		});
	  // if (cy.zoom() <= 0.5){
	  // 	node.style("visibility", "hidden")
	  // }
	});

	// MOUSEOVER/MOUSEOUT changes
	cy.autoungrabify(true);
	cy.on("tapdragover tapdragout",
			"node[role='protein_main']", function(event) {
		var node = event.target;
		node.toggleClass("hl");
		node.connectedEdges().toggleClass("hl");
		node.neighborhood().toggleClass("hl2");
		if (event.type=="tapdragover"){
			cy.autoungrabify(false);
		} else if (event.type=="tapdragout"){
			cy.autoungrabify(true);
		}
		});

	cy.on("tapdragover tapdragout",
			"node[role='domain'], node[role='elm'], node[role='3dlm'], "+
			"node[role='iprets'], node[role ^='uni'], node[role ^='mod']", function(event) {
		var node = event.target;
		node.toggleClass("hl");
		node.connectedEdges().toggleClass("hl");
		node.neighborhood().toggleClass("hl2");
		// if (event.type=="tapdragover"){
		// 	node.ungrabify();
		// } else if (event.type=="tapdragout") {
		// 	node.grabify();
		// }
	});

	cy.on("mouseover mouseout", "edge[role$='_interaction']", function(event) {
		var edge = event.target;
		edge.toggleClass("hl");
		edge.connectedNodes().toggleClass("hl2");
		edge.connectedNodes().parent().toggleClass("hl2");
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
		["#toggle_box", "#toggle_doms", "#toggle_pp_int"].forEach(function(ele) {
			$(ele).prop("checked", true);
		});
		["#toggle_elms", "#toggle_confirmed_elms", "#toggle_iprets", "#toggle_dom_int",
		 "#toggle_idom_int", "#toggle_elmdom_int", "#toggle_elmdom_pred_int",
		 "#toggle_prets_int", "#toggle_phos", "#toggle_acet", "#toggle_uni_regions",
		 "#toggle_uni_variants", "#toggle_uni_mutagen", "#toggle_input_mut",
		 "#toggle_cosmic_mut"].forEach(function(ele) {
			$(ele).prop("checked", false);
		});

		["#elm_confirmed_container", "confirmed_dmi3d_container",
		 "#lo_slider_container"].forEach(function(ele) {
			$(ele).hide();
		});

		$("#lo_slider").prop("disabled", true);
	});

	// BUTTON: Center the graph on the page
  	$("#center").click(function(){
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
		var z = cy.zoom() + 0.2;
		cy.zoom(z);
	});
	$("#zoom_out").click(function(){
		var z = cy.zoom() - 0.2;
		cy.zoom(z);
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
	// BUTTON: Get snapshot as SVG
  	$("#dl_svg").click(function(filename) {
		var svgContent = cy.svg({scale: 1, full: true});
		var blob = new Blob([svgContent], {type:"image/svg+xml;charset=utf-8"});
		saveAs(blob, "graph.svg");
	});


	/// TOGGLE TOOLS:
  	// BUTTON: Toggle box surrounding the whole protein
	$("#toggle_box").click(function(){
		var eles = cy.$("node[role='protein_main']");
    	var eles2 = cy.$("node[role='protein_main']:selected");
    	var checked = document.getElementById("toggle_box").checked;
		if (checked) {
			eles.style("border-opacity", "1");
		} else {
			eles.style("border-opacity", "0");
			eles2.style("border-opacity", "1");
		}
  	});

 	// BUTTON: Toggle ALL protein elements
	$("#toggle_all_ele").change(function(){
		var checked = document.getElementById("toggle_all_ele").checked;
		if (checked) {
			$("#toggle_doms").prop("checked", true).change();
			$("#toggle_elms").prop("checked", true).change();
			$("#toggle_iprets").prop("checked", true).change();
		} else {
			$("#toggle_doms").prop("checked", false).change();
			$("#toggle_elms").prop("checked", false).change();
			$("#toggle_iprets").prop("checked", false).change();
		}
	});

  	// BUTTON: Toggle domains
  	$("#toggle_doms").change(function(){
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
			$("#lo_slider").prop("disabled", true);
			$("#lo_slider_container").hide();
			$("#toggle_elmdom_int").prop("checked", false);
			$("#toggle_elmdom_pred_int").prop("checked", false);
    	}
  	});

  	// BUTTON: Toggle elms
 	$("#toggle_elms").change(function(){
		var nodes = cy.$("node[role='elm']");
		var edges = nodes.connectedEdges();
    	var checked = document.getElementById("toggle_elms").checked;
    	if (checked) {
			nodes.style("display", "element");
			$("#elm_confirmed_container").show();
    	} else {
			nodes.style("display", "none");
			edges.style("display", "none");
			$("#elm_confirmed_container").hide();
			$("#toggle_confirmed_elms").prop("checked", false);
			$("#toggle_elmdom_int").prop("checked", false);
			$("#toggle_elmdom_pred_int").prop("checked", false);
		}
	});
	
	$("#toggle_confirmed_elms").change(function(){
		var nodes = cy.$("node[role='elm'][status!='TP']");
		var edges = nodes.connectedEdges();
		var checked = document.getElementById("toggle_confirmed_elms").checked;
		var checked2 = document.getElementById("toggle_elmdom_int").checked;
		if (checked) {
			nodes.style("display", "none");
			edges.style("display", "none");
		} else {
			nodes.style("display", "element");
			if (checked2) {
				edges.style("display", "element")
			}
		}
	});

	// BUTTON: Toggle InterPreTs regions
	$("#toggle_iprets").change(function(){
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
	$("#toggle_all_int").change(function(){
		var nodes = cy.$("node[role='domain'], node[role='elm'], node[role='iprets']");
		var edges = nodes.connectedEdges();
		var checked = document.getElementById("toggle_all_int").checked;
		if (checked) {
			nodes.style("display", "element");
			edges.style("display", "element");
			["#toggle_pp_int", "#toggle_dom_int", "#toggle_idom_int",
			"#toggle_elmdom_int", "#toggle_elmdom_pred_int",
			"#toggle_prets_int"].forEach(function(ele) {
				$(ele).prop("checked", true).change();
			});
		} else {
			edges.style("display", "none");
			["#toggle_pp_int", "#toggle_dom_int", "#toggle_idom_int",
			"#toggle_elmdom_int", "#toggle_elmdom_pred_int",
			"#toggle_prets_int"].forEach(function(ele) {
				$(ele).prop("checked", false).change();
			});
		}
	});

  	// BUTTON: Toggle protein-protein interactions
  	$("#toggle_pp_int").change(function(){
    	var edges = cy.$("edge[role='prot_prot_interaction']");
    	var checked = document.getElementById("toggle_pp_int").checked;
		if (checked) {
      		edges.style("display", "element");
		} else {
     	 	edges.style("display", "none");
    	}
  	});

	// BUTTON: Toggle domain-domain interactions (3did)
	$("#toggle_dom_int").change(function(){
		var edges = cy.$("edge[role='DOM_interaction']");
		var checked = document.getElementById("toggle_dom_int").checked;
		if (checked) {
			toggle_int_above_pval("DOM_interaction");
			$("#toggle_doms").prop("checked", true).change();
			// makeEchart3();
		} else {
			edges.style("display", "none");
			// makeEchart3();
		}
	});

	// BUTTON: Toggle domain-domain interactions (predicted)
	$("#toggle_idom_int").change(function(){
		var edges = cy.$("edge[role='iDOM_interaction']");
		var checked = document.getElementById("toggle_idom_int").checked;
		if (checked) {
			toggle_int_above_pval("iDOM_interaction");
			$("#toggle_doms").prop("checked", true).change();
			$("#lo_slider").prop("disabled", false);
			$("#lo_slider_container").show();
		} else {
			edges.style("display", "none");
			$("#lo_slider").prop("disabled", true);
			$("#lo_slider_container").hide();
		}
	});

	// BUTTON: Toggle elm-domain interactions
	$("#toggle_elmdom_int").change(function(){
		var edges = cy.$("edge[role='ELM_interaction']");
		var checked = document.getElementById("toggle_elmdom_int").checked;
		if (checked) {
			toggle_int_above_pval("ELM_interaction");
			$("#toggle_elms").prop("checked", true).change();
		} else {
			edges.style("display", "none");
			$("#toggle_elms").prop("checked", false).change();
		}
	});

	// BUTTON: Toggle DMI (3did)
	$("#toggle_dmi3d_int").change(function(){
		var edges = cy.$("edge[role='DMI_interaction']");
		var nodes = cy.$("node[role='3dlm']");
		var checked = this.checked;
		if (checked) {
			nodes.style("display", "element");
			edges.style("display", "element");
			$("#confirmed_dmi3d_container").show();
			$("#confirmed_dmi3d_toggle").prop("checked", true).change();
		} else {
			nodes.style("display", "none");
			edges.style("display", "none");
			$("#confirmed_dmi3d_container").hide();
		}
	});

	$("#confirmed_dmi3d_toggle").change(function(){
		var nodes = cy.$("node[role='3dlm'][status!='TP']");
		var edges = nodes.connectedEdges();
		var checked = this.checked;
		if (checked) {
			nodes.style("display", "none");
			edges.style("display", "none");
		} else {
			nodes.style("display", "element");
			edges.style("display", "element")
		}
	});

	// BUTTON: Toggle elm-domain interactions (predicted)
	$("#toggle_elmdom_pred_int").change(function(){
		var edges = cy.$("edge[role='iELM_interaction']");
		var checked = document.getElementById("toggle_elmdom_pred_int").checked;
		if (checked) {
			edges.style("display", "element");
			$("#toggle_elms").prop("checked", true).change();
		} else {
			edges.style("display", "none");
		}
 	});

	// BUTTON: Toggle InterPreTS interactions
	$("#toggle_prets_int").change(function(){
		var edges = cy.$("edge[role='INT_interaction']");
		var checked = document.getElementById("toggle_prets_int").checked;
		if (checked) {
			edges.style("display", "element");
			$("#toggle_iprets").prop("checked", true).change();
		} else {
			edges.style("display", "none");
		}
	});

	// BUTTON: Toggle ALL protein modifications
	$("#toggle_all_mod").change(function(){
		var checked = document.getElementById("toggle_all_mod").checked;
		if (checked) {
			["#toggle_input_mut",
			"#toggle_phos", "#toggle_acet"].forEach(function(ele) {
				$(ele).prop("checked", true).change();
			});
			toggle_cosmic_muts_by_count();
		} else {
			["#toggle_input_mut", "#toggle_cosmic_mut",
			"#toggle_phos", "#toggle_acet"].forEach(function(ele) {
				$(ele).prop("checked", false).change();
			});
		}
	});

	// BUTTON: Toggle mutations/PTMs
	$("#toggle_input_mut").change(function(){
		toggle_nodes(this, "node[role='mod_input']");
	});
	$("#toggle_cosmic_mut").change(function(){
		var nodes = cy.$("node[role='mod_cosmic']");
		if (this.checked) {
			toggle_cosmic_muts_by_count();
		} else {
			nodes.style("display", "none");
		}	
	});
	$("#toggle_phos").change(function(){
		toggle_nodes(this, "node[role='mod_phos']");
	});
	$("#toggle_acet").change(function(){
		toggle_nodes(this, "node[role='mod_acet']");
	});

	// BUTTON: Toggle UniProt Regions
	$("#toggle_uni_regions").change(function(){
		toggle_nodes_and_connectedEdges(this, "node[role='uni_region']");
	});

	// BUTTON: Toggle UniProt Variants
	$("#toggle_uni_variants").change(function(){
		toggle_nodes_and_connectedEdges(this, "node[role='uni_var']");
	});

	// BUTTON: Toggle UniProt Mutagenesis
	$("#toggle_uni_mutagen").change(function(){
		toggle_nodes_and_connectedEdges(this, "node[role='uni_mtg']");
	});

	// SLIDERS
	var sliderLO = document.getElementById("lo_slider");
	var valueLO = document.getElementById("lo_value");
	valueLO.innerHTML = sliderLO.value;
	sliderLO.oninput = function() {
	  	valueLO.innerHTML = this.value;
		var edges1 = cy.$("edge[role='iDOM_interaction'][lo>="+this.value+"]");
		edges1.style("display", "element");
		var edges2 = cy.$("edge[role='iDOM_interaction'][lo<"+this.value+"]");
		edges2.style("display", "none");
		$("#toggle_idom_int").prop("checked", true);
	};
	// var sliderPval = document.getElementById("pval_slider");
	// var valuePval = document.getElementById("pval_value");
	// valuePval.innerHTML = sliderPval.value;
	// sliderPval.oninput = function() {
	//   valuePval.innerHTML = this.value;
	// 	var edges1 = cy.$("edge[role='iDOM_interaction'][p_val<="+this.value+"]");
	// 	edges1.style("display", "element");
	// 	var edges2 = cy.$("edge[role='iDOM_interaction'][p_val>"+this.value+"]");
	// 	edges2.style("display", "none");
	// 	$("#toggle_idom_int").prop("checked", true);
	// };

	var cutoffPval = document.getElementById("pval_cutoff");
	cutoffPval.oninput = function() {
		if (document.getElementById("toggle_dom_int").checked) {
			cy.$("edge[role='DOM_interaction'][p_val<="+this.value+"]").style("display", "element");
			cy.$("edge[role='DOM_interaction'][p_val>"+this.value+"]").style("display", "none");
		}
		if (document.getElementById("toggle_idom_int").checked) {
			cy.$("edge[role='iDOM_interaction'][p_val<="+this.value+"]").style("display", "element");
			cy.$("edge[role='iDOM_interaction'][p_val>"+this.value+"]").style("display", "none");
		}
		if (document.getElementById("toggle_elmdom_int").checked) {
			cy.$("edge[role='ELM_interaction'][p_val<="+this.value+"]").style("display", "element");
			cy.$("edge[role='ELM_interaction'][p_val>"+this.value+"]").style("display", "none");
		}
	}


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

	// QTIPs on node/edge click
	// cy.nodes().on("tap", function( event ){
	cy.on("click", "node[role='protein_main']", function(event) {
		var node = event.target;
		var role = node.data("role");
		if("blast" in node.data()){
			var label = node.data("label");
			var length = node.data("length");
			var link = node.data("link");
			var blast = node.data("blast");
			var blast_hits = "<span class='tipProt'>Blast hits:</span><br>\n";
			blast_hits += "<div style='margin: 0 auto;'>";
			blast.forEach(function(hit) {
				var coor = hit.split("|")[0];
				var acc = hit.split("|")[1];
				var gene = hit.split("|")[2];
				var info = hit.split("|")[3];
				blast_hits+="\t<span style='float:left;width:80px;'><b>"+coor+"</b></span>";
				blast_hits+=" | <a href='https://www.uniprot.org/uniprot/"+acc+"'>"+acc+"</a>";
				blast_hits+=" <b>"+gene+"</b> "+info+"<br>\n";
			});
			blast_hits += "</div>"
			var qtip_content = "<span class='tip'>User-input sequence<br>\n"+
							"\t<span class='tipProt'>Label</span> | <b>"+label+"</b><br>\n" +
							"\t<span class='tipProt'>Length</span> | "+
							"<a href='"+link+"' target='_blank'_>" +
							length+" AA <i class='fas fa-external-link-alt fa-xs'></i></a><br>"+
							blast_hits+
							"</span>";
		} else {
			var gene = node.data("label");
			var des = node.data("des");
			var acc = node.data("protein");
			var length = node.data("length");
			var qtip_content = "<span class='tip'>"+
							"<span class='tipProt'>Gene</span> | <b>"+gene+"</b><br>" +
							"<span class='tipProt'>Protein</span> | <b>"+des+"</b><br>" +
							"<span class='tipProt'>Accession</span> | " +
							"<a href='https://www.uniprot.org/uniprot/"+acc+"'>" +
							acc+" <i class='fas fa-external-link-alt fa-xs'></i>" +
							"</a><br>" +
							"<span class='tipProt'>Length</span> | "+
							"<a href='https://www.uniprot.org/uniprot/"+acc+".fasta'>" +
							length+" AA <i class='fas fa-external-link-alt fa-xs'></i>"+
							"</span>";
		}
		node.qtip({
			content: qtip_content,
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
	});

	cy.on("click", "node[role='domain']", function(event) {
		var node = event.target;
		var name = node.data("label");
		var acc = node.data("acc");
		var tp = node.data("type");
		var des = node.data("des");
		var start = node.data("start");
		var end = node.data("end");
		var prot = node.data("protein");
		var e_val = Number.parseFloat(node.data("e_val")).toExponential(2);
		var qtip_content = "<span class='tip' style='color: #074987;'>" +
							"<b>Source: <a href='https://pfam.xfam.org'>Pfam</a></b>" +
							"</span><br>" +
							"<span class='tip'>" +
								"<span class='tipPfam'>Family</span> | " +
									"<b><i>"+name+"</i></b> (<a href='https://pfam.xfam.org/family/"+acc+"'>"+acc+" <i class='fas fa-external-link-alt fa-xs'></i></a>)<br>" +
								"<span class='tipPfam'>Type</span> | <b>"+tp+"</b><br>"+
								"<span class='tipPfam'>Description</span> | <b>"+des+"</b><br>"+
								"<span class='tipPfam'>Start - End</span> | " +
									"<b>"+start+"</b> - <b>"+end+"</b>" +
									" (<a href='https://pfam.xfam.org/protein/"+prot+"'>"+prot+" <i class='fas fa-external-link-alt fa-xs'></i></a>)<br>"+
								"<span class='tipPfam'>E-value</span> | "+e_val+
							"</span>";
		node.qtip({
			content: qtip_content,
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
		var node 	= event.target;
		var id 		= node.data("id");
		var label 	= node.data("label");
		var acc 	= node.data("acc");
		var name 	= node.data("name");
		var regex 	= node.data("regex");
		var start 	= node.data("start");
		var end 	= node.data("end");
		var seq 	= node.data("seq");
		var prot 	= node.data("protein");
		var status	= node.data("status");
		var i = 0;
		var tot = 0;
		var des0 = [];
		node.data("des").split(" ").forEach(function(char) {
			i += char.length+1;
			tot += char.length+1;
			if (i>=60 && tot<node.data("des").length) {
				des0.push(char+"<br>");
				i = 0;
			} else {
				des0.push(char);
			}
		});
		var des = des0.join(" ");
		var phos = "";
		var clas = "";
		if (node.data("phos").includes("yes")) {
			phos = "<span class='pos'><i class='fas fa-plus-circle'></i> "+
							"Required phosphosites found in motif</span>";
			clas = "pos";
		} else if (node.data("phos").includes("no")) {
			phos = "<span class='neg'><i class='fas fa-minus-circle'></i></i> "+
							"Required phosphosites not found in motif</span>";
			clas = "neg";
		}

		var edit_seq = "";
		var pos = start;
		seq.split("").forEach(function(res) {
			if (node.data("phos").includes(pos)) {
				edit_seq += "(<span class='"+clas+"'><b>"+res+"</b></span>)";
			} else {
				edit_seq += res;
			}
			pos += 1;
		});

		var qtip_content = 	"<span class='tip' style='color: #7f7c7b;'>\n" +
								"<b>Source: <a href='https://elm.eu.org'>ELM</a></b>\n" +
							"</span><br>\n" +
							"<span class='tip'>\n"+
								"<span class='tipELM'>Identifier</span> | "+
									"<a href='http://elm.eu.org/elms/"+label+"'>"+
											"<span style='color: #11249b;'>"+label+
											"</span> <i class='fas fa-external-link-alt fa-xs'></i>"+
									"</a><br>\n"+
								"<span class='tipELM'>Accession</span> | "+acc+"<br>\n"+
								"<span class='tipELM'>Class</span> | <b>"+name+"</b><br>\n"+
								"<span class='tipELM'>Description</span> | "+des+"<br>\n"+
								"<span class='tipELM'>Subsequence</span> | " +
									"<b>"+start+"</b> - <i>"+edit_seq+"</i> - <b>"+end+"</b>" +
									" (<a href='http://elm.eu.org/instances/"+label+"/"+prot+"/'>"+prot+
									" <i class='fas fa-external-link-alt fa-xs'></i></a>)<br>\n"+
								"<span class='tipELM'>Status</span> | "+status+"<br>"+
								"<span class='tipELM'>Observations</span> | "+phos+"<br>"+
							"</span><br>\n"+

							"<div class='row'>\n"+
								"<div class='col-sm-4 text-center'>\n"+
									"<button class='txt-btn' onclick=removeThisEle(\""+id+"\") >\n"+
										"Remove\n"+
									"</button>\n"+
								"</div>\n"+
								"<div class='col-sm-4 text-center'>\n"+
									"<button class='txt-btn' onclick=removeAllinProt(\""+label+"\",\""+prot+"\") >\n"+
										"Remove all in this protein\n"+
									"</button>\n"+
								"</div>\n"+
								"<div class='col-sm-4 text-center'>\n"+
									"<button class='txt-btn' onclick=removeAll(\""+label+"\") >\n"+
										"Remove all in network\n"+
									"</button>\n"+
								"</div>\n"+
							"</div>";

		node.qtip({
			content: qtip_content,
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
		var qtip_content = "<span class='tip' style='color: #7b241c;'>" +
							"<b>Predicted with <a href='http://www.russelllab.org/cgi-bin/tools/interprets.pl/interprets.pl'>InterPreTS</a></b>" +
							"</span><br>" +
							"<span class='tip'>"+
								"<span class='tipInP'>PDB Template ID</span> | " +
								"<a href='https://www.rcsb.org/structure/"+pdb+"'><b>"+pdb+"</b> <i class='fas fa-external-link-alt fa-xs'></i></a> chain "+chain+"<br>"+
								"<span class='tipInP'>Alignment</span> |  "+
								"Protein/<b>"+start+"-"+end+"</b>; Template/"+pdb_start+"-"+pdb_end+"<br>"+
								"<span class='tipInP'>Alignment score</span> | BLAST e-val=<i>"+ev+"</i>, <i>"+pcid+"%</i> id<br>"+
							"</span>";
		node.qtip({
			content: qtip_content,
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

	cy.on("click","node[role='mod_cosmic']", function(event) {
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
		var qtip_content = "<span class='tip' style='color: #28446f;'>" +
								"<b><a href='https://cancer.sanger.ac.uk/cosmic'>COSMIC <i class='fas fa-external-link-alt fa-xs'></i></a></b><br>" +
								"<a href="+gene_link+">View gene "+gene+" <i class='fas fa-external-link-alt fa-xs'></i></a>"+
							"</span><br>\n"+html;

		node.qtip({
			content: qtip_content,
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
		var qtip_content = "<span class='tip' style='color: #33a1c2;'>" +
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
							"</span><br>";

		edge.qtip({
			content: qtip_content,
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
			doms.push("<span style='color: "+node.data("color")+";'>"+
							"<b>"+node.data("label")+"</b></span>"
							+" (<a href='https://pfam.xfam.org/family/"+node.data("acc")+"'>"+node.data("acc")+" <i class='fas fa-external-link-alt fa-xs'></i></a>)");
			prots.push(node.parent().data("label"));
		});
		var ds    = edge.data("ds");
		var pdb_n = edge.data("pdb_n");
		var p_val = Number.parseFloat(edge.data("p_val")).toExponential(2);

		var qtip_content = 	"<span class='tip' style='color: #16a085;'>" +
								"<b>Domain-Domain Interaction</b>" +
							"</span><br>" +
							"<span class='tip'>"+
								"<span class='tip3did'>Interacting Domains</span> | "+				 		
									doms.join(" - ")+
							"</span><br>"+
							"<span class='tip'>"+
								"<span class='tip3did'>Interacting Proteins</span> | "+prots.join(" - ")+
							"</span><br>"+
							"<span class='tip'>"+
								"<span class='tip3did'>Sources</span> | "+ds+
							"</span><br>"+
							"<span class='tip'>"+
								"<span class='tip3did'># PDB structures</span> | "+
								"<a href='https://3did.irbbarcelona.org/dispatch.php?type=interaction&type1=domain&type2=domain&value1="+pfams[0]+"&value2="+pfams[1]+"'>"+
									pdb_n+" <i class='fas fa-external-link-alt fa-xs'></i>"+
								"</a>"+
							"</span><br>"+
							"<span class='tip'>"+
								"<span class='tip3did'>Interaction P-value</span> | "+p_val+
							"</span>";
		edge.qtip({
			content: qtip_content,
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
		var ds = edge.data("ds");
		var lo = Number(edge.data("lo")).toFixed(2);
		var p_val = Number.parseFloat(edge.data("p_val")).toExponential(2);
		var qtip_content = 	"<span class='tip' style='color: #d4ac0d;'>" +
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
							"</span><br>"+
							"<span class='tip'>"+
								"<span class='tipIdom'>Interaction P-value</span> | "+p_val+
							"</span>";
		edge.qtip({
			content: qtip_content,
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
		var p_val = Number.parseFloat(edge.data("p_val")).toExponential(2);
		var qtip_content = 	"<span class='tip' style='color: #b95db9;'>" +
								"Source: <a href='http://elm.eu.org/downloads.html#interactions'><b>ELM</b> (interactions)</a></b>" +
							"</span><br>" +
							"<span class='tip'>"+
								"<span class='tipELMint'>Interacting Elements</span> | "+
									doms.join(" - ")+
								"</span><br>"+
							"<span class='tipELMint'>Interacting Proteins</span> | "+prots.join(" - ")+
							"</span><br>"+
							"<span class='tipELMint'>Interaction P-value</span> | "+p_val+
							"</span>";
		edge.qtip({
			content: qtip_content,
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
		var pdb = edge.data("pdb");
		var z = edge.data("z-score");
		var pval = edge.data("p-value");
		var chains = [];
		var prots = [];
		nodes.forEach(function(node) {
			chains.push(node.data("pdb").split("|")[2]);
			prots.push(node.parent().data("label"));
		});
		var qtip_content = 	"<span class='tip' style='color: #7b241c;'>" +
								"<b>Predicted with <a href='http://www.russelllab.org/cgi-bin/tools/interprets.pl/interprets.pl'>InterPreTS</a></b>" +
							"</span><br>" +
							"<span class='tip'>"+
								"<span class='tipInP'>PDB Template ID</span> |  "+
								"<a href='https://www.rcsb.org/structure/"+pdb+"'><b>"+pdb+"</b> <i class='fas fa-external-link-alt fa-xs'></i></a><br>"+
								"<span class='tipInP'>PDB Template Chains</span> | "+chains.join(" - ")+"<br>"+
								"<span class='tipInP'>Interacting Proteins</span> | "+prots.join(" - ")+"<br>"+
								"<span class='tipInP'>Z-score</span> | "+z+"<br>"+
								"<span class='tipInP'>P-value</span> | "+pval+
							"</span>";
		edge.qtip({
			content: qtip_content,
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
