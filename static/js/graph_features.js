"use strict";

$(document).ready(function(){

	//// CYTOSCAPE.JS FEATURES /////////////////////////////////////////////////

	//set original position
	cy.nodes().forEach(function(n){
		var positions = n.position();
		//~ console.log( positions );

		n.data("orgPos", {
			x: positions.x,
			y: positions.y
			});
		//console.log(n);
	});

	// // Lock PTM & Mutation nodes
	// cy.$('node[role="domain"]').ungrabify();
	// cy.$('node[role="elm"]').ungrabify();
	// cy.$('node[role="phosphorylation"]').grabify();
	// cy.$('node[role="acetylation"]').grabify();
	// cy.$('node[role="mutation"]').grabify();
	// cy.$('node[role="position"]').grabify();


  // MOUSEOVER
  cy.on('mouseover mouseout','node[role=\"whole\"]', function(event) {
    var node = event.target;
		// var edges = node.connectedEdges();
		// var neigh = node.neighborhood();
		node.toggleClass("highlight");
		node.connectedEdges().toggleClass("highlight");
		node.neighborhood().toggleClass("highlight2");
  });

	cy.on('mouseover','node[role=\"domain\"]', function(event) {
		var node = event.target;
		var edges = node.connectedEdges();
		var neigh = node.neighborhood();
		node.toggleClass("highlight");
		edges.toggleClass("highlight");
		neigh.toggleClass("highlight2");
		node.ungrabify();
	});

	cy.on('mouseout','node[role=\"domain\"]', function(event) {
		var node = event.target;
		var edges = node.connectedEdges();
		var neigh = node.neighborhood();
		node.toggleClass("highlight");
		edges.toggleClass("highlight");
		neigh.toggleClass("highlight2");
		node.grabify();
	});

	cy.on('mouseover mouseout','node[role=\"elm\"]', function(event) {
		var node = event.target;
		var edges = node.connectedEdges();
		var neigh = node.neighborhood();
		node.toggleClass("highlight");
		edges.toggleClass("highlight");
		neigh.toggleClass("highlight2");
	});


	cy.on('mouseover','node[role=\"phosphorylation\"]', function(event) {
		var node = event.target;
		node.toggleClass("highlight");
		node.ungrabify();
	});
	cy.on('mouseout','node[role=\"phosphorylation\"]', function(event) {
		var node = event.target;
		node.toggleClass("highlight");
		node.grabify();
	});

	cy.on('mouseover','node[role=\"acetylation\"]', function(event) {
		var node = event.target;
		node.toggleClass("highlight");
		node.ungrabify();
	});
	cy.on('mouseout','node[role=\"acetylation\"]', function(event) {
		var node = event.target;
		node.toggleClass("highlight");
		node.grabify();
	});

	cy.on('mouseover','node[role=\"mutation\"]', function(event) {
		var node = event.target;
		node.toggleClass("highlight");
		node.ungrabify();
	});
	cy.on('mouseout','node[role=\"mutation\"]', function(event) {
		var node = event.target;
		node.toggleClass("highlight");
		node.grabify();
	});

	cy.on('mouseover mouseout','edge[role=\"prot_prot_interaction\"]', function(event) {
		var edge = event.target;
		var nodes = edge.connectedNodes();
		edge.toggleClass("highlight");
		nodes.toggleClass("highlight2");
	});

	cy.on('mouseover mouseout','edge[role=\"DOM_interaction\"]', function(event) {
		var edge = event.target;
		var nodes = edge.connectedNodes();
		edge.toggleClass("highlight");
		nodes.toggleClass("highlight2");
	});

	cy.on('mouseover mouseout','edge[role=\"iDOM_interaction\"]', function(event) {
		var edge = event.target;
		var nodes = edge.connectedNodes();
		edge.toggleClass("highlight");
		nodes.toggleClass("highlight2");
	});

	cy.on('mouseover mouseout','edge[role=\"ELM_interaction\"]', function(event) {
		var edge = event.target;
		var nodes = edge.connectedNodes();
		edge.toggleClass("highlight");
		nodes.toggleClass("highlight2");
	});

	// CLICKS
	cy.on('click','edge[role=\"prot_prot_interaction\"]', function(event) {
			var edge = event.target;
			var nodes = edge.connectedNodes();
			var genes = []
			nodes.forEach(function(node) {
				genes.push(node.data("label"));
			});

			var ds = edge.data("ds");
			var links = edge.data("links").split(";");
			var all_links = []
			links.forEach(function(link) {
				all_links.push("<a href=\"https://thebiogrid.org/interaction/"+link+"\">"+link+"</a>");
			});

			edge.qtip({
				content: "<b>"+ds+" interaction</b><br>"+
								 "<b>"+genes.join(" - ")+"</b><br>"+
								 all_links.join(", "),
				position: {
           my: 'top center',
           at: 'bottom center'
       	},

				style: {
					classes: 'qtip-bootstrap',
					tip: {
						width: 16,
						height: 8
			    }
			  }
			});
	});

	// CLICKS
	cy.on('click','edge[role=\"DOM_interaction\"]', function(event) {
			var edge = event.target;
			var nodes = edge.connectedNodes();
			var doms = []
			nodes.forEach(function(node) {
				doms.push(node.data("label"));
			});

			var ds = edge.data("ds");
			// var links = edge.data("links").split("; ");
			// var all_links = []
			// links.forEach(function(link) {
			// 	all_links.push("<a href=\"https://thebiogrid.org/interaction/"+link+"\">"+link+"</a>");
			// });

			edge.qtip({
				content: "<b>"+ds+" interaction</b><br>"+
								 "<b>"+doms.join(" - ")+"</b><br>"+
								 "<a href=\"https://3did.irbbarcelona.org/dispatch.php?type=interaction&type1=domain&type2=domain&value1="+doms[0]+"&value2="+doms[1]+"\">link</a>",
								 // all_links.join(", "),
				position: {
           my: 'top center',
           at: 'bottom center'
       	},

				style: {
					classes: 'qtip-bootstrap',
					tip: {
						width: 16,
						height: 8
			    }
			  }
			});
	});


	cy.on('click','node[role=\"whole\"]', function(event) {
		var node = event.target;
		var gene = node.data("label");
		var des = node.data("des");
		var acc = node.data("protein");
		var length = node.data("length");
		node.qtip({
		  content: "<a style=\"color:#17A589;\">Gene</a> | " +
								"<b>"+gene+"</b><br>" +
								"<a style=\"color:#17A589;\">Protein</a> | " +
								"<b>"+des+"</b><br>" +
								"<a style=\"color:#17A589;\">UniProt</a> | " +
								"<a href=\"https://www.uniprot.org/uniprot/"+acc+"\">"+acc+"</a>"+
								"<br><a style=\"color:#17A589;\">Length</a> | "+length+"<br>",
		  position: {
		    my: 'top center',
		    at: 'bottom center'
		  },
		  style: {
				classes: 'qtip-bootstrap',
		    tip: {
		      width: 16,
		      height: 8
		    }
		  }
		});
	});

	cy.on('click','node[role=\"domain\"]', function(event) {
		var node = event.target;
		var name = node.data("label");
		var acc = node.data("acc");
		var des = node.data("des");
		var start = node.data("start");
		var end = node.data("end");
		node.qtip({
		  content: "<span style='color:#074987;'><b><i>"+name+"</i></b> "+
							 "(<a href=\"https://pfam.xfam.org/family/"+acc+"\">"+acc+"</a>)</span><br>"+
							 "<span style='background-color:#074987; color:white;'><b> "+des+" </b></span><br>"+
							 "<b><i>"+start+"-"+end+"</b></i><br>",

		  position: {
		    my: 'top center',
		    at: 'bottom center'
		  },
		  style: {
				classes: 'qtip-bootstrap',
		    tip: {
		      width: 16,
		      height: 8
		    }
		  }
		});
	});

  // BUTTON: Center the graph on the page
  $("#center").click(function(){
    // cy.center();
		cy.fit();
  });

  // BUTTON: Enables/Disables the option of clicking and draggin the graph nodes
  var lock_action = 0;
  $("#lock").click(function(){
    if ( lock_action == 0 ) {
      cy.autolock(true);
      //$("#lock i").removeClass("fa-unlock").addClass("fa-lock");
      $("#lock i").toggleClass("fa-unlock fa-lock");
      lock_action = 1;
    } else {
      cy.autolock(false);
      // $("#lock i").removeClass("fa-lock").addClass("fa-unlock");
      $("#lock i").toggleClass("fa-lock fa-unlock");
      lock_action = 0;
    }
  });

  // BUTTON: reset graphic
   $("#reset").click(function(){
     //all_pos = cy.position();
     //nodes.positions(all_pos);
     console.log("In reset function..");
     //cy.layout();
     //console.log( cy.nodes() );
     cy.nodes().forEach(function(n){
		 var p = n.data('orgPos');
		 n.position({ x: p.x, y: p.y });
     });
     cy.reset();
     cy.fit();
   });

	 // Manually adjust zoom
	 $("#zoom_in").click(function(){
		 var z = cy.zoom() + 0.2
		 cy.zoom( z )
	 });
	 $("#zoom_out").click(function(){
		 var z = cy.zoom() - 0.2
		 cy.zoom( z )
	 });


  // BUTTON: Toggle position labels
  //var coord_action = 0;
  $("#coord").click(function(){
    var eles = cy.$('node[role = "dom_pos"]');
    var checked = document.getElementById("coord").checked;

    if (checked) {
			eles.style("text-opacity", 1);
			eles.style("background-color", "#FFF");
      //coord_action = 1;
    } else {
      //coord_action = 0;
      eles.style("text-opacity", 0);
      eles.style("background-color", "#A9A9A9");

    }
  });

  // BUTTON: Toggles box surrounding the whole protein
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

  // BUTTON: Toggle regions:
  // Domains
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
  // ELMs
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
  // LMD2-LMs
  $("#toggle_newlms").click(function(){
    var nodes = cy.$('node[role="newLM"]');
		var edges = nodes.connectedEdges();
    var checked = document.getElementById("toggle_newlms").checked;
    if (checked) {
      nodes.style("visibility", "visible");
    } else {
			nodes.style("visibility", "hidden");
			edges.style("visibility", "hidden");
			$("#toggle_lmd2_int").prop("checked", false);
    }
  });
  // InterPreTS regions
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

  // BUTTON: Toggle interactions:
  // Protein-Protein
  //var pp_int_action = 0;
  $("#toggle_pp_int").click(function(){
    var edges = cy.$('edge[role="prot_prot_interaction"]');
    var checked = document.getElementById("toggle_pp_int").checked;
		if (checked) {
      edges.style("visibility", "visible");

    } else {
      edges.style("visibility", "hidden");
    }
  });

  // Dom-Dom (i)
  //var dom_int_action = 0;
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

  // Dom-Dom (ii)
  //var idom_int_action = 0;
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

  // ELM-Dom
  //var elm_int_action = 0;
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

  // LMD2-int
  //var lmd2_int_action = 0;
  $("#toggle_lmd2_int").click(function(){
    var edges = cy.$('edge[role="LMD2_interaction"]');
		var nodes = edges.connectedNodes();
    var checked = document.getElementById("toggle_lmd2_int").checked;
		if (checked) {
      edges.style("visibility", "visible");
			nodes.style("visibility", "visible");
			$("#toggle_newlms").prop("checked", true);
    } else {
      edges.style("visibility", "hidden");
    }
  });
  // InterpreTS
  //var prets_int_action = 0;
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

	// BUTTON: Shows/Hides phosphorylation sites
  $("#phospho").click(function(){
    var eles = cy.$('node[role="phosphorylation"]');
    var checked = document.getElementById("phospho").checked;
		if (checked) {
      eles.style("visibility", "visible");
    } else {
      eles.style("visibility", "hidden");
    }
  });
  // Button: Shows/Hides acetylation sites
  $("#acet").click(function(){
    var eles = cy.$('node[role="acetylation"]');
    var checked = document.getElementById("acet").checked;
    if (checked) {
      eles.style("visibility", "visible");
    } else {
      eles.style("visibility", "hidden");
    }
  });
  // Button: Shows/Hides mutations
  $("#mutation").click(function(){
    var eles = cy.$('node[role="mutation"]');
    var checked = document.getElementById("mutation").checked;
    if (checked) {
      eles.style("visibility", "visible");
    } else {
      eles.style("visibility", "hidden");
    }
  });


  // BUTTON: Download graph as...
  $("#dl_png").click(function(){
	  // console.log("in PNG function..");
		var image = cy.png()
		var iframe = "<iframe src='"+image+`' frameborder='0'
			style='border:0; top:0px; left:0px; bottom:0px; right:0px; width:100%; height:100%;'
			allowfullscreen></iframe>`;
  	var win = window.open();
	  win.document.write(iframe);

  });
    // BUTTON: Download graph as...
  $("#dl_jpg").click(function(){
	  // console.log("in JPG function..");
		var image = cy.jpg()
		var iframe = "<iframe src='"+image+`' frameborder='0'
			style='border:0; top:0px; left:0px; bottom:0px; right:0px; width:100%; height:100%;'
			allowfullscreen></iframe>`;
  	var win = window.open();
	  win.document.write(iframe);
  });
    // BUTTON: Download graph as...
  $("#dl_json").click(function(){
	  // console.log("in JSON function..");
	  var jsonBlob = new Blob([ JSON.stringify( cy.json() ) ], { type: 'application/javascript;charset=utf-8' });
	  saveAs( jsonBlob, 'graph.json' );
  });


});
