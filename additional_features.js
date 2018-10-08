"use strict";

$(document).ready(function(){

	//populate datatable
	var table = $('#interactionTable').DataTable( {
		dom: "lBfrtip",
		buttons: [
            'copy', 'csv', 'excel', 'pdf', 'print'
        ],
		"ajax": ints_json,
		"columnDefs": [
			  //adds links to columns with UniProt accession
				{
					targets: [1, 3],
					"render": function ( data, type, row, meta ) {
						return '<a href=http://www.uniprot.org/uniprot/'+data+'  target="_blank">'+data+'</a>';
						}
				},
				{
					targets: [0, 1, 2, 3],
					"width": "120px"
				},
				//gives text wrap width
				{
					targets: [-2],
					"render": function (data, type, row, meta) {
						return "<div class='text-wrap width-300'>" + data + "</div>";
						}
				},
				//hides column
				{
					targets: [-2],
					"visible": false,
					"searchable": false
				}
			],
		"initComplete": function(settings){
			$('#interactionTable thead th, #interactionTable tfoot th').each(function () {
			   var td = $(this);
			   var headerText = td.text();
			   var headerTitle= "";
			   if ( headerText == "Protein (A)" ) headerTitle =  "Gene Name of Interacting Protein A";
				 else if ( headerText == "Accession (A)" ) headerTitle =  "UniProt Accession of Interacting Protein A";
				 else if ( headerText == "Protein (B)" ) headerTitle =  "Gene Name of Interacting Protein B";
				 else if ( headerText == "Accession (A)" ) headerTitle =  "UniProt Accession of Interacting Protein A";
			   else if ( headerText == "Type" ) headerTitle =  "Type of Interaction";
			   else if ( headerText == "F (A)" ) headerTitle =  "Interacting Element of Protein A";
				 else if ( headerText == "F (A)" ) headerTitle =  "Interacting Element of Protein B";
				 else if ( headerText == "Source" ) headerTitle =  "Source of interaction information (DB, tool, etc)";
			   td.attr('title', headerTitle);
			});

      /* Apply the tooltips */
      $('#example thead th[title]').tooltip(
	      {
	         "container": 'body',
	         "delay": 0,
					 "track": true,
					 "fade": 250
      });

      //show only when hovering within Type column
			$('#interactionTable tr td:nth-child(3)').each(function () {
         var td = $(this);
         var headerText = td.text();
         var headerTitle= "";
         //console.log( td.text() );
         if ( headerText == "PROT::PROT" ) headerTitle =  "Protein-Protein interaction; supported by experimental evidence from BioGRID";
         else if (headerText == "DOM::DOM" ) headerTitle = "Domain-Domain interaction; derived via HiRes 3D structure of interaction domains from 3did";
         else if (headerText == "iDOM::iDOM" ) headerTitle = "Domain-Domain interaction; predicted from protein signatures";
         else if (headerText == "DOM::ELM" ) headerTitle = "Domain-Linear Motif interaction; from ELM database";
				 else if (headerText == "ELM::DOM" ) headerTitle = "Domain-Linear Motif interaction; from ELM database";
         else if (headerText == "InterPreTS" ) headerTitle = "Region-Region interaction; based on homology with a structure of interacting proteins, predicted with InterPreTS";

         td.attr('title', headerTitle);
      });

      /* Apply the tooltips */
      $('#example thead th[title]').tooltip(
      {
         "container": 'body',
         "delay": 0,
				 "track": true,
				 "fade": 250
      });
    }
	});


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


  // Lock PTM & Mutation nodes
  cy.nodes().nonorphans()
    .on('grab', function(){
      cy.$('node[role = "pp_mod"]').ungrabify();
      cy.$('node[role = "ac_mod"]').ungrabify();
      cy.$('node[role = "mutation"]').ungrabify();
      cy.$('node[role = "position"]').ungrabify();
    })
    .on('free', function(){
      cy.$('node[role = "pp_mod"]').grabify();
      cy.$('node[role = "ac_mod"]').grabify();
      cy.$('node[role = "mutation"]').grabify();
      cy.$('node[role = "position"]').grabify();
    });

  // Mouseover nodes display info
  cy.on('mouseover','node[role = \"whole\"]', function(event) {
    var node = event.cyTarget;
    node.qtip({
      content: 'hello',
      show: {
        event: event.type,
        ready: true
      },
      hide: {
        event: 'mouseout unfocus'
      }
    }, event);
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

  // BUTTON: Remove/Restore selected elements from the graph
  // $("#remove").click(function(){
  //   var rem = cy.$(':selected').remove();
  // });
  //
  // $("#restore").click(function(){
  //   rem.restore();
  // })

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
    var eles = cy.$('node[role = "whole"]');
    var eles2 = cy.$('node[role = "whole"]:selected');
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
    var eles = cy.$('node[role = "domain"]');
    var checked = document.getElementById("toggle_doms").checked;
    if (checked) {
			eles.style("visibility", "visible");
			// eles.style("display","element");
    } else {
			eles.style("visibility", "hidden");
			// eles.style("display","none");
    }
  });
  // ELMs
  $("#toggle_elms").click(function(){
    var eles = cy.$('node[role = "elms"]');
    var checked = document.getElementById("toggle_elms").checked;
    if (checked) {
      eles.style("visibility", "visible");
    } else {
      eles.style("visibility", "hidden");
    }
  });
  // LMD2-LMs
  $("#toggle_newlms").click(function(){
    var eles = cy.$('node[role = "newLM"]');
    var checked = document.getElementById("toggle_newlms").checked;
    if (checked) {
      eles.style("visibility", "visible");
    } else {
			eles.style("visibility", "hidden");
    }
  });
  // InterPreTS regions
  $("#toggle_iprets").click(function(){
    var eles = cy.$('node[role = "iprets"]');
    var checked = document.getElementById("toggle_iprets").checked;
    if (checked) {
      eles.style("visibility", "visible");
    } else {
			eles.style("visibility", "hidden");
    }
  });

  // BUTTON: Toggle interactions:
  // Protein-Protein
  //var pp_int_action = 0;
  $("#toggle_pp_int").click(function(){
    var eles = cy.$('edge[role = "prot_prot_interaction"]');
    var checked = document.getElementById("toggle_pp_int").checked;
		if (checked) {
      eles.style("visibility", "visible");
    } else {
      eles.style("visibility", "hidden");
    }
  });
  // Dom-Dom (i)
  //var dom_int_action = 0;
  $("#toggle_dom_int").click(function(){
    var eles = cy.$('edge[role = "DOM_interaction"]');
    var checked = document.getElementById("toggle_dom_int").checked;
		if (checked) {
      eles.style("visibility", "visible");
    } else {
      eles.style("visibility", "hidden");
    }
  });
  // Dom-Dom (ii)
  //var idom_int_action = 0;
  $("#toggle_idom_int").click(function(){
    var eles = cy.$('edge[role = "iDOM_interaction"]');
    var checked = document.getElementById("toggle_idom_int").checked;
		if (checked) {
      eles.style("visibility", "visible");
    } else {
      eles.style("visibility", "hidden");
    }
  });
  // ELM-Dom
  //var elm_int_action = 0;
  $("#toggle_elmdom_int").click(function(){
    var eles = cy.$('edge[role = "ELM_interaction"]');
    var checked = document.getElementById("toggle_elmdom_int").checked;
		if (checked) {
      eles.style("visibility", "visible");
    } else {
      eles.style("visibility", "hidden");
    }
  });
  // LMD2-int
  //var lmd2_int_action = 0;
  $("#toggle_lmd2_int").click(function(){
    var eles = cy.$('edge[role = "LMD2_interaction"]');
    var checked = document.getElementById("toggle_lmd2_int").checked;
		if (checked) {
      eles.style("visibility", "visible");
    } else {
      eles.style("visibility", "hidden");
    }
  });
  // InterpreTS
  //var prets_int_action = 0;
  $("#toggle_prets_int").click(function(){
    var eles = cy.$('edge[role = "INT_interaction"]');
    var checked = document.getElementById("toggle_prets_int").checked;
		if (checked) {
      eles.style("visibility", "visible");
    } else {
      eles.style("visibility", "hidden");
    }
  });

	// BUTTON: Shows/Hides phosphorylation sites
  $("#phospho").click(function(){
    var eles = cy.$('node[role = "pp_mod"]');
    var checked = document.getElementById("phospho").checked;
		if (checked) {
      eles.style("visibility", "visible");
    } else {
      eles.style("visibility", "hidden");
    }
  });
  // Button: Shows/Hides acetylation sites
  $("#acet").click(function(){
    var eles = cy.$('node[role = "ac_mod"]');
    var checked = document.getElementById("acet").checked;
    if (checked) {
      eles.style("visibility", "visible");
    } else {
      eles.style("visibility", "hidden");
    }
  });
  // Button: Shows/Hides mutations
  $("#mutation").click(function(){
    var eles = cy.$('node[role = "mutation"]');
    var checked = document.getElementById("mutation").checked;
    if (checked) {
      eles.style("visibility", "visible");
    } else {
      eles.style("visibility", "hidden");
    }
  });


  // BUTTON: Download graph as...
  $("#dl_png").click(function(){
	  console.log("in PNG function..");
		var image = cy.png()
		var iframe = "<iframe src='"+image+`' frameborder='0'
		style='border:0; top:0px; left:0px; bottom:0px; right:0px; width:100%; height:100%;'
		allowfullscreen></iframe>`;
  	var win = window.open();
	  win.document.write(iframe);

  });
    // BUTTON: Download graph as...
  $("#dl_jpg").click(function(){
	  console.log("in JPG function..");
	  // var b64key = 'base64,';
	  // var b64 = cy.jpg().substring( content.indexOf(b64key) + b64key.length );
		// var b64toBlob = require('b64-to-blob')
	  // var imgBlob = base64ToBlob( b64, 'image/png' );
		// var imgBlob = saveBase64ImagePng( b64,'graph.png')
	  // saveAs( imgBlob, 'graph.png' );
		var image = cy.jpg()
		var iframe = "<iframe src='"+image+`' frameborder='0'
			style='border:0; top:0px; left:0px; bottom:0px; right:0px; width:100%; height:100%;'
			allowfullscreen></iframe>`;
  	var win = window.open();
	  win.document.write(iframe);
  });
    // BUTTON: Download graph as...
  $("#dl_json").click(function(){
	  console.log("in JSON function..");
	  var jsonBlob = new Blob([ JSON.stringify( cy.json() ) ], { type: 'application/javascript;charset=utf-8' });
	  saveAs( jsonBlob, 'graph.json' );
  });


});
