"use strict";

$(document).ready(function(){
	$('#interaction-table').DataTable( {
		ajax: ints_json
	});
// 	//populate datatable
// 	var table = $('#interaction-table').DataTable( {
// 		dom: "lBfrtip",
// 		buttons: ['copy', 'csv', 'excel', 'pdf', 'print'],
// 		"ajax": ints_json,
// 		"scrollX": true,
// 		"columnDefs": [
// 			//adds links to columns with UniProt accession
// 			{
// 				targets: [1, 3],
// 				"render": function ( data, type, row, meta ) {
// 					return '<a href=http://www.uniprot.org/uniprot/'+data+'  target="_blank">'+data+'</a>';
// 					}
// 			},
// 			// {
// 			// 	targets: [0, 1, 2, 3],
// 			// 	"width": "120px"
// 			// },
// 			//gives text wrap width
// 			// {
// 			// 	targets: [-2],
// 			// 	"render": function (data, type, row, meta) {
// 			// 		return "<div class='text-wrap width-300'>" + data + "</div>";
// 			// 		}
// 			// },
// 			//hides column
// 			// {
// 			// 	targets: [-2],
// 			// 	"visible": false,
// 			// 	"searchable": false
// 			// }
// 		],
// 		"initComplete": function(settings){
// 			$('#interaction-table thead th, #interaction-table tfoot th').each(function () {
// 				var td = $(this);
// 			   	var headerText = td.text();
// 				var headerTitle= "";
// 			   	if ( headerText == "Protein (A)" ) headerTitle =  "Gene Name of Interacting Protein A";
// 				else if ( headerText == "Accession (A)" ) headerTitle =  "UniProt Accession of Interacting Protein A";
// 				else if ( headerText == "Protein (B)" ) headerTitle =  "Gene Name of Interacting Protein B";
// 				else if ( headerText == "Accession (B)" ) headerTitle =  "UniProt Accession of Interacting Protein A";
// 			   	else if ( headerText == "InteractionType" ) headerTitle =  "Type of Interaction";
// 			   	else if ( headerText == "F (A)" ) headerTitle =  "Interacting Element of Protein A";
// 				else if ( headerText == "Start-End (A)" ) headerTitle =  "Start and End Positions of Element of Protein A";
// 				else if ( headerText == "Mutations (A)" ) headerTitle =  "User Input-Mutations within Element of Protein A";
// 				else if ( headerText == "F (B)" ) headerTitle =  "Interacting Element of Protein B";
// 				else if ( headerText == "Start-End (B)" ) headerTitle =  "Start and End Positions of Element of Protein B";
// 				else if ( headerText == "Mutations (B)" ) headerTitle =  "User Input-Mutations within Element of Protein B";
// 				else if ( headerText == "DataSource" ) headerTitle =  "Source of interaction information (DB, tool, etc)";
// 			   	td.attr('title', headerTitle);
// 			});

			
// 			//show only when hovering within Type column
// 			$('#interaction-table tr td:nth-child(3)').each(function () {
// 				var td = $(this);
// 				var headerText = td.text();
// 				var headerTitle= "";
// 				if ( headerText == "PROT::PROT" ) headerTitle =  "Protein-Protein interaction; supported by experimental evidence from BioGRID";
// 				else if (headerText == "DOM::DOM" ) headerTitle = "Domain-Domain interaction; derived via HiRes 3D structure of interaction domains from 3did";
// 				else if (headerText == "iDOM::iDOM" ) headerTitle = "Domain-Domain interaction; predicted from protein signatures";
// 				else if (headerText == "DOM::ELM" ) headerTitle = "Domain-Linear Motif interaction; from ELM database";
// 				else if (headerText == "ELM::DOM" ) headerTitle = "Domain-Linear Motif interaction; from ELM database";
// 				else if (headerText == "InterPreTS" ) headerTitle = "Region-Region interaction; based on homology with a structure of interacting proteins, predicted with InterPreTS";
// 				td.attr('title', headerTitle);
// 			});
			
// 			/* Apply the tooltips */
// 			$('#interaction-table thead th[title]').tooltip({
// 				"container": 'body',
// 				"delay": 0,
// 				"track": true,
// 				"fade": 250
// 			});

// 	    }
// 	});
});
