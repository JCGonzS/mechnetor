$(document).ready(function() {
	
	//~ Custom rule to parse protein input
	
	//~ jQuery.validator.addMethod("domain", function(value, element) {
		//~ return this.optional(element) || /^http:\/\/mycorporatedomain.com/.test(value);
		//~ }, "Please specify the correct domain for your documents");

	
	
	// validate input form on keyup and submit
	$("#form").validate({
		rules: {
			query_prots: "required",
			max_prots: {
				required: "#show_max_prots:checked",
				digits: true
			},
			query_lmd2: {
				required: "#show_query_lmd2:checked"
			}
			
		},
		messages: {
			query_prots: "Please enter your proteins interactions.",
			max_prots: "Please enter a valid number to limit the number of proteins to.",
			query_lmd2: "Please enter your LMD2 data."
		}
	});



    $("#show_max_prots").click(function() {
		var input = document.getElementById('max_prots');
		if (this.checked) {
			input.hidden = false;
			input.focus();
			input.value=10;
		} else {
			input.hidden=true;
		}
	});
	
	$("#show_query_lmd2").click(function() {
		var input = document.getElementById('query_lmd2');
		if (this.checked) {
			input.hidden = false;
			input.focus();
			//input.value="see Example 3 format";
			
		} else {
			input.hidden = true;
		}
	});


	

});
