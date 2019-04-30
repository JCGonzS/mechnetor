$(document).ready(function() {

	//~ Custom rule to parse protein input

	//~ jQuery.validator.addMethod("domain", function(value, element) {
		//~ return this.optional(element) || /^http:\/\/mycorporatedomain.com/.test(value);
		//~ }, "Please specify the correct domain for your documents");

	// validate input form on keyup and submit
	$("#form").validate({
		rules: {
			prots_input: "required",
			// max_prots: {
			// 	required: "#show_max_prots:checked",
			// 	digits: true
			// },
			// query_lmd2: {
			// 	required: "#show_query_lmd2:checked"
			// }
		},
		messages: {
			prots_input: "Please enter your proteins in the input box.\nIf you selected a local file, you need to upload it first, the submit.",
			max_prots: "The maximum number of interactors will be defaulted to 5.",
			// query_lmd2: "Please enter your LMD2 data."
		}
	});

  $("#all_prots").click(function() {
		var input = document.getElementById('max_prots');
		if (this.checked) {
			input.disabled = true;
		} else {
			input.disabled= false;
		}
	});
	//
	// $("#show_query_lmd2").click(function() {
	// 	var input = document.getElementById('query_lmd2');
	// 	if (this.checked) {
	// 		input.hidden = false;
	// 		input.focus();
	// 		//input.value="see Example 3 format";
	// 	} else {
	// 		input.hidden = true;
	// 	}
	// });

});
