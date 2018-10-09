function prot_file_upload() {
	
	var fileInput = document.querySelector('#uploadProtFile');
	var textFromFile = "";
	var file = fileInput.files[0];
	var reader = new FileReader();
	reader.onload = function(e) {
		textFromFile = e.target.result;
		document.getElementById("query_prots").value = textFromFile;
	}
	reader.readAsText(file);	
};

function mut_file_upload() {
	
	var fileInput = document.querySelector('#uploadMutFile');
	var textFromFile = "";
	var file = fileInput.files[0];
	var reader = new FileReader();
	reader.onload = function(e) {
		textFromFile = e.target.result;
		document.getElementById("query_muts").value = textFromFile;
	}
	reader.readAsText(file);	
};
