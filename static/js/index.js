$(document).ready(function() {

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
			prots_input: "Please enter your proteins in the input box.\nIf you selected a local file, you need to upload it first, then submit.",
			max_prots: "The maximum number of interactors will be defaulted to 5.",
			// query_lmd2: "Please enter your LMD2 data."
		}
	});

	function resetInput(id) {
		document.getElementById(id).value="";
	}

	$("#all_intrs").click(function() {
		var input = document.getElementById('add_n_interactors');
		if (this.checked) {
			input.disabled = true;
		} else {
			input.disabled= false;
		}
	});

	// Examples
	$("#ex1").click(function () {
    $("#prots_input").val(`TCF3
CBFA2T3
ID3
PRKAR2A
HIF1A`);
	$("#muts_input").val(`CBFA2T3/D105G
CBFA2T3/N495D
CBFA2T3/A490P
TCF3/N554K
TCF3/N555D
TCF3/N555T
TCF3/V560L
TCF3/V560E
TCF3/D564V
TCF3/D564E
TCF3/E624G
ID3/C10STOP
ID3/M44V
ID3/C47STOP
ID3/Y48C
ID3/Y48STOP
ID3/L54M
ID3/P56S
ID3/P56A
ID3/P56L
ID3/P56H
ID3/L64F
ID3/L64V
ID3/L70V
ID3/Y76D
ID3/Q81STOP
ID3/Q81H
ID3/V82A
ID3/Q100STOP`);
	  $("#add_n_interactors").val(0);
	  $("#species").val("HUMAN");
  });

	$("#ex2").click(function () {
		$("#prots_input").val(`ARID1A
SIN3A
RB1
HDAC2
ACTL6A
`);
		$("#muts_input").val("");
		$("#add_n_interactors").val(0);
		$("#species").val("HUMAN");
	});

	$("#ex3").click(function () {
		$("#prots_input").val(`TP53
TIN1
`);
		$("#muts_input").val("");
		$("#add_n_interactors").val(2);
		$("#species").val("HUMAN");
//  		$('#query_lmd2').val(`
// #AC(a)	ID(a)	GN(a)	AC(b)	ID(b)	GN(b)	Type	F(a):F(b)	start	end
// O60504	VINEX_HUMAN	SORBS3	P08247	UniProt_ID	SYP	LMD2	PPxxPx[KR]	236	242
// O60504	VINEX_HUMAN	SORBS3	Q13153	UniProt_ID	PAK1	LMD2	PPxxPx[KR]	12	18
// O60504	VINEX_HUMAN	SORBS3	Q13153	UniProt_ID	PAK1	LMD2	ExEDD	175	179
// O60504	VINEX_HUMAN	SORBS3	Q13153	UniProt_ID	PAK1	LMD2	PPxxxxRP	187	194
// O60504	VINEX_HUMAN	SORBS3	P22681	UniProt_ID	CBL	LMD2	PPxxPx[KR]	544	550
// O60504	VINEX_HUMAN	SORBS3	P22681	UniProt_ID	CBL	LMD2	PPxPxxR	579	585
// O60504	VINEX_HUMAN	SORBS3	P22681	UniProt_ID	CBL	LMD2	PPxPxxR	823	829
// O60504	VINEX_HUMAN	SORBS3	Q8TDM6	UniProt_ID	DLG5	LMD2	PPxPxxR	1012	1018
// O60504	VINEX_HUMAN	SORBS3	Q8TDM6	UniProt_ID	DLG5	LMD2	PPxxPx[KR]	1166	1172
// O60504	VINEX_HUMAN	SORBS3	Q8TDM6	UniProt_ID	DLG5	LMD2	[KR]xPxxAP	1111	1117
// O60504	VINEX_HUMAN	SORBS3	Q8TDM6	UniProt_ID	DLG5	LMD2	[KR]P[ILMVF]PxxP	870	876
// O60504	VINEX_HUMAN	SORBS3	O00401	UniProt_ID	WASL	LMD2	PPxPxxR	278	284
// O60504	VINEX_HUMAN	SORBS3	O00401	UniProt_ID	WASL	LMD2	PPxPxxR	299	305
// O60504	VINEX_HUMAN	SORBS3	O00401	UniProt_ID	WASL	LMD2	PPxPxxR	310	316
// O60504	VINEX_HUMAN	SORBS3	O00401	UniProt_ID	WASL	LMD2	PPxPxxR	323	329
// O60504	VINEX_HUMAN	SORBS3	O00401	UniProt_ID	WASL	LMD2	PPxPxxR	335	341
// O60504	VINEX_HUMAN	SORBS3	O00401	UniProt_ID	WASL	LMD2	ExEDD	489	493`);
	});

	$("#ex4").click(function() {
		$("#prots_input").val(`#Drosophila's homeoproteins
EXD_DROME
UBX_DROME
TIN_DROME`);
		$("#muts").val(" ");
		$("#max_prots").val(0);
		$("#species").val("DROME");
	});

// 	$("#ex5").click(function() {
// 		$("#prots_input").val(`>BCR-ABL1
// MVDPVGFAEAWKAQFPDSEPPRMELRSVGDIEQELERCKASIRRLEQEVNQERFRMIYLQ
// TLLAKEKKSYDRQRWGFRRAAQAPDGASEPRASASRPQPAPADGADPPPAEEPEARPDGE
// GSPGKARPGTARRPGAAASGERDDRGPPASVAALRSNFERIRKGHGQPGADAEKPFYVNV
// EFHHERGLVKVNDKEVSDRISSLGSQAMQMERKKSQHGAGSSVGDASRPPYRGRSSESSC
// GVDGDYEDAELNPRFLKDNLIDANGGSRPPWPPLEYQPYQSIYVGGMMEGEGKGPLLRSQ
// STSEQEKRLTWPRRSYSPRSFEDCGGGYTPDCSSNENLTSSEEDFSSGQSSRVSPSPTTY
// RMFRDKSRSPSQNSQQSFDSSSPPTPQCHKRHRHCPVVVSEATIVGVRKTGQIWPNDGEG
// AFHGDADGSFGTPPGYGCAADRAEEQRRHQDGLPYIDDSPSSSPHLSSKGRGSRDALVSG
// ALESTKASELDLEKGLEMRKWVLSGILASEETYLSHLEALLLPMKPLKAAATTSQPVLTS
// QQIETIFFKVPELYEIHKEFYDGLFPRVQQWSHQQRVGDLFQKLASQLGVYRVLGYNHNG
// EWCEAQTKNGQGWVPSNYITPVNSLEKHSWYHGPVSRNAAEYLLSSGINGSFLVRESESS
// PGQRSISLRYEGRVYHYRINTASDGKLYVSSESRFNTLAELVHHHSTVADGLITTLHYPA
// PKRNKPTVYGVSPNYDKWEMERTDITMKHKLGGGQYGEVYEGVWKKYSLTVAVKTLKEDT
// MEVEEFLKEAAVMKEIKHPNLVQLLGVCTREPPFYIITEFMTYGNLLDYLRECNRQEVNA
// VVLLYMATQISSAMEYLEKKNFIHRDLAARNCLVGENHLVKVADFGLSRLMTGDTYTAHA
// GAKFPIKWTAPESLAYNKFSIKSDVWAFGVLLWEIATYGMSPYPGIDLSQVYELLEKDYR
// MERPEGCPEKVYELMRACWQWNPSDRPSFAEIHQAFETMFQESSISDEVEKELGKQGVRG
// AVSTLLQAPELPTKTRTSRRAAEHRDTTDVPEMPHSKGQGESDPLDHEPAVSPLLPRKER
// GPPEGGLNEDERLLPKDKKTNLFSALIKKKKKTAPTPPKRSSSFREMDGQPERRGAGEEE
// GRDISNGALAFTPLDTADPAKSPKPSNGAGVPNGALRESGGSGFRSPHLWKKSSTLTSSR
// LATGEEEGGGSSSKRFLRSCSASCVPHGAKDTEWRSVTLPRDLQSTGRQFDSSTFGGHKS
// EKPALPRKRAGENRSDQVTRGTVTPPPRLVKKNEEAADEVFKDIMESSPGSSPPNLTPKP
// LRRQVTVAPASGLPHKEEAGKGSALGTPAAAEPVTPTSKAGSGAPGGTSKGPAEESRVRR
// HKHSSESPGRDKGKLSRLKPAPPPPPAASAGKAGGKPSQSPSQEAAGEAVLGAKTKATSL
// VDAVNSDAAKPSQPGEGLKKPVLPATPKPQSAKPSGTPISPAPVPSTLPSASSALAGDQP
// SSTAFIPLISTRVSLRKTRQPPERIASGAITKGVVLDSTEALCLAISRNSEQMASHSAVL
// EAGKNLYTFCVSYVDSIQQMRNKFAFREAINKLENNLRELQICPATAGSGPAATQDFSKL
// LSSVKEISDIVQR

// BCR
// ABL1
// 		`);
// 		$("#muts_input").val("");
// 		$("#add_n_interactors").val(0);
// 		$("#species").val("HUMAN");
// 	});


});
