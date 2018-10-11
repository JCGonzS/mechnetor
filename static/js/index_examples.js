document.getElementById("ex1").addEventListener('click', function () {
      // var name = document.getElementById('name');
      // name.value = "Example 1";
      var text1 = document.getElementById('query_prots');
      text1.value = `TCF3
CBFA2T3
ID3`;
      var text2 = document.getElementById('query_muts');
      text2.value = `O75081/D105G CBFA2T3
P15923/N554K TCF3
P15923/N555D TCF3
P15923/N555T TCF3
P15923/V560L TCF3
P15923/V560E TCF3
P15923/D564V TCF3
P15923/D564E TCF3
P15923/E624G TCF3
Q02535/C10STOP ID3
Q02535/M44V ID3
Q02535/C47STOP ID3
Q02535/Y48C ID3
Q02535/Y48STOP ID3
Q02535/L54M ID3
Q02535/P56S ID3
Q02535/P56A ID3
Q02535/P56L ID3
Q02535/P56H ID3
Q02535/L64F ID3
Q02535/L64V ID3
Q02535/L70V ID3
Q02535/Y76D ID3
Q02535/Q81STOP ID3
Q02535/Q81H ID3
Q02535/V82A ID3
Q02535/Q100STOP ID3`;
  });

document.getElementById("ex2").addEventListener('click', function () {
      // var name = document.getElementById('name');
      // name.value = "Example 2"
      var text1 = document.getElementById('query_prots');
      text1.value = `SMARCA4
SMARCA2
WWOX
DPF2
SMARCB1
SMARCC1
SMARCC2
SMARCE1
ARID1B
ARID1A
SMARCD1
RB1
CBX5
ACTL6A
SIN3A
HDAC2
TP53
POLR2A
SMARCD3
TAT`;
      var text2 = document.getElementById('query_muts');
      text2.value = `P51532/T308M SMARCA4
P51532/V1404I SMARCA4
P51532/V1552I SMARCA4`;
});

document.getElementById("ex3").addEventListener('click', function () {
      // var name = document.getElementById('name');
      // name.value = "Example 3 with LMD2 data"
      var text1 = document.getElementById('query_prots');
      text1.value = `SORBS3
PAK1
CBL
SYP
DLG5
WASL`;
//       var text3 = document.getElementById('query_lmd2');
//       text3.value = `#Example
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
// O60504	VINEX_HUMAN	SORBS3	O00401	UniProt_ID	WASL	LMD2	ExEDD	489	493`;
});
