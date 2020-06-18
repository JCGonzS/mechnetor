"use strict";

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
cy.nodes().on("tap", function( e ){
    var node = e.target;
    var role = node.data("role");
    if( role == "protein_main" ){
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
                "<a href='https://www.uniprot.org/uniprot/"+acc+"'>" +
                acc+" <i class='fas fa-external-link-alt fa-xs'></i>" +
                "</a><br>" +
                "<span class='tipProt'>Length</span> | "+
                "<a href='https://www.uniprot.org/uniprot/"+acc+".fasta'>" +
                length+" AA <i class='fas fa-external-link-alt fa-xs'></i>"+
                "</span>",
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
    }
});

cy.nodes().on("tap", function( e ){
    var node = e.target;
    var role = node.data("role");
    if( role == "user_seq" ){
        var label = node.data("label");
        var length = node.data("length");
        var link = node.data("link");
        var blast = node.data("blast");
        var blast_hits = "<span class='tipProt'>Blast hits:</span><br>\n";
        blast_hits += "<div style='margin: 0 auto;'>";
        blast.forEach(function(hit) {
            var coor = hit.split("|")[0]
            var acc = hit.split("|")[1]
            var gene = hit.split("|")[2]
            var info = hit.split("|")[3]
            blast_hits+="\t<span style='float:left;width:80px;'><b>"+coor+"</b></span>"
            blast_hits+=" | <a href='https://www.uniprot.org/uniprot/"+acc+"'>"+acc+"</a>"
            blast_hits+=" <b>"+gene+"</b> "+info+"<br>\n"
        });
        blast_hits += "</div>"
        node.qtip({
            content:
                "<span class='tip'>User-input sequence<br>\n"+
                "\t<span class='tipProt'>Label</span> | <b>"+label+"</b><br>\n" +
                "\t<span class='tipProt'>Length</span> | "+
                "<a href='"+link+"' target='_blank'_>" +
                length+" AA <i class='fas fa-external-link-alt fa-xs'></i></a><br>"+
                blast_hits+
                "</span>",
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
    }
});

cy.on("click", "node[role='domain']", function(event) {
    var node = event.target;
    var name = node.data("label");
    var acc = node.data("acc");
    var des = node.data("des");
    var start = node.data("start");
    var end = node.data("end");
    var prot = node.data("protein");
    var e_val = Number.parseFloat(node.data("e_val")).toExponential(2);
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
                    " (<a href='https://pfam.xfam.org/protein/"+prot+"'>"+prot+" <i class='fas fa-external-link-alt fa-xs'></i></a>)<br>"+
                "<span class='tipPfam'>E-value</span> | "+e_val+
            "</span>",
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
    var label = node.data("label");
    var acc 	= node.data("acc");
    var name 	= node.data("name");
    var regex = node.data("regex");
    var start = node.data("start");
    var end 	= node.data("end");
    var seq 	= node.data("seq");
    var prot 	= node.data("protein");
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

    node.qtip({
        content:
            "<span class='tip' style='color: #7f7c7b;'>\n" +
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
            "</div>",
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

    node.qtip({
        content:
            "<span class='tip' style='color: #7b241c;'>" +
            "<b>Predicted with <a href='http://www.russelllab.org/cgi-bin/tools/interprets.pl/interprets.pl'>InterPreTS</a></b>" +
            "</span><br>" +
            "<span class='tip'>"+
                "<span class='tipInP'>PDB Template ID</span> | " +
                "<a href='https://www.rcsb.org/structure/"+pdb+"'><b>"+pdb+"</b> <i class='fas fa-external-link-alt fa-xs'></i></a> chain "+chain+"<br>"+
                "<span class='tipInP'>Alignment</span> |  "+
                "Protein/<b>"+start+"-"+end+"</b>; Template/"+pdb_start+"-"+pdb_end+"<br>"+
                "<span class='tipInP'>Alignment score</span> | BLAST e-val=<i>"+ev+"</i>, <i>"+pcid+"%</i> id<br>"+
            "</span>",
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

    node.qtip({
        content:
            "<span class='tip' style='color: #28446f;'>" +
                "<b><a href='https://cancer.sanger.ac.uk/cosmic'>COSMIC <i class='fas fa-external-link-alt fa-xs'></i></a></b><br>" +
                "<a href="+gene_link+">View gene "+gene+" <i class='fas fa-external-link-alt fa-xs'></i></a>"+
            "</span><br>\n"+html,
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

        edge.qtip({
            content:
                "<span class='tip' style='color: #16a085;'>" +
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
                "</span>",

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
        var ds 		= edge.data("ds");
        var lo 		= Number(edge.data("lo")).toFixed(2);
        var p_val = Number.parseFloat(edge.data("p_val")).toExponential(2);

        edge.qtip({
            content:
                "<span class='tip' style='color: #d4ac0d;'>" +
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
                "</span>",

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

        edge.qtip({
            content:
                "<span class='tip' style='color: #b95db9;'>" +
                    "Source: <a href='http://elm.eu.org/downloads.html#interactions'><b>ELM</b> (interactions)</a></b>" +
                "</span><br>" +
                "<span class='tip'>"+
                    "<span class='tipELMint'>Interacting Elements</span> | "+
                        doms.join(" - ")+
                    "</span><br>"+
                "<span class='tipELMint'>Interacting Proteins</span> | "+prots.join(" - ")+
                "</span><br>"+
                "<span class='tipELMint'>Interaction P-value</span> | "+p_val+
                "</span>",
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
    var pdb = edge.data("pdb").split("|")[1];
    var z = edge.data("z-score");
    var pval = edge.data("p-value");
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
                "<a href='https://www.rcsb.org/structure/"+pdb+"'><b>"+pdb+"</b> <i class='fas fa-external-link-alt fa-xs'></i></a><br>"+
                "<span class='tipInP'>PDB Template Chains</span> | "+chains.join(" - ")+"<br>"+
                "<span class='tipInP'>Interacting Proteins</span> | "+prots.join(" - ")+"<br>"+
                "<span class='tipInP'>Z-score</span> | "+z+"<br>"+
                "<span class='tipInP'>P-value</span> | "+pval+
            "</span>",
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
