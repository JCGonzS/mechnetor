import sys
import os
import json
from flask import Flask, render_template, request, url_for, redirect, jsonify

from mechnetor_app import app
from mechnetor_app.lib.utils import print_log, get_unique_random_identifier
from mechnetor_app.config import get_paths, get_templates
from mechnetor_app import mechnetor


path = get_paths()
template = get_templates()
sys.stdout = open(path["log"], 'a')
print_log("test", "testing")


def get_stats_for_charts(stats_file):
    with open(stats_file, "r") as f:
        d = json.load(f)
    return d["table_columns"], d["not_found"], d["no_int_prots"]


@app.route("/")
@app.route("/index")
def load_index():
    maintenance = False
    if maintenance:
        return render_template(template["maintenance"])
    else:
        return render_template(template["index"])


@app.route("/help")
def load_help():
    return render_template(template["help"])


@app.route("/output", methods=["POST", "GET"])
def load_output():
    if request.method == "POST":
        input_d = {}
        # input_d["query_name"] = request.form["query_name"]
        input_d["prots_input"] = request.form["prots_input"]
        input_d["muts_input"] = request.form["muts_input"]
        input_d["sps"] = request.form["species"]
        input_d["make_graph"] = True

        input_d["hide"] = False
        if request.form.get("hide_no_int"):
            input_d["hide"] = True

        input_d["add_interactors"] = 0
        if request.form.get("all_interactors"):
            input_d["add_interactors"] = "all"
        elif request.form.get("add_n_interactors"):
            input_d["add_interactors"] = int(request.form["add_n_interactors"])

        input_d["only"] = False
        if request.form.get("only_int_pairs"):
            input_d["only"] = True   

        job_id = get_unique_random_identifier(path["jobs"])
        job_dir = path["jobs"]+"job_"+job_id+"/"
        print_log(job_id, "New Job {}".format(job_id))
        try:
            os.mkdir(job_dir)
        except OSError:
            print_log(job_id,
                      "Creation of the directory {} failed".format(job_dir))

        with open(job_dir+"input_"+job_id+".json", "w") as out:
            json.dump(input_d, out)

        return redirect(url_for("run_job", job_id=job_id))


@app.route('/job/<job_id>', methods=['GET', 'POST'])
def run_job(job_id):
    job_dir = path["jobs"]+"job_"+job_id+"/"
    graph_json = "graph_elements_"+job_id+".json"
    table_file = "interaction_table_"+job_id+".json"
    stats_file = "req_parameters_"+job_id+".json"

    if not (os.path.isfile(job_dir+graph_json) and
            os.path.isfile(job_dir+table_file) and
            os.path.isfile(job_dir+stats_file)):

        with open(job_dir+"input_"+job_id+".json", "rt") as f:
            d = json.load(f)

        if d["sps"]=="HUMAN":
            db = "piv_human"
        else:
            db = "piv_all"
        # try:
        error = mechnetor.main(
                        INPUT_1=d["prots_input"],
                        INPUT_2=d["muts_input"],
                        ORG=d["sps"],
                        ADDITIONAL_INTERACTORS=d["add_interactors"],
                        ONLY_INT_PAIRS=d["only"],
                        MAIN_OUTPUT_DIR=path["jobs"],
                        CUSTOM_ID=job_id,
                        DATA_DIR=path["data"],
                        BLASTDB_DIR="/net/home.isilon/ds-russell/blastdb/",
                        PSQL_USER="bq_jgonzalez", PSQL_DB=db,
                        MAKE_NETWORK=d["make_graph"],
                        HIDE_NO_INT=d["hide"],
                        TABLE_FORMAT="json"
                        )
        # except:
            # print "except error"
        #     return render_template(error_template)
        if error:
            if error == 1:
                return render_template(template["error1"])
            elif error == 2:
                return render_template(template["error2"])
    else:
        print_log(job_id, "Files exist")

    # Read stats file
    (table_columns, not_found, no_int_prots) = get_stats_for_charts(job_dir+stats_file)

    return render_template(
        template["results"],
        graph_json="jobs/"+"job_"+job_id+"/"+graph_json,
        ints_json="jobs/"+"job_"+job_id+"/"+table_file,
        stats_json="jobs/"+"job_"+job_id+"/"+stats_file,
        not_found=not_found,
        no_int_prots=no_int_prots,
        table_columns=table_columns
    )

