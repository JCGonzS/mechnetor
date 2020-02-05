import sys, os, time, string, random, datetime, json
from jc_app import app
from flask import Flask, render_template, request, url_for, redirect, jsonify
from pymongo import MongoClient
from flask_debugtoolbar import DebugToolbarExtension
from flask_debugtoolbar_lineprofilerpanel.profile import line_profile
import pipeline

main_dir = "/var/www/flask_apps/jc_test/jc_app/"
output_dir = main_dir+"static/jobs/"
log_file = main_dir+"static/jobs/log.txt"
sys.stdout = open(log_file, 'a')

# Jinja Templates
index_template = "index.html.jinja2"
help_template = "features.html.jinja2"
error_template = "input_error.html.jinja2"
results_template = "results.html.jinja2"

client = MongoClient('localhost', 27017)

def print_log(job_id, msg):
    st = "[{}]".format(datetime.datetime.now())+" [JOB ID: "+job_id+"]"
    print st, msg

def get_unique_random_identifier(output_dir):
    flag = 0
    while flag == 0:
        ide = ''.join(random.choice(
                string.ascii_uppercase + string.ascii_lowercase + string.digits
                ) for _ in range(8))
        if not os.path.exists(output_dir+"job_"+ide+"/"):
            flag = 1
    return ide

def get_stats_for_charts(prot_ints_number, int_types_number):
    prots_labels, prots_counts = [], []
    for k in sorted(prot_ints_number, key=prot_ints_number.get):
        prots_labels.append(str(k))
        prots_counts.append(prot_ints_number[k])

    colors = {
        "PROT::PROT": "#5F6A6A",
        "DOM::DOM": "#16A085",
        "iDOM::iDOM": "#D4AC0D",
        "ELM::DOM": "#AF7AC5"
    }
    names = {
        "PROT::PROT": "Binary",
        "DOM::DOM": "Domain-Domain",
        "iDOM::iDOM": "Domain-Domain (inferred)",
        "ELM::DOM": "Linear motifs - domain"
    }
    int_types_labels, int_types_series = [], []
    for k in sorted(int_types_number, key=int_types_number.get):
        if k == "iELM::DOM":
            continue
        int_types_labels.append(str(k))
        int_types_series.append({
                                    "value": int_types_number[k],
                                    "itemStyle": {"color": colors[k]}
                                })

    return prots_labels, prots_counts, int_types_labels, int_types_series

@app.route("/")
@app.route("/index")
def index():
    return render_template(index_template)

@app.route("/help")
def help():
	return render_template(help_template)

@line_profile
@app.route("/output", methods = ["POST", "GET"])
def output():
    if request.method == "POST":
        input_d = {}
        #query_name = request.form["query_name"]
        input_d["prots_input"] = request.form["prots_input"]
        input_d["muts_input"] = request.form["muts_input"]
        input_d["sps"] = request.form["species"]
        input_d["add_interactors"] = "all"
        if "add_intrs" in request.form:
            input_d["add_interactors"] = int(request.form["add_intrs"])
        input_d["hide"] = False
        if request.form.get("hide_no_int"):
            input_d["hide"] = True
        input_d["make_graph"] = True

        job_id = get_unique_random_identifier(output_dir)
        job_dir = output_dir+"job_"+job_id+"/"
        try:
            os.mkdir(job_dir)
        except OSError:
            print_log(job_id,
                      "Creation of the directory {} failed".format(job_dir))

        with open(job_dir+"input_"+job_id+".json", "w") as out:
            json.dump(input_d, out)

        return redirect(url_for("run_job", job_id=job_id))# prots_input=prots_input))

@app.route('/job/<job_id>', methods=['GET', 'POST'])
def run_job(job_id):
    job_dir = output_dir+"job_"+job_id+"/"
    graph_json = "graph_elements_"+job_id+".json"
    table_file = "interaction_table_"+job_id+".json"
    first_run = False

    if not (os.path.isfile(job_dir+graph_json) and os.path.isfile(job_dir+table_file)):

        with open(job_dir+"input_"+job_id+".json", "r") as f:
            d = json.load(f)

        # try:
        results = pipeline.main(
                        INPUT_1=d["prots_input"],
                        INPUT_2=d["muts_input"],
                        SP=d["sps"],
                        ADDITIONAL_INTERACTORS=d["add_interactors"],
                        MAIN_OUTPUT_DIR=output_dir,
                        CUSTOM_ID=job_id,
                        BLASTDB_DIR="/net/home.isilon/ds-russell/blastdb/",
                        CLIENT=client,
                        MAKE_NETWORK=d["make_graph"],
                        HIDE_NO_INT=d["hide"],
                        TABLE_FORMAT="json"
                        )
        first_run = True
        # except:
            # print "except error"
        #     return render_template(error_template)

        if results == "Error1":
            return render_template(error_template)
        else:
            param = {} # parameters
            (param["not_found"], param["no_int_prots"],
                param["prot_ints_number"], param["int_types_number"]) = results

    # Save required info in extra file
    if first_run:
        with open(job_dir+"req_parameters_"+job_id+".json", "w") as out:
            json.dump(param, out)
    else:
        with open(job_dir+"req_parameters_"+job_id+".json", "r") as f:
            param = json.load(f)

    (prots_labels, prots_counts,
    int_types_labels, int_types_series) = get_stats_for_charts(
                                                     param["prot_ints_number"],
                                                     param["int_types_number"])

    return render_template(results_template,
                       graph_json="jobs/"+"job_"+job_id+"/"+graph_json,
                       ints_json="jobs/"+"job_"+job_id+"/"+table_file,
                       not_identified=param["not_found"],
                       no_int_prots=param["no_int_prots"],
                       chart1_labels=prots_labels,
                       chart1_series=prots_counts,
                       chart2_labels=int_types_labels,
                       chart2_series=int_types_series
                       )
