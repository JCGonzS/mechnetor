import sys, os, time, string, random, datetime, json
from flask import Flask, render_template, request, url_for, redirect, jsonify
from flask_debugtoolbar import DebugToolbarExtension
from flask_debugtoolbar_lineprofilerpanel.profile import line_profile
from importlib.machinery import SourceFileLoader
from piv2_app import app 
from . import mechnetor

main_dir = "/var/www/flask_apps/piv2/piv2_app/"
data_dir = main_dir+"static/data/"
output_dir = main_dir+"static/jobs/"
log_file = main_dir+"static/jobs/log.txt"
sys.stdout = open(log_file, 'a')

# Jinja Templates
index_template = "index.html"
maintenance_template = "maintenance_index.html"
help_template = "features.html"
error1_template = "input_error.html"
error2_template = "toobig_error.html"
results_template = "results.html"

def print_log(job_id, msg):
    st = "[{}]".format(datetime.datetime.now())+" [JOB ID: "+job_id+"]"
    print(st, msg)

def get_unique_random_identifier(output_dir):
    flag = 0
    while flag == 0:
        ide = ''.join(random.choice(
                string.ascii_uppercase + string.ascii_lowercase + string.digits
                ) for _ in range(8))
        if not os.path.exists(output_dir+"job_"+ide+"/"):
            flag = 1
    return ide

def get_stats_for_charts(stats_file):
    with open(stats_file, "r") as f:
        d = json.load(f)
    return (d["table_columns"], d["not_found"], d["no_int_prots"])

@app.route("/")
@app.route("/index")
def index():
    maintenance = False
    if maintenance:
        return render_template(maintenance_template)
    else:
        return render_template(index_template)

@app.route("/help")
def help():
	return render_template(help_template)

@app.route("/output", methods = ["POST", "GET"])
# @line_profile
def output():
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

        job_id = get_unique_random_identifier(output_dir)
        job_dir = output_dir+"job_"+job_id+"/"
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
# @line_profile
def run_job(job_id):
    job_dir = output_dir+"job_"+job_id+"/"
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
                        MAIN_OUTPUT_DIR=output_dir,
                        CUSTOM_ID=job_id,
                        DATA_DIR=data_dir,
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
            if error==1:
                return render_template(error1_template)
            elif error==2:
                return render_template(error2_template)
    else:
        print_log(job_id, "Files exist")

    # Read stats file
    (table_columns, not_found, no_int_prots) = get_stats_for_charts(job_dir+stats_file)
    return render_template(results_template,
                       graph_json="jobs/"+"job_"+job_id+"/"+graph_json,
                       ints_json="jobs/"+"job_"+job_id+"/"+table_file,
                       stats_json="jobs/"+"job_"+job_id+"/"+stats_file,
                       not_found=not_found,
                       no_int_prots=no_int_prots,
                       table_columns=table_columns
                       )
