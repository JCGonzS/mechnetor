import mechnetor_app.main

# ap = argparse.ArgumentParser()
# ap.add_argument("-p", "--prots", required=True,
#                 help="File with protein input as protein identifier or FASTA sequences (Required)")
# ap.add_argument("-m", "--muts", required=False,
#                 help="File with mutation input in Mechismo format")
# ap.add_argument("-sp", "--species", required=False, default="HUMAN",
#                 help="Organism (default: Hsa)")
# ap.add_argument("-ai", "--add_ints", required=False, default=0,
#                 help="Number of additional interactors per protein. Any number or 'all' (default: 0)")
# ap.add_argument("-id", "--job_id", required=False, default=False,
#                 help="Custom job name (random by default)")
# ap.add_argument("-log", required=False, default=False,
#                 help="Print log file")
# args = vars(ap.parse_args())
#
# prots = open_file(args["prots"]).read()
# muts = ""
# if args["muts"]:
#     muts = open_file(args["muts"]).read()
#
# if args["log"]:
#     sys.stdout = open(args["log"], 'w')
#
# if args["species"] == "HUMAN":
#     # db = "piv_human"
#     db = "mechnetor_human"
# else:
#     # db = "piv_all"
#     db = "mechnetor_all"
#
# main(INPUT_1=prots, INPUT_2=muts, ORG=args["species"],
#      ADDITIONAL_INTERACTORS=args["add_ints"],
#      CUSTOM_ID=args["job_id"], PSQL_DB=db,
#      MAKE_NETWORK=True, TABLE_FORMAT="json", CMD_LINE=True,
#      RUN_IPRETS=False)
#
# sys.exit()
