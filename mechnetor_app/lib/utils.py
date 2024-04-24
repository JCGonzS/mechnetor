import re
import os
import gzip
import datetime
import random
import string


def open_file(input_file, mode="rt"):
    """ Open file Zipped or not
    """

    if re.search(".gz$", input_file):
        infile = gzip.open(input_file, mode)
    else:
        infile = open(input_file, mode)
    return infile


def print_log(job_id, msg):
    st = "[{}]".format(datetime.datetime.now())+" [JOB ID: "+job_id+"]"
    print(st, msg)
    return


def get_unique_random_identifier(output_dir):
    flag = 0
    while flag == 0:
        ide = ''.join(random.choice(
                string.ascii_uppercase + string.ascii_lowercase + string.digits
                ) for _ in range(8))
        if not os.path.exists(output_dir+"job_"+ide+"/"):
            flag = 1
    return ide
