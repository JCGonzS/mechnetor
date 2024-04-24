import os
import sys
from dotenv import dotenv_values


HERE_PATH = os.path.dirname(os.path.abspath(__file__))

# Load configuration parameters from .env file
config = dotenv_values(dotenv_path=f"{HERE_PATH}/.env")


def get_paths():
    path = {
        "data": f"{HERE_PATH}/static/data/",
        "jobs": f"{HERE_PATH}/static/jobs/",
        "log": f"{HERE_PATH}/log.txt"
    }
    return path


def get_templates():
    template = {
        "index": "index.html",
        "maintenance": "maintenance_index.html",
        "help": "features.html",
        "error1": "input_error.html",
        "error2": "toobig_error.html",
        "results": "results.html"
    }
    return template
