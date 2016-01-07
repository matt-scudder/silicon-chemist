#!/usr/bin/python
#coding=utf-8
"""
This is the server for SiGC. It should contain a minimum of logic
not related to routes for the website, and offload as much
processing as possible to the other modules.
"""
import traceback #never know when you'll need this.
from flask import Flask
from flask import render_template
from .. import sic #where the magic happens - note the relative import!

app = Flask(__name__)


@app.route("/")
def jsme_display():
	"""
	Displays the JSME interface for creating reactant and product sets.
	"""
	return False

@app.route("/submit_reaction")
def accept_reaction():
	"""
	Handles reactant and product input from users and runs SiC³ on them.
	"""
	return False

if __name__ == "__main__":
	app.run(debug=True) #NOTE: Change this to port 80, host 0.0.0.0 and debug=False for final release
