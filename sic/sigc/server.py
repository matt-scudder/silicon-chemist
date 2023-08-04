"""
This is the server for SiGC. It should contain a minimum of logic
not related to routes for the website, and offload as much
processing as possible to the other modules.
"""

import os

from flask import Flask, render_template, request
from flask_debugtoolbar import DebugToolbarExtension

from sic import sic

app = Flask(__name__)
app.config.from_prefixed_env()
if (app.config["SECRET_KEY"]):
        toolbar = DebugToolbarExtension(app)
@app.route("/")
def jsme_display():
        """
        Displays the JSME interface for creating reactant and product sets.
        """

        return render_template("reactions/reaction_sketch.html") #add in template vars as needed...

@app.route("/submit_reaction",methods=["POST"])
def accept_reaction():
        """
        Handles reactant and product input from users and runs SiC³ on them.
        """

        if request.method == "POST":
            return sic.find_mechanism(request.json["reactants"],request.json["products"],solv=request.json["solvent"]).replace("\n","<br>")
