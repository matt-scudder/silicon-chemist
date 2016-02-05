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
from flask import request
import sic #no relative import - this should only be run by runserver.py

app = Flask(__name__)


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
            return sic.find_mechanism(request.json["reactants"],request.json["products"],solv=request.json["solvent"])
