#!/bin/bash
#Automates the "build" procedure required due to 1. lack of Netbeans 2. issues with CDK versioning

cd silly-build
rm -r aamtool
rm -r uk/ac/ebi/reactionblast
rm -r uk/ac/ebi/centres
cp -r ../../build/classes/aamtool ./
cp -r ../../build/classes/uk/ac/ebi/reactionblast uk/ac/ebi
cp -r ../../build/classes/uk/ac/ebi/centres uk/ac/ebi
jar cfe usableRDT.jar aamtool.ReactionDecoder .
cp usableRDT.jar ../usableRDT.jar
