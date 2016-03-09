#!/bin/bash
# Because the people behind AAMTool are using Netbeans, they don't have to think about how to
# really build their stuff on the command line. Therefore, this script exists
# to merge the CDK library into the ReactionDecoder jar, so that when you run it with java -jar,
# it doesn't complain about missing classes because it decides that the classpath of the jar
# is the only classpath in the universe.

#first compile and build jar
ant compile
ant jar 

#first wipe the temp directory
rm -rf dist/tmp-merge/*

#now copy the new jars which you hopefully made, into that directory
cp dist/ReactionDecoder.jar dist/tmp-merge/ReactionDecoder.jar
cp dist/cdk-1.4.19.jar dist/tmp-merge/cdk-1.4.19.jar

#go to the other directory and extract both jars
cd dist/tmp-merge
jar xvf ReactionDecoder.jar
jar xvf cdk-1.4.19.jar

#delete original jars

#rm ReactionDecoder.jar
#rm cdk-1.4.19.jar

#create new jar with aamtool entrypoint
jar cvfe TrueReactionDecoder.jar aamtool.ReactionDecoder .

#copy it back
cp TrueReactionDecoder.jar ../

#celebrate
echo "Successfully created jar."
