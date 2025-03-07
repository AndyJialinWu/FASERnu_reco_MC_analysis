#!/bin/bash


inputFileName=PhysicsNTUP_ML.root
outputFileName=TMVAOutput.root


# Run the TMVA classification macro
ScriptName=TMVAClassification.C
root -l $ScriptName'("'${inputFileName}'", "'${outputFileName}'")'
 


