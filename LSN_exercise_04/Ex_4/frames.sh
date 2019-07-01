#!/bin/bash

# Load config values
file="run.dat" #the file where you keep your string name

name=$(cat "$file")        #the output of 'cat $file' is assigned to the $name variable

mkdir

echo $name

