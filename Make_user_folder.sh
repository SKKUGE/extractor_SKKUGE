#!/bin/zsh

user=KJS
project=PELIB

[ ! -d ./Input ] && { mkdir ./Input; }
[ ! -d ./User ] && { mkdir ./User; }
[ ! -d ./Output ] && { mkdir ./Output; }

[ ! -d ./Input/${user} ] && eval mkdir ./Input/${user};
[ ! -d ./Input/${user}/FASTQ ] && eval mkdir ./Input/${user}/FASTQ;
[ ! -d ./Input/${user}/FASTQ/${project} ] && eval mkdir ./Input/${user}/FASTQ/${project};
[ ! -d ./Input/${user}/Reference ] && eval mkdir ./Input/${user}/Reference;
[ ! -d ./Input/${user}/Reference/${project} ] && eval mkdir ./Input/${user}/Reference/${project};

[ ! -d ./User/${user} ] && eval mkdir ./User/${user};

# TODO: What is the usage?
# > ./User/${user}/${project}.txt




