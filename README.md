# GCB535
This repository contains the course materials for GCB 535 taught at the University of Pennsylvania.

## Sage Math Cloud
Computing needs for the course were handled via SageMathCloud. We used a medium course plan, and requested disk space upgrades to support the hard drive space required for the class. If you use SageMathCloud then make sure that you don't 'assign' items to students until you are absolutely sure they are ready for your course.

## Setup
Some exercises require setup. For example, code to simulate data used in the k-means clustering exercise (26_Machine_Learning_I_assignment/kmeans-population.csv). Code to generate such files will be found in folders prefixed with the word "SETUP".

## Python Programming Exercises
If you are only interested in the python programming portions, check into @sarahmid's repository sarahmid/python-for-genomics-miniseries. The exercises in this course share a common origin with those exercises.

## Data on SageMath Cloud
For modules that have big data files (e.g. the ChIP-seq modules), it is best to create two separate assignments in SMC: one that contains the data itself and is never 
collected by the TAs for grading, and another one that contains the assignment notebooks and is collected by the TAs. This is because files are copied whenever assignments
are collected and so if an assignment contains large data files that never change, they will be duplicated for every student in the class, wasting a lot of space. 