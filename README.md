# TextLDS
Code for David Belanger and Sham Kakade "A Linear Dynamical System Model For Text." ICML 2015. 

## Setup
1) As is standard in many NLP applications, you'll have to convert all of your words to ints, and keep a file with the word-to-index mapping somewhere. This needs to be done in a separate script. 

2) Convert your text data into a list of Matlab .mat files, where each .mat file contains a lot of text (potentially for many documents concatenated together). The file should contain a single field, called 'ints' that's just a 1-dimensional array of ints. It's length will be the number of tokens observed. 


3) Then, call:

getCoocurrenceAtLags(workingDir,id,V,filenameList)

workingDir is where to write output data to. 

id is a name for the dataset.

V is the size of the vocabulary. 

filenameList is a list of the files made in step (2). 

4) Then, call 

getWhitenedCoocurrenceAtLags(dataDir,kappa)

kappa is a smoothing term for whitening the data. It should be interpreted as a pseudocount for each word. 

5) The main function for training is learn_word_lds.m

You'll specify the path to where the data was generated in (4) as the  'datasetname' argument. 

6) Use Evaluation/onehot_Steady_State_KF.m to perform likelihood computation, filtering, and smoothing on held-out data. The 'y_' argument is a 1d dense array of int ids for the words. If you have a multi-sentence document, either concatenate all sentences together or call this function once for every sentence. In our experiments, we did the former. 

