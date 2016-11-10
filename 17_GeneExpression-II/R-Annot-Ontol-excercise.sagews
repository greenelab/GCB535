︠80e0975f-47ca-4c44-82f3-a113242d9e7bas︠
%auto
%default_mode r
%python
r.set("workdir", "'%s'" % os.path.abspath(os.curdir))
%r
setwd(workdir) 
︡dad8273e-3a98-4d8a-825f-8afb45a4e3b2︡{"auto":true}︡
︠3b12e69e-86a7-4689-957d-46755154f42f︠
# This is your assignment and homework. (R/Bioconductor for Gene-expression Analysis: Annotation, QC, Ontology). For this assignment, we're going to be reusing
# some of the old code you built from yesterday's assignment, and add a couple of features.
︠d05f02ed-c13a-4e9e-959e-dbc8591802f0︠
#1a. First, let's copy over the relavent bits of the pipeline you built yesterday over to this worksheet.
# - load the libraries we will use for our analysis here [Using R]
### NEW STEP: load the biomaRt library, which we will be using to extract annotations from ensembl
### ENTER YOUR CODE BELOW



︡4f267b8f-7f98-4399-bcd1-750f7628c46c︡︡{"done":true}
︠5d0b542c-d6b1-43d7-bc26-08bb4fdbfb00︠
# - load phenotype data into R
### NEW STEP: For this analysis, let's compare NASH vs. CONTROL. We'll have to make a new .csv phenotype file for that!
### ENTER YOUR CODE BELOW



︡b5015cdb-b813-4557-bef1-3e743a981bd9︡︡{"done":true}
︠68afbcaf-6ec0-4dcb-98a2-72f62794d55f︠
# - Next, we need to load our CEL file data into R.
### NEW STEP: rather than redownload the data, let's just change the path to where the data lives (Hint: use UNIX. ".." moves up a directory.)
### ENTER YOUR CODE BELOW



︡641d3256-a531-47f0-b050-6f24e2754b95︡︡{"done":true}
︠65cf7844-b73b-40c5-bd57-bf16e0619126︠
#1b. Continue your pipeline
# - load CEL data into R
### ENTER YOUR CODE BELOW 



# - normalize the data using rma
### ENTER YOUR CODE BELOW 



# - plot the normalized output
### ENTER YOUR CODE BELOW [execute]




︡dcfdd2e3-c657-40c8-a13e-376330ee89fe︡︡{"done":true}
︠af8adee0-42bf-4dd2-a122-1a14db8f6262︠
# 1c. Continue your pipeline and create the design matrix for analysis. Make sure the CONTROLS are reference (=0), and NASH is the treatment (=1)
### ENTER YOUR CODE BELOW [execute]



︡e0d58ab2-f536-4ebe-9452-169fcf172dd1︡︡{"done":true}
︠10c637ae-ce42-449f-a568-84a28b347bc3︠
# 1d. Now, perform gene expression analysis, obtain the top 100 results, print the top five hits, and output results to "mytop100results_nashvctrl.txt".
### ENTER YOUR CODE BELOW [execute]




︡a229f36a-574b-417e-905c-3be1ca849996︡︡{"done":true}
︠7ef2ec31-d49e-4435-93a1-8bb5ea7361de︠
# 2a. OK, now let's obtain some annotation information
# - create a variable called mymart, which loads the appropriate library from the ensembl database for mouse.
# - obtain a list of the affy probe ids for the top 100 results, stored in a variable called pidsTophits
# - print the first 5 entries from pidsTophits

###ENTER YOUR CODE BELOW [execute]




︡68a91c19-8481-43f3-b37e-c725f125d5ce︡︡{"done":true}
︠d66fbd4b-494f-4ecf-9d42-206889fcce26︠
# 2b. Continue obtaining annotation information
### note that getting all of the annotations will take some time!
# - obtain the following annotations from ensembl stored in a new variable called myannot for your list of affy probe ids:
# - probeids (for the array technology used for this experiment), chromosome, start_position, end_position, gene symbol (mouse), entrez gene id, Description, Hugo gene symbol
# - print out to screen for the first 5 entries: probeid, chr, start, end, mgi, and description

###ENTER YOUR CODE BELOW [execute]




︡86b574dd-154f-4e57-ab6c-39b0fca6f7c9︡︡{"done":true}
︠4e0a36e6-3a3d-46be-888d-32c5fab58add︠
# 2c. Write two annotations to file:
# - one which reports everything but the description, without quotes, and separated by commas (useful for data parsing) to a file called tt-annot-forparse.csv
# - one which reports everything, with quotes, separated by commas (useful for reading in excel) to file called tt-annot.csv

###ENTER YOUR CODE BELOW [execute]




︡af312a95-a1f3-430e-a6dd-7f07bb67c3fe︡︡{"done":true}
︠cf0cdf72-268f-48d4-9abd-fc1db7ffe345︠
# 3a. Now that we have annotations, we need to merge them together with the table that contained our association information. 
# However, notice that the annot matrix lists the probe ids in a different order than the statistical association table
#
# - create a new column that binds the pidsTophits variable with the tt variable (hint: cbind). Store into a new variable called tt_ids.
# - rename the column containing pidTophits to affyprobeid in the variable tt_ids (hint: colnames).

###ENTER YOUR CODE BELOW [execute]




︡840707b5-d63b-4c31-a78c-2bcf421b6304︡︡{"done":true}
︠243501f5-b87f-4fd7-b226-34162af7e9cb︠
# 3b. Now, we need to prepare and merge the annot file with tt_ids:
# - rename the column in annot that contains your probeid to "affyprobeid"
# - output the first 5 entries in annot and make sure the name has been changed
# - Use R to merge tt_ids together with annot, by affyprobeid, into a new table called tt_plus_annot

###ENTER YOUR CODE BELOW [execute]



︡dd1295b0-cb72-4a02-96ba-ef3a9ac167d9︡︡{"done":true}
︠41532509-998e-4ab5-a2c5-3999ff3367e5︠
# 3c. Let's resort the table to get the most significant results on top. You can use the order() command for this
# - use the R function order to sort the the table by adj.P.Val, and store the output to tt_annot_s
# - output the results of the first 5 row of tt_plus_annot

###ENTER YOUR CODE BELOW [execute]




︡0b6cdb13-9732-44f5-bb9d-5b411c1a64de︡︡{"done":true}
︠a87e89e6-20fa-4a11-9834-8c2e82266aaf︠
# 3d. What are the genes listed in the top 20 results? 
#
# YOUR ANSWER HERE >>> 



︠6429f400-bb16-4fee-9456-618d55266e77︠
%python 
##DO NOT DELETE the line above here
### 4a. OK, let's use python on our top 100 results, and get a list of the genes that were found from our analysis.
# Use python to print out the list genes from the top 100 results file tt-annot-forparse.csv,to a file called 'mygenelist.txt'

### ENTER YOUR CODE BELOW [execute]



︡c5421f9d-3dc0-4bda-88e9-1ab99db8840c︡︡{"done":true}
︠79c332b1-0911-45c3-ad1f-c2c6d73dd427︠
#4b. Do we have unique entries in our list? If so, use UNIX to determine which entries were duplicated and how many duplicates were present

### YOUR ANSWER HERE >>>


#4c. Use UNIX to create a list of only unique entries. Save this result to a file called "mygenelist_uniq.txt"


### ENTER UNIX COMMAND HERE >>>


︠26774586-32f1-4fce-93f5-d1e7a7990fd6︠
# 5. Now, let us use WebGestalt to perform some ontology analysis using the unique set of genes created above.

# a. Download the list of genes from sagemathcloud to your computer, and copy-paste the list of genes into WebGestalt. 
# Does that tool find all of the genes that you provided? Why or why not?

# YOUR ANSWER HERE >>>

# b. How would you address the issue in part a?

# YOUR ANSWER HERE >>>

# c. What background gene list will you select? Why is this important?

# YOUR ANSWER HERE >>>

# d. Selection some enrichment analyses and perform. Do you find anything associated at compelling levels of statistical support? if so, what? Any results that make sense for the comparison that you made? 

# YOUR ANSWER HERE >>> 

# e. Here, we selected the top 100 genes for enrichement analysis. Was this choice warrented? Why or why not? Having performed this analysis, what might you do differently?

# YOUR ANSWER HERE >>> 
︠2df21eb3-0aca-48fe-81f6-27621907a723︠
# 6. Run the the sessionInfo() command to document the packages you used to make the analysis work for posterity and reproducibility.

###ENTER YOUR CODE HERE [execute]



︡05b62a80-dd24-4f88-94f7-ae074f1018d6︡︡{"done":true}
︠e8da7a8e-5253-4cd3-a311-380a8475a476︠









