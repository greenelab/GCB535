︠80e0975f-47ca-4c44-82f3-a113242d9e7bas︠
%auto
%default_mode r
%python
r.set("workdir", "'%s'" % os.path.abspath(os.curdir))
%r
setwd(workdir)
︡5a069e15-72f2-433b-9642-d2046b9ee829︡{"auto":true}︡
︠841e5c4d-4e90-4eb2-b4e4-e687e9391b47︠
# This is your in-class and homework assignment. The objective of this work is to implement your own pipeline for the analysis of some new array data.
# By the end of the assignment, you will have constructed your own gene expression analysis pipeline.

# Let's take a look at a new data set to analyze: GSE35961

# 1a. Describe in your own words, the treatment groups, number of samples, and subjects that were characterized in this experiment.

# YOUR ANSWER HERE>> 

# 1b. What is the citation attached to this paper? 

# YOUR ANSWER HERE>> 

# 1c. What is NASH? What is Metformin? What is the hypothesis that this experiment was designed to test? 

# YOUR ANSWER HERE>> 

︠374ec30e-ddf8-4b70-b898-4c0ec0fe5b11︠
### 2. Let's use R to download this data set, use UNIX to prepare our associated input files and organize our directory.

## 2a. load the libraries we will use for our analysis here [Using R]

### ENTER YOUR CODE BELOW [execute]




︡edad7093-3c28-4acb-81b8-992651bd6301︡︡{"done":true}
︠d91f91d8-f3ee-4773-ba0f-128b0d22b628︠
## 2b. Download the data set GSE35961 from GEO [Using R]

### ENTER YOUR CODE BELOW [execute]




︡8912d24e-1c4f-4bc5-82c5-036b105995aa︡︡{"done":true}
︠4c10cbcb-e6b6-4555-ab1c-e3eef5b4f896o︠
#2c. Now, we need to process the data that we downloaded.
# - Expand the GSE35961_RAW.tar archive [UNIX]
# - Uncompress all of the .gz files [UNIX]
# - Delete the GSE35961_RAW.tar file [UNIX]

### ENTER YOUR UNIX COMMANDS for C BELOW [No need to execute in the worksheet]



︠52c85c3a-3d9a-4abc-b0a6-5141c823a8cd︠
#3. Next, we need to prepare our phenotype file for analysis. For this assignment, I would like to compare the samples of NASH to the NASH treated with metformin.
# - Prepare a phenotype.csv file which references the appropriate CEL files for analysis
# - Load this phenotype file into R and store it into a variable called phenoData
# - return a summary for the variable phenoData using the pData command

### ENTER YOUR CODE BELOW [execute]




︡7f12a280-b273-447e-b89d-621b840ccb7b︡︡{"done":true}
︠f792cf40-4580-4140-ab64-c8897493140e︠
# 4. Next, we need to load our CEL file data into R.
# - read a list of cel files into an object call celFilelist
# - create a new variable called celFiles which contains only the CEL files that you want to analyze (ie present in the phenoData variable in #3)
# - report the contents of the celFiles variables
# - read intensity data from the list of celFiles into a variable called affyRaw 

### ENTER YOUR CODE BELOW [execute]




︡22504e5d-0daf-4a9c-85f3-65417a88e244︡︡{"done":true}
︠65cf7844-b73b-40c5-bd57-bf16e0619126︠
# 5. Let's now take a look at the data
# 5a create a boxplot of the intensity data for all CEL files loaded

### ENTER YOUR CODE BELOW [execute]




︡80de9e50-77a9-4984-b91e-2802a59e9b95︡︡{"done":true}
︠545ed812-0b26-4f63-a339-4ddeea32f2bd︠
# 5b. create a histogram of the intensity data for all CEL files loaded

### ENTER YOUR CODE BELOW [execute]




︡083cd4de-dd61-4609-8ae1-151a5227b0a3︡︡{"done":true}
︠d36103f5-bef4-442f-ac7d-2a7de6bf205f︠
# 5c. Now, normalize the intensity data using RMA, and store the normalized data in a variable called genenorm
# 5d. Then plot a boxplot for the normalized intensity data

### ENTER YOUR CODE BELOW [execute]




︡67affa32-6fa1-4c18-ac60-0bdff3e1bf71︡︡{"done":true}
︠af8adee0-42bf-4dd2-a122-1a14db8f6262︠
# 6. Prepare your design matrix for the analysis
# - create list of your treated/untreated sample and store in a variable called group, and print
# - create the design matrix using the group variable, and store in a variable called design, and print

### ENTER YOUR CODE BELOW [execute]




︡d6999d2c-edaa-433a-8eaf-743c2e0bac77︡︡{"done":true}
︠10c637ae-ce42-449f-a568-84a28b347bc3︠
# 7. Now, perform gene expression analysis
# - analyze the normalized data using the lmFit() function, your design matrix, and store in a variable called fit
# - apply empirical bayes correction using eBayes(), store in a variable called efit
# - get the top 250 results, sorted by P-value and store in a table called tt
# - print the top 5 results to this notebook

### ENTER YOUR CODE BELOW [execute]




︡c7644ff8-320c-4547-b702-0e18ba4e99d5︡︡{"done":true}
︠e182a2fa-a1bf-4299-bcf8-366096bb13ae︠
# 8a. Generate some summary outputs for your top results
# - make a table of your top 250 results to a file called mytop250results.txt

### ENTER YOUR CODE BELOW [execute]




︡bf9f7fa8-f92b-405b-bf12-ab48bfbf28d0︡︡{"done":true}
︠050e3e40-2df6-485f-8104-50a90376dd73︠
# 8b. Make a volcano plot of all of your results.

### ENTER YOUR CODE BELOW [execute]




︡7ba92652-d194-4f5a-838d-5ff0e82aa38b︡︡{"done":true}
︠e78cc857-da76-4c8d-a350-1a85760f451f︠
# 9. One final item in your pipeline to finish it off complete.
# - run the the sessionInfo() command to document the packages you used to make the analysis work for posterity and reproducibility.

###ENTER YOUR CODE BELOW [execute]




︡3bce1577-7187-49f0-9d4d-7cce219c2c30︡︡{"done":true}
︠b6d2f2f8-b64b-48e5-bc1a-50e186b640a2︠
# 10. Some final statistical interpretation questions about the analysis you performed here:

# a. How many probes did you analyze in total? (Hint: look at the eFit object)

# YOUR ANSWER HERE>> 

# b. How many probes returned an adjusted p-value less than 0.05? 

# YOUR ANSWER HERE>> 

# c. Are you convinced that any probe is differently expressed between NASH and NASH+Metformin conditions? If yes, which ones, and why? if not, why not?

# YOUR ANSWER HERE>> 









