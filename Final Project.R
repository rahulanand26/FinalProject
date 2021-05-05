# Final Project Code
# Rahul Anand and Sai Pola

# Please be patient with this code, it will take a long time to run due to it
# having to go through 1,250,018 elements in the PDB database. Even when broken 
# into subsets, the code still takes many hours to run in total on my machine.



# Import the text file containing all protein sequences. The file is broken down as follows:
# "> Protein name and code"
# "Protein sequence"
x <- readLines('C:/Users/rahul/Downloads/pdb_seqres.txt')
head(x)



# Vector to store the percentage of prolines in residues
proline_count <- c()

# Vector to store the sequences
sequences <- c()

q <- 0

# Separates sequences from massive file
for (line in x){
  
  # progress counter
  q <- q + 1
  print(q)
  
  
  # Skip the lines with protein names/codes
  if (startsWith(line, ">") == TRUE){
    next
  }
  
  else{
    # Calculate number of prolines / length of sequence and add to vector
    a <- lengths(regmatches(line, gregexpr("P", line)))
    b <- a / nchar(line)
    proline_count <- c(proline_count, b)
    
    # Add the sequences to a list of sequences separate from the names
    sequences <- c(sequences, line)
    
  }
}



# PLOT A VIOLIN PLOT OF PROLINE RATIOS
install.packages('vioplot')
library(vioplot)
vioplot(proline_count, col = 'green', ylab = 'Proline ratio', main = 'Proline Composition Ratio of PDB Sequences')



# Separate PRPs and Non-PRPs
PRPs <- c()
nonPRPs <- c()
for (i in 1:length(sequences)){
  
  # progress counter
  print(i)
  
  # found a PRP
  if (proline_count[i] >= 0.25){
    # add this polypeptide to the vector of PRPs
    PRPs <- c(PRPs, sequences[i])
  }
  
  # not a PRP
  else{
    # add this polypeptide to the vector of non-PRPs
    nonPRPs <- c(nonPRPs, sequences[i])
  }
}



# Isolate all unique PRPs (remove repeats from the PDB as Matthew instructed)
unique_prps <- unique(PRPs, incomparables = FALSE)
# We find 624 unique PRPs out of the original 1875 identified

# Isolate all unique non-PRPs (remove repeats from the PDB as Matthew instructed)
unique_non_prps <- unique(nonPRPs, incomparables = FALSE)
# We find 137080 unique non-PRPs out of the original 623134 identified



# Here we are getting our reference sequence from colostrum (the substance in
# breastmilk known for its health benefits and proline-rich polypeptide abundance)
# The reason why we are comparing every sequence to this rather than comparing all 
# versus all is because with over 62,000 total protein chains, we don't have the 
# computer power or time to execute every possible combination into RMSD. Keep in
# mind that RMSD allows only comparison of 2 sequences at a time.

refseq_fasta <- readLines('C:/Users/rahul/Downloads/lactotransferrin.fasta')
# Lactotransferrin is our reference polypeptide as it is common in colostrum.
lactotransferrin <- refseq_fasta[2]

install.packages('Biostrings')
library(Biostrings)
?PairwiseAlignments

# Get length of lactotransferrin
length_lacto <- nchar(lactotransferrin)



# Create Vector to hold all PRP pairwise alignment objects
PRP_PAs <- c()

# counter
q <- 0
# iterate through all the unique PRPs
# We use global alignment in order to use the entire sequences
for (s in unique_prps){
  
  # the strings are identical lengths so we can do pairwise alignment
  if (nchar(s) == length_lacto){
    PRP_PAs <- append(PRP_PAs, PairwiseAlignments(s, lactotransferrin, type = "global"))
  }
  
  # the PRP string is shorter than lactotransferrin so we need to increase string 
  # to match
  else if (nchar(s) < length_lacto){
    # pad the string with dashes
    s <- str_pad(s, length_lacto, side = "right", pad = "-")
    
    # run pairwise alignment now that they are the same length
    # append it to the vector
    PRP_PAs <- append(PRP_PAs, PairwiseAlignments(s, lactotransferrin, type = "global"))
  }
  
  # the PRP string is longer than lactotransferrin so we increase the lactotransferrin
  # length to match
  else{
    # Create a padded version of the lactotransferrin sequence to run pairwise alignment
    lacto_padded <- str_pad(lactotransferrin, nchar(s), side = "right", pad = "-")
    
    # run pairwise alignment now that they are the same length
    # append it to the vector
    PRP_PAs <- append(PRP_PAs, PairwiseAlignments(s, lacto_padded, type = "global"))
  }
  
  # progress counter
  q <- q + 1
  print(q)
}


# Now we will run through all the pairwise alignment objects between PRPs and 
# lactotransferrin to get our results of Percent Sequence Identity (PID)
# function from Biostrings package
?pid

# Create a vector to hold all the PIDs 
# PID = 100 * (identical positions) / (aligned positions + internal gap positions)
PRP_PIDs <- c()

for (i in 1:length(PRP_PAs)){
  # Append the found percent sequence identity to the vector
  PRP_PIDs <- append(PRP_PIDs, pid(PRP_PAs[i]))
  
  # progress counter
  print(i)
 
}



# we will now do the same 2 steps we did for PRPs, but now for the non-PRPs
# Create Vector to hold all non-PRP pairwise alignment objects
non_PRP_PAs <- c()

# counter
q <- 0
# iterate through all the unique non-PRPs
# We use global alignment in order to use the entire sequences
for (s in unique_non_prps){
  
  # the strings are identical lengths so we can do pairwise alignment
  if (nchar(s) == length_lacto){
    non_PRP_PAs <- append(non_PRP_PAs, PairwiseAlignments(s, lactotransferrin, type = "global"))
  }
  
  # the non-PRP string is shorter than lactotransferrin so we need to increase string 
  # to match
  else if (nchar(s) < length_lacto){
    # pad the string with dashes
    s <- str_pad(s, length_lacto, side = "right", pad = "-")
    
    # run pairwise alignment now that they are the same length
    # append it to the vector
    non_PRP_PAs <- append(non_PRP_PAs, PairwiseAlignments(s, lactotransferrin, type = "global"))
  }
  
  # the PRP string is longer than lactotransferrin so we increase the lactotransferrin
  # length to match
  else{
    # Create a padded version of the lactotransferrin sequence to run pairwise alignment
    lacto_padded <- str_pad(lactotransferrin, nchar(s), side = "right", pad = "-")
    
    # run pairwise alignment now that they are the same length
    # append it to the vector
    non_PRP_PAs <- append(non_PRP_PAs, PairwiseAlignments(s, lacto_padded, type = "global"))
  }
  
  # progress counter
  q <- q + 1
  print(q)
  
}


# Now we will run through all the pairwise alignment objects between non-PRPs and 
# lactotransferrin to get our results of Percent Sequence Identity (PID)
# function from Biostrings package
# Create a vector to hold all the PIDs 
# PID = 100 * (identical positions) / (aligned positions + internal gap positions)
non_PRP_PIDs <- c()

for (i in 1:length(non_PRP_PAs)){
  # Append the found percent sequence identity to the vector
  non_PRP_PIDs <- append(non_PRP_PIDs, pid(non_PRP_PAs[i]))
  
  # progress counter
  print(i)
  
}

# Remove the identical lactotransferrin instance from the PID list because comparing
# it to itself will only skew the data
non_PRP_PIDs <- non_PRP_PIDs[ non_PRP_PIDs != 100 ]


# Make plots below
# Plot a violin plot of the PIDs
library(vioplot)
vioplot(PRP_PIDs, non_PRP_PIDs, col = 'red', ylab = 'Percent Sequence Identity', xlab = ' PRP PIDs                                                                         non-PRP PIDS ', main = 'Percent Sequence Identity to Lactotransferrin in Colostrum between PRPs and non-PRPs')    

# Get Mean, Median, Min, Max
summary(PRP_PIDs)
summary(non_PRP_PIDs)


