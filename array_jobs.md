### Useful linux/bash commands for running array jobs

Just as a note, Linux is an [operating system](https://en.wikipedia.org/wiki/Operating_system) and Bash is a [command-line interface](https://en.wikipedia.org/wiki/Command-line_interface)

UF HiPerGator has a lot of resources online for running jobs, for example for [submitting an array job](https://help.rc.ufl.edu/doc/SLURM_Job_Arrays). I want to add onto these resources by sharing some helpful commands/lines when I’ve run array jobs

## The commands
To start, I have found these commands to be particularly helpful when running jobs on the cluster:

```{r analysis, results="markup"}# ls
echo
cut
head
tail
rev
rm
````

_Man_ is also a good one to know; Less so within the scripts but it’s “an interface to the on-line reference manuals”. On the cluster, you can type “man” then any of these commands to get the description and the flag options (e.g. “man ls”)

Here is what each of these commands is defined to do:
```{r analysis, results="markup"}
ls – list directory contents
echo – display a line of text 
cut -  remove sections from each line of files
head – output the first part of files
tail – output the last part of files
rev – reverse lines of a file or files
rm – remove files or directories
````
Also, included will be the flags one can use and their descriptions (if it’s a required flag with or without a default, if it needs to be a numerical value, etc.)

## Running array job
As mentioned above, HiPerGator has online resources for [submitting an array job](https://help.rc.ufl.edu/doc/SLURM_Job_Arrays). Here I will provide a quick overview.

Array jobs are useful for running the same (or similar) commands on multiple files. To run an array job, we specify this in our sbatch script using the following line:

```{r analysis, results="markup"}
#SBATCH --array
```

If you want to run the same commands on the first two files (by default in [lexigraphical order] https://en.wikipedia.org/wiki/Lexicographical_order) in the same directory, you add this line to your sbatch script: 

```{r analysis, results="markup"}
#SBATCH --array 1-2
````

If you want to run jobs 1-100 you’d use:

```{r analysis, results="markup"}
#SBATCH --array 1-100
````

If you only want to run jobs 1, 5, 17:

```{r analysis, results="markup"}
#SBATCH --array 1,5,17
````

If you want to run jobs on 1-100, but only submit 10 jobs at a time (useful when sharing an account with others):

```{r analysis, results="markup"}
#SBATCH --array 1-100%10
````

## One input, one output

Let’s say we want to build gene trees on 100 gene alignments using RAxML. For the array we'll use:

```{r analysis, results="markup"}
#SBATCH --array 1-100%10
````
So, we’re running 100 files and submitting them 10 at a time.

Here’s our RAxML command:

```{r analysis, results="markup"}
raxmlHPC-PTHREADS-SSE3 -f a -x $RANDOM -N 100 -T 5 -p $RANDOM -s INPUT_FILE_VARIABLE_HERE -n OUTPUT_FILE_VARIABLE_HERE -m GTRGAMMA
````

What we’ll need to do to is make variables for each of the input alignments and to output trees. For RAxML, the flags we need to make variables for are “-s” and “-n”, respectively. Now let's write something to call our particular input files. Let’s say our input alignment names are “cluster1.pep.filtered.fa, cluster2.pep.filtered.fa,… cluster100.pep.filtered.fa”. What we want to do is make a variable that will loop through each of our input alignments that we designate in our array. Here’s how we can make our input variable:

```{r analysis, results="markup"}
input_aln=`ls *.filtered.fa | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
````

In short, this variable captures the first 100 files (because that's what defined in our array line) that end in ".filtered.fa". More specifically, what we are first doing is using _ls_ to list everything in our directory that ends in “.filtered.fa”. We then pipe this to _head_ (the “|” is the pipe command) and we give the variable [(SLURM generates this one)](https://slurm.schedmd.com/documentation.html) “$SLURM_ARRAY_TASK_ID” to the flag “-n” to get the first file our array (in this case file 1) and then pipe this to _tail_ to get our last file (file 100). So, now we have our input variable for RAxML:

```{r analysis, results="markup"}
raxmlHPC-PTHREADS-SSE3 -f a -x $RANDOM -N 100 -T 5 -p $RANDOM -m GTRGAMMA -s $input_aln -n OUTPUT_FILE_VARIABLE_HERE
````

Where each one of the alignment files will be passed to our RAxML commands one at a time. It’s important to remember when calling a defined variable to start with a “$” then give the variable name. Next we’ll want to make an output variable. We could very simply add the file extension “.tre” to our input file names, like so:

```{r analysis, results="markup"}
out=`echo $input_aln`.tre
````

Where _echo_ prints our alignment name (e.g. cluster1.pep.filtered.fa), to which at the end we are adding on “.tre”. Our output file name would be “cluster1.pep.filtered.fa.tre”. Or, we can use other bash commands to edit our file name and then make an out variable. This becomes useful when running multiple steps and simply tacking on file extensions gets unwieldly. For example, we can do:

```{r analysis, results="markup"}
sample_name=`echo $input_aln | cut -d '.' -f 1`
out=`echo $sample_name`.tre
````

Here we are making a variable called “sample_name”, which is using _echo_ to print our input alignment name, then use _cut_ on our alignment name at the delimiter (-d) “.” and select the first field (-f). So, for “cluster1.pep.filtered.fa” the variable “sample_name” will keep only “cluster1”.

If we wanted our “sample_name” to be “cluster1.pep” we would use:

```{r analysis, results="markup"}
sample_name=`echo $input_aln | cut -d '.' -f 2`
````

If instead we want to remove the “.fa” (sometimes file names don’t have the same number of delimiters) we could use the _rev_ command like so: 

```{r analysis, results="markup"}
sample_name=`echo $input_aln | rev | cut -d '.' -f 2- | rev`
````

This would give us “cluster1.pep.filtered”. What _rev_ does is literally reverse the order of the text. So, for “cluster1.pep.filtered.fa”, _rev_ will print “af.deretlif.pep.1retsulc”. Using “-f 2-“ will _cut_ the “.fa” and the second _rev_ will return the file name to the original order. 

So, if we make our sample name by using: 

```{r analysis, results="markup"}
sample_name=`echo $input_aln | cut -d '.' -f 1`
````

Our output tree variable can be written like so:

```{r analysis, results="markup"}
out=`echo $sample_name`.tre
````

Which will name our tree “cluster1.tre”. Now we can add our output variable to our RAxML command:

```{r analysis, results="markup"}
raxmlHPC-PTHREADS-SSE3 -f a -x $RANDOM -N 100 -T 5 -p $RANDOM -m GTRGAMMA -s $input_aln -n $out
````

One other thing that might be helpful in the command above is the variable “$RANDOM”. It’s a random number generator, which is ideal for building trees or in MCMC analyses (rather than hard coding a number like ‘1234567’… that’s not quite random…)

So, our overall sbatch script might look something like this:

```{r analysis, results="markup"}
#!/bin/bash
#SBATCH --job-name=RAXML
#SBATCH --output=raxml_bash_%j.out
#SBATCH --error=raxml_bash_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sarah.carey@ufl.edu
#SBATCH --time=24:00:00
#SBATCH --nodes=1-1
#SBATCH --ntasks=5
#SBATCH --mem-per-cpu=1G
#SBATCH --qos=mcdaniel-b
#SBATCH --array 1-100%10

module load gcc/5.2.0
module load raxml/8.2.10

input_aln=`ls *.filtered.fa | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
sample_name=`echo $input_aln | cut -d '.' -f 1`
out=`echo $sample_name`.tre

raxmlHPC-PTHREADS-SSE3 -f a -x $RANDOM -N 100 -T 5 -p $RANDOM -s $input_aln -n $out -m GTRGAMMA
````




