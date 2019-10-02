# Useful linux/bash commands for running array jobs

Just as a note, Linux is an [operating system](https://en.wikipedia.org/wiki/Operating_system) and Bash is a [command-line interface](https://en.wikipedia.org/wiki/Command-line_interface)

UF HiPerGator has a lot of resources online for running jobs, for example for [submitting an array job](https://help.rc.ufl.edu/doc/SLURM_Job_Arrays). I want to add onto these resources by sharing some helpful commands/lines when I’ve run array jobs

## The commands
To start, I have found these commands to be particularly helpful when running jobs on the cluster:

```{r analysis, results="markup"}# ls
ls
echo
cut
head
tail
rev
rm
````

_Man_ is also a good one to know; Less so within the scripts but it’s “an interface to the on-line reference manuals”. On the cluster, we can type “man” then any of these commands to get the description and the flag options (e.g. “man ls”)

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

If we want to run the same commands on the first two files (by default in [lexigraphical order](https://en.wikipedia.org/wiki/Lexicographical_order)) in the same directory, we add this line to the sbatch script: 

```{r analysis, results="markup"}
#SBATCH --array 1-2
````

If we want to run jobs 1-100 we’d use:

```{r analysis, results="markup"}
#SBATCH --array 1-100
````

If we only want to run jobs 1, 5, 17:

```{r analysis, results="markup"}
#SBATCH --array 1,5,17
````

If we want to run jobs on 1-100, but only submit 10 jobs at a time (useful when sharing an account with others):

```{r analysis, results="markup"}
#SBATCH --array 1-100%10
````

## One input, one output

Let’s say we want to build gene trees on 100 gene alignments using [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/). For the array line we'll use:

```{r analysis, results="markup"}
#SBATCH --array 1-100%10
````
So, we’re running 100 files and submitting them 10 at a time.

Here’s our RAxML command:

```{r analysis, results="markup"}
raxmlHPC-PTHREADS-SSE3 -f a -x $RANDOM -N 100 -T 5 -p $RANDOM -m GTRGAMMA -s INPUT_FILE_VARIABLE_HERE -n OUTPUT_FILE_VARIABLE_HERE
````

What we’ll need to do to is make variables for each of the input alignments and to output trees. For RAxML, the flags we need to make variables for are “-s” and “-n”, respectively. Now let's write something to call our particular input files. Let’s say our input alignment names are “cluster1.pep.filtered.fa, cluster2.pep.filtered.fa,… cluster100.pep.filtered.fa”. What we want to do is make a variable that will loop through each of our input alignments that we designate in our array. Here’s how we can make our input variable:

```{r analysis, results="markup"}
input_aln=`ls *.filtered.fa | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
````

In short, this variable captures the first 100 files (because that's what defined in our array line) that end in ".filtered.fa". More specifically, what we are first doing is using _ls_ to list everything in our directory that ends in “.filtered.fa”. We then pipe this to _head_ (the “|” is the pipe command) and we give the variable “$SLURM_ARRAY_TASK_ID” ([SLURM generates this one)](https://slurm.schedmd.com/documentation.html); it's the index value of your array) to the flag “-n” to list 1 to N files (index value 1 will list the first file, index value 2 will list the first two files, etc.) which we then pipe to _tail_ to get our last file in the list (for index value 1 we will end up with file 1, for index value 2 with the second file in the list, etc.). So, now we have our input variable for RAxML where each one of the alignment files will be passed to our RAxML commands one at a time.:

```{r analysis, results="markup"}
raxmlHPC-PTHREADS-SSE3 -f a -x $RANDOM -N 100 -T 5 -p $RANDOM -m GTRGAMMA -s $input_aln -n OUTPUT_FILE_VARIABLE_HERE
````

It’s important to remember when calling a defined variable to start with a “$” then give the variable name. 

Next we’ll want to make an output variable. We could very simply add the file extension “.tre” to our input file names, like so:

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

If instead we want to remove the “.fa” (sometimes file names don’t have the same number of delimiters) we could use _rev_ like so: 

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

## Two inputs, multiple outputs  

Now, let’s say we want to map paired-end genome resequence data using [HISAT2](https://ccb.jhu.edu/software/hisat2/index.shtml). But, because HISAT2 only outputs samples in the uncompressed sam format, we also want to run [Samtools](http://www.htslib.org/doc/samtools.html) to convert to a compressed bam (bams use a lot less storage than sams) and index the file to look at it in [IGV](https://software.broadinstitute.org/software/igv/). And we want to do this all the in the same script, so we don't have to submit multiple jobs. 

But, for paired-end data we now have two input files. So, how can we make our input variable(s)? Why math, of course! If the fastq files are in the same directory, usually read 1 will precede read 2, so the input will be in an order like:

1 Ala_M1_21.1.1_R1.fq.gz<br/>
2 Ala_M1_21.1.1_R2.fq.gz<br/>
3 Uco_F1_15.12.12_R1.fq.gz<br/>
4 Uco_F1_15.12.12_R2.fq.gz<br/>

So, first we make variables that will grab our paired index values (e.g. index values 1 and 2, and index values 3 and 4):

```{r analysis, results="markup"}
two=$(($SLURM_ARRAY_TASK_ID*2))
one=$(($two-1))
````

So, for example, for the first array index, "$SLURM_ARRAY_TASK_ID" will equal 1. Therein, "$two" will equal 2, because 1 * 2 = 2. Then "$one" will equal 1, because 2 - 1 = 1. For array index 2, "$SLURM_ARRAY_TASK_ID" will equal 2. Therein, "$two" will equal 4, because 2 * 2 = 4 and "$one" will equal 3, because 4 - 1 = 3.

Then we use these numbers to get read 1 and 2, where our fastq files end in ".fq.gz":
```{r analysis, results="markup"}
left_read=`ls *fq.gz | head -n $one | tail -n 1`
right_read=`ls *fq.gz | head -n $two | tail -n 1`
````
So, "$one" will always grab left reads (R1) and "$two" will always grab right reads. One thing to note, is even though we (in this example) have four input files, two of them go into each command (read 1 and read 2), so the array will only be for two runs:

```{r analysis, results="markup"}
#SBATCH --array 1-2
````

Then, like above, we can make a sample name:
```{r analysis, results="markup"}
sample_name=`echo $left_read | rev | cut -d '_' -f 3- | rev`
````

And make the output for samtools:
```{r analysis, results="markup"}
sam_out=`echo $sample_name`.sam
````

And now we can run HISAT2, where "-1" and "-2" take the inputs read 1 and read 2, respectively and "-S" is our sam output:

```{r analysis, results="markup"}
hisat2 -x Ceratodon_purpureus_mainGenome.fasta -1 $left_read -2 $right_read -p10 -S $sam_out
````

But, the great thing about bash variables is we can pass them to other commands. In this case, we can send our "$sam_out" to Samtools to convert it to bam. First we’ll make our output variable for our bam file (it's *strongly advised* to make a new file to write output to rather than overwritting an existing file):

```{r analysis, results="markup"}
bam_out=`echo $sample_name`.sorted.bam
````

Then we can run Samtools. Here we’ll sort the sam file and output it as bam. Then we’ll remove our rather large sam file using _rm_ and finally index the bam file:

```{r analysis, results="markup"}
samtools sort $sam_out -o $bam_out
rm $sam_out
samtools index $bam_out
```

All in all, our sbatch script might look something like this:

```{r analysis, results="markup"}
#!/bin/bash
#SBATCH --job-name=HISAT
#SBATCH --output=hisat_%j.out
#SBATCH --error=hisat_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sarah.carey@ufl.edu
#SBATCH --time=12:00:00
#SBATCH --nodes=1-1
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=1G
#SBATCH --qos=mcdaniel-b
#SBATCH --array 1-2

module load gcc/5.2.0
module load hisat2/2.1.0
module load samtools/1.8

two=$(($SLURM_ARRAY_TASK_ID*2))
one=$(($two-1))

left_read=`ls *fq.gz | head -n $one | tail -n 1`
right_read=`ls *fq.gz | head -n $two | tail -n 1`

sample_name=`echo $left_read | rev | cut -d '_' -f 3- | rev`
echo $sample_name
sam_out=`echo $sample_name`.sam
bam_out=`echo $sample_name`.sorted.bam

hisat2 -x Ceratodon_purpureus_mainGenome.fasta -1 $left_read -2 $right_read -p10 -S $sam_out

samtools sort $sam_out -o $bam_out
rm $sam_out
samtools index $bam_out
````


