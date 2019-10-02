### Useful linux/bash commands for running array jobs

Just as a note, Linux is an [operating system](https://en.wikipedia.org/wiki/Operating_system) and Bash is a [command-line interface](https://en.wikipedia.org/wiki/Command-line_interface)

UF HiPerGator has a lot of resources online for running jobs, for example for [submitting an array job](https://help.rc.ufl.edu/doc/SLURM_Job_Arrays). I want to add onto these resources by sharing some helpful commands/lines when I’ve run array jobs

## the commands
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

## running array job
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



