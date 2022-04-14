# Using TACC cluster

For those of you who have Knights Landing access through TACC
on the Stampede cluster, here are some simple instructions
to complete these exercises using that resource.



## Logging in to Stampede's KNL cluster

To access the KNL cluster you must log in to a
different login node, this one is

```
ssh username@login-knl1.stampede.tacc.utexas.edu
```

You will be prompted for both your normal password
and for your multifactor authorization passcode
which you have set up. 



## Submitting KNL jobs
Like most clusters, Stampede operates on a job management system which you must
utilize in order to use a node. An example job submission script is provided in this
directory titled "run.sh." Job submission scripts are simply shell scripts
which the job scheduler executes, but you supply extra parameters by
specially formatted comments.

Comments are supplied on the run.sh script describing what the parameters mean. For 
the most part you should not need to modify the script much except to change what
executable to run.

