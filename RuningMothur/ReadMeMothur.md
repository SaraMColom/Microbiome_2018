# ReadMe for runing mothur
This subfolder was created by Sara Colom to run mothur for sequence processing and OTU construction via batch script on the University of Michigan Great Lakes-ARC-TS (Slurm cluster, campus-wide computing cluster).
See website: https://arc-ts.umich.edu/greatlakes/user-guide/ for detailed instruction.

# General steps:
1.Sign into the University of Michigan Great Lakes system, replace InsertUsernameHere with your username.

`ssh -l InsertUsernameHere greatlakes.arc-ts.umich.edu`

2. Navigate to your scratch file replace InsertUsernameHere with your username.
 
`cd /scratch/lsa_root/lsa/InsertUsernameHere/`

3. Create slurm--see link here: https://arc-ts.umich.edu/greatlakes/slurm-user-guide/

4. Load mothur version 1.43

`module load Bioinformatics`
`module load mothur/1.43.0`

5. Unzip fastaQ files, replace 'FileNameHere' with your file name(s)

`gunzip FileNameHere`

6. Convert fastaQ to fasta files

5. Run job, replace JobName with your job's name.

`sbatch JobName`

## Important note:
Try runing a subset of your data locally first before sending a job. 
You can avoid losing time/money if troubleshooting is required.
Make sure your fasta files are not too long, too many '-' or other seperators can lead to bugs in mothur.
File names should be conventional for MiSeq sequencing output. 
