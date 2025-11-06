
## FAQ

<details>
  <summary>Expand frequently asked questions section</summary>

This section lists a number of typically occurring problems that are usually not addressable by programming/bug fixing. Please look here first if an error you encountered is already listed.

#### .mfc files

- **Q**: What is a \*.mfc file?
**A**: mfc stands for 'MATAFILER control file' and refers to the original name of MG-TK: MATAFILER. .mfc files are used to store the call to MG-TK, but could as well be named \*.sh or \*.something

#### Assembly too large??

- **Q**: I have an enormous assembly (e.g. several soils together) and suddenly all kind of things break.
**A**: Once your assembly has too many bp or entries, mappers and the sam format can reach it limits (usually defined by INT32). The first thing to try is to use bowtie2 (-mapper 1) as mapper and set the large mem flag (-mapperLargeRef 1). You can still run into trouble with the sam format not supporting that many different fasta entries. Here, you can try to reduce the number of actual entries in the fasta file, by increasing -assemblyScaffMinSize from the default 500, to e.g. 1000. That means that all contigs <1000bp would not be used in the mapping, gene catalog, binning etc.

#### Setting environmental variables

**Problem:** If you have issue to source or define variables from/inside of MFF micromamba environment.

**Solution:** Adding 2 files into your MFF micromamba environment will help you resolve the issue (talk to Klara or Jogi).

#### sbatch error: Batch job submission failed: Job dependency problem

**Problem:** Usually some jobs run but some jobs are hanging out in queue and never finish. In that case we can use `scontrol show job [ID]` to see which are the jobs and dependencies. 

**Solution:** If the dependencies are never fulfilled then we can delete all those jobs using `scancel`, after that rerun MG-TK. What MG-TK does, is to pick up where it left of - the pipeline is designed in a way that it creates `stone` files for finished processes so the pipeline knows where to continue. Sometimes files are also still in scratch dir and not copied to final dir. You just need to re-run command.

You can also delete **all** jobs where the dependency failed, saving you a lot of time (command: `squeue -u $USER | grep 'dencyNev' | cut -f11 -d' ' | xargs  -t -i scancel {}`) or ask MG-TK to do this in case the max number of jobs is reached via the flag `-killDepNever 1`.

#### local tempdir on node not writable

**Problem:** Sometimes you will encounter an error where a job starts for ~1 sec on a node and immediately aborts, the error log (.etxt) showing something about not being able to create/write to a local dir (e.g. `/nbi/local/tmp/12312421/MG-TK`). This is usually the local SSD space not being available (for various reasons not related to MG-TK). Such nodes act as kind of 'honeypots', accepting a lot of jobs and killing them immediately. This can seriously harm your performance in getting jobs done.

**Solution:**  Therefore it is important to a) let you local sysadming know that the SSD is no longer available on said node (node name is always printed as first line in the MG-TK .otxt logs for a job) and b) you can exclude this node from MG-TK using the `-excludeNodes [nodename1,nodename2]` MG-TK flag.

#### Recursion error while running contig stats

**Problem:** While calculating contig stats for large samples, you may encounter `RecursionError: maximum recursion depth exceeded in comparison` error. 

**Solution:** Increase the recursion depth in `extract_gtdb_mg.py` by using `sys.setrecursionlimit()` function, e.g. include `sys.setrecursionlimit(1500)` at the beginning of the python script. 

#### Unusually high quality values

**Problem:** Some samples show unusually high quality of over 80

**Solution:** Some sequencing systems omit detailed quality values over a certain threshold. In this case, the quality values given are among the highest possible (>80) and do not represent the actual quality of sequences. To solve this quality values have to be calculated manually.

#### GeneCat stops without producing MGS and no error messages appear

**Problem:** GeneCat.pl stops without error messages but no MGS are produced. Close examination shows missing genes among some bins, but these genes are present in the assembly.

**Solution:** This error can occur when previous GeneCat runs stop unexpectedly or fail due to previous issues with the assembly. In this case, files could be created but not completed. The MGS's are then not able to finish due to the expectation of some genes being present but the previous run had stopped before these could be written. It is best to make sure all assemblies are completed and then start a fresh GeneCat run.

#### Automatic installation with the installer script does not finish due to several dependency issues

**Problem:** Parts of the installer using micromamba do not finish due to dependency issues

**Solution:** This can be resolved stepwise. First, make sure the environment where the problem occurs is created with the right name in micromamba. Then install packages (or dependencies) that have issues manually with micromamba. Try conda_forge first, then bioconda for the -c parameter. Restart the installer and note down any further issues. When a problem occurs with a package that is already installed, it can help to remove it and then reinstall it manually. If manual installation does not resolve conflicts, remove version numbers from the yml file of problematic packages and start this process again. It is important to note this somehow, in order to troubleshoot later on if any incompatible version of a package was installed this way.

#### I want to get MAGs of a specific MGS?

**Problem:** I want to extract all MAGs from an MG-TK output that are of a particular MGS (metagenomic species). 

**Solution:** We have an accociated repository named [MAGRec-TK](https://github.com/huminfo8/magrec-tk/tree/main) which you can install and run on the MG-TK output. 

#### Unrecognised CONDcmd (or other similar error during installation)

**Problem:** On startup there are errors, originating in the MG-TK config file (for example, CONDcmd not recognised). 

**Solution:** MG-TK map and config files are very strictly tab delimited. The whitespace delimiting fields (in the mapping file, or the MG-TK config file) must be a tab, any other kind of whitespace character will cause an error. In vim you can view whitespace using the command `:set list`, and tab appears as "^I". 

</details>
