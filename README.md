# Data

```
https://www.ebi.ac.uk/genomes/phage.html
```


# NCBI toolkit

NCBI = ncbi_cxx--22_0_0

Configure toolkit to build only what is required for BLAST project.
```
NCBI/configure --with-projects=NCBI/scripts/projects/blast/project.lst
```

Compile toolkit for configured project.
```
cd NCBI/GCC700-DebugMT64/build && make all_p
```

Get source code for BLAST app
```
mkdir foo
NCBI/scripts/common/new_project.sh foo app/blast NCBI/GCC700-DebugMT64/build/
cd foo
```

At this point you need to configure the makefiles so that they point to to the built NCBI toolkit.

Compile your BLAST app.
```
make -f Makefile.foo_app
./foo -help
```


# HPC

```
http://eresearch-info.qut.edu.au/Guides/HPC/Connecting/

/work/cripsr_bio
```






# Archived code
Archived code at: https://github.com/benjaminsaljooghi/CRISPR-archive


