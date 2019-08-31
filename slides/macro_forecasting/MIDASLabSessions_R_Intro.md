MIDAS Lab Sessions: SETUP
========================================================
autosize: true

    
Setup required package for running examples
========================================================
First you need to install devtools to be able to source GitHub repository files

```r
#install.packages("devtools")
require(devtools)
```
Install required package


```r
install_github("jstriaukas/MIDASLec")
```
Check if it loads up


```r
require(MIDASLec)
```

Alternative
========================================================

Go to

<center>

**https://github.com/jstriaukas/MIDASLec**

</center>




Download all files to your PC (press Clone or download -> Download ZIP)

Once you have files on your PC, store it on some path

Folder *example_code* has files for running all example code we will cover

Before every session run:



```r
path <- "your.path.to.folder"
source(path)
```


