beagle 4.0/4.1 source code

Fork of "Unofficial beagle 4.0/4.1 source code imported from http://faculty.washington.edu/browning/beagle/beagle.html"

Compile for ubuntu
------------------

Get java 8:

    sudo apt-add-repository ppa:webupd8team/java
    sudo apt-get update
    sudo apt-get install oracle-java8-installer

Setup to use the correct java version:

    sudo update-alternatives --config java

compile:

    ant


Run using the new jar:

    java -jar build/jar/beagle.jar

Updates:
--------

The following changes are made: A new input parameter _excludefromref_ to exclude samples from the reference vcf file. This has been done to use a sample filter on the reference vcf file, previously _excludesamples_ was used to filter the input gt, gl and/or gtgl files and the reference vcf file.

