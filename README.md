#### About Gcoder

Gcoder is a simple python script that codes gaps in a multiple sequence alignment according to the method of Simmons and Ochoterena (2000). It works exclusively through the command line and requires the ElParsito.py module, which can be found in the ElConcatenero package (https://github.com/ODiogoSilva/ElConcatenero.git).

It currently supports as input file formats:

- FASTA
- Phylip
- Nexus

The current supported output file formats are Nexus only. This is because this script was made with MrBayes analyses in mind. Therefore, it already makes the appropriate changes in the Nexus header so that MrBayes can recognize the binary state block as restriction site data).

There is no need for instalation. The only requirement is that the ElParsito.py module must be on the same directory as Gcoder.py (or you can use any other way to let the main script know where the module is). I do recommend, to make it easier to call the program, that you add it to your $PATH variable. 

Finally, please note that Gcoder.py is far from immune to bugs and crashes, and I'll be happy to know about them through my e-mail (o.diogosilva@gmail.com) so that I can fix them. Any suggestions or comments are also welcome.

#### Options

Gcoder.py has the following options (which can also be consulted by typing "Gcoder.py -h" in the command line):

  -h, --help            **show this help message and exit**
  
  -if *{fasta,nexus,phylip}* **Format of the input file (default is 'nexus')**
  
  -in *INFILE*            **Input file. Only one allowed!**
  
  -o *OUTFILE*            **Name of the output file (default is 'Outfile')**
								
#####Note:

The order of the options does not matter.
								
#### Usage
