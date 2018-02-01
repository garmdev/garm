

---[[ Copyright




                              +----------------+
                              | zt version 1.1 |
                              +----------------+

                       copyright 2001-2007 Eric Bonnet 
                       mailto:eric.bonnet@psb.ugent.be






  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


---[[ Files

  zt.pdf        zt description & user guide
  zt            solaris binary executable file, compiled with GNU gcc.
  readme.txt    this file.
  gene.txt      test file.
  env.txt       test file.
  geo.txt       test file.
  assoc.txt     test file.
  gond.txt      test file.
  pres.txt      test file.
  zt.c          source file. 
  rr.c          source file. 
  rr.h          source file (header).
  gpl.txt       GNU General Public License


---[[ Installation 

  First create a directory where you want to store the soft.
  Then extract the zip archive in this directory.


---[[ Usage

  zt is a command line tool, so you'll have to open a dos box first. You can
  do that with the command Start>Run and then you type "cmd". A dos box should
  appear. Then you'll have to move to the directory where you've extracted all
  the files, using the command cd, for example "cd c:\zt". Then you can
  execute the program zt.exe with the different options:
  
  Simple Mantel test:      
  zt.exe -s file1 file2 number_of_randomizations

  Partial Mantel test:     
  zt.exe -p file1 file2 file3 number_of_randomizations

  Options:
  -l display licence.
  -r partial Mantel test with raw option.
  -e force exact permutations procedure.
  -h display this help.

  Number of randomizations are allowed between 99 and 999999999.

  Examples:

  simple Mantel test:
  zt.exe -s assoc.txt gond.txt 10000

  partial mantel test:
  zt.exe -p gene.txt env.txt geo.txt 10000

  For further details, see zt.pdf


---[[ Data files format

  Input files must be flat ASCII files.
  Matrix are half-matrix (lower part of the square matrix),
  without (!) diagonal terms.
  The first number is the size of the complete matrix.
  For Example if you have the following half-matrix data:

0.2 
0.3 0.4 

  The size of this matrix is 3 x 3.
  So, the correct data file for zt will be:

3
0.2 
0.3 0.4 

  Be careful, floating point must be dot(.) not comma(,).
  Try to avoid trailing spaces as well. Using the Notepad is recommended to
  create input files.


---[[ Compilation

  Compilers are not installed by default on Windows, so first you will have to
  install one, like for example cygwin (http://www.cygwin.com).

  There are only three source files, so compilation should be straightforward.

  For example with gcc:

  gcc -o zt zt.c rr.c -lm


---[[ Fee

  There's no fee for this soft. This software is free (as in free beer) but is 
  also emailware: if you like it (or not) please send me a small email !

---[[ Thanks

  Since 2002, many people dropped me a line, with many useful comments.
  Thanks to you all !!!

  Leena Uimaniemi
  Junfeng Pang
  Joerns Fickel
  Justin Quirouette
  Hillary Mrosso
  Dr. W. M. Hartmann
  Santiago A. Catalano
  Dan Bolnick
  Geoffrey Bell
  Wilbert Heeringa
  Nidal Odat
  Mark E. Berres
  Joah Madden
  Faith M. Walker
  Christophe Bouget
  D. Ashley Robinson
  Jeffrey S. McKinnon
  Juan Núñez-Farfán
  Karen Carney
  Masahiro Nakaoka
  Jason B. MacKenzie
  Fernanda Michalski
  Dr. Kristina Sefc
  Miguel Olvera Vargas
  Lindie Janse van Rensburg
  Sylvain Gerber
  Susanne Hauswaldt
  Sergios-Orestis Kolokotronis
  Richard Case
  xianyun zheng
  M.A. Spoljaric
  Sylvie Cousin
  Eliot McIntire
  Mahmood Sasa
  Stephen Gregory
  Tawfiq Froukh
  Andrew J. Bohonak
  Aaron Wagner
  Sarah Fishman
  Raymond Wan
