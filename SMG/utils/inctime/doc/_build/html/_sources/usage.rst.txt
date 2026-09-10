Usage
=====

In this page is given an overview on how to use inctime.

Download
--------

The code is hosted in the Redmine portal at CPTEC. To checkout the code, use the following command:

.. code-block:: shell

  $ svn checkout https://svn.cptec.inpe.br/gdad/jgerd/tags/inctime

Compile
-------

All is needed is a fortran compiler. To compile the code, enter into the :code:`src` directory and type :code:`make`. The :code:`Makefile` will also compile the associated modules.

Once the compilation is done, an executable called :code:`inctime` is created.

Use
---

The way to use this program is through the command line as follows:

.. code-block:: shell

   $ ./inctime [yyyymmddhh, yyyymmdd] [<+,->nynmndnhnnns] [Form. output]

The inctime parameters are as follows:

  * [yyyymmddhh, yyyymmdd] 

    - Initial Time			

  * [<+,->nynmndnhnnns]    

    - ( -) calculate the passed date
    - ( +) calculate the future date (default)		
    - (ny) Number of year (default is 0)
    - (nm) Number of months (default is 0)
    - (nd) Number of days (default is 0)
    - (nh) Number of hours (default is 0)
    - (nn) Number of minutes (default is 0)
    - (ns) Number of seconds (default is 0)

  * [ Form. Output ]       

    - Format to output date. is a template Format
      The format descriptors are similar to those
      used in the GrADS:

      - "%y4"  substitute with a 4 digit year
      - "%y2"  a 2 digit year
      - "%m1"  a 1 or 2 digit month
      - "%m2"  a 2 digit month
      - "%mc"  a 3 letter month in lower cases
      - "%Mc"  a 3 letter month with a leading letter in upper case
      - "%MC"  a 3 letter month in upper cases
      - "%d1"  a 1 or 2 digit day
      - "%d2"  a 2 digit day
      - "%h1"  a 1 or 2 digit hour
      - "%h2"  a 2 digit hour
      - "%h3"  a 3 digit hour (?)
      - "%n2"  a 2 digit minute
      - "%e"   a string ensemble identify
      - "%jd"  a julian day without hours decimals
      - "%jdh" a julian day with hour decimals
      - "%jy"  a day of current year without hours decimals
      - "%jyh" a day of current year with hours decimals

.. seealso:: It is possible to use words to compose the output format.

More examples:

.. code-block:: shell

   $ ./inctime 2001091000 +1d %d2/%m2/%y4
   $ ./inctime 2001091000 +48h30n %h2Z%d2%MC%y4
   $ ./inctime 2001091000 -1h30n 3B42RT.%y4%m2%d2%h2.bin
   $ ./inctime 2001091000 -2h45n 3B42RT.%y4%m2%d2%h2.bin
   $ ./inctime 2001091000 -1y3m2d1h45n ANYTHING.%y4%m2%d2%h2.ANYTHING
