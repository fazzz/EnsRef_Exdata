# EnsRef_Exdata

###Installation
``
./configure CPPFLAGS="-I_LBFGSINSTALLDIR_/include" LDFLAGS="-L_LBFGSINSTALLDIR_/lib"
make
make install
``

###How to use
``
ensResEx _experimentaldatafilename_ _simulationdatafilename_ _weightfilename_ _observalfilename_
``

###Options
[--Del]
[--Expm]
[--Expm2]
[--h] help `show the help message`
