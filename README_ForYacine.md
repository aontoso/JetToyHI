# Setting up software and run jet analysis

```sh
cd <dirInWhichYouWantToInstall>
```

## Install ROOT
The easiast is to just grep a precompiled version from the root website (take ROOT6)
* https://root.cern.ch/content/release-61404

## Install PYTHIA8
```sh
wget http://home.thep.lu.se/~torbjorn/pythia8/pythia8235.tgz
tar xvfz pythia8235.tgz
cd pythia8235
./configure
make
PYTHIA=$PWD
cd ..
```

## Install fastjet

```sh
curl -O http://fastjet.fr/repo/fastjet-3.3.0.tar.gz 
tar zxvf fastjet-3.3.0.tar.gz
cd fastjet-3.3.0/

./configure --prefix=$PWD/../fastjet-install
make
make check
make install
FASTJET=$PWD/../fastjet-install
cd ..

export FJ_CONTRIB_VER=1.026 
curl -Lo source.tar.gz http://fastjet.hepforge.org/contrib/downloads/fjcontrib-"$FJ_CONTRIB_VER".tar.gz
tar xzf source.tar.gz
cd fjcontrib-"$FJ_CONTRIB_VER"
./configure --fastjet-config=$FASTJET/bin/fastjet-config --prefix=`$FASTJET/bin/fastjet-config --prefix`
make 
make install 
make fragile-shared #make shared library
make fragile-shared-install
cd ..
```

## Jet workshop software
```sh
git clone https://github.com/mverwe/JetToyHI.git
cd JetToyHI
echo `$FASTJET/bin/fastjet-config --prefix` > .fastjet
echo `$PYTHIA` > .pythia8
git pull --rebase origin sharedLayers

```

```sh
cd PU14
echo `$FASTJET/bin/fastjet-config --prefix` > .fastjet
./mkmk
make
cd ..

scripts/mkcxx.pl -f -s -1 -r -8 '-IPU14' -l '-LPU14 -lPU14 -lz'
make
```

Now you are done installing software. Let's generate 10 pythia events and run a simple jet analysis.
```sh
./runCreatePythiaEvents -nev 10 -pthat 120 -tune 14
./runSimpleJetAnalysis -hard PythiaEventsTune14PtHat120.pu14  -nev 10
```

You will have produced a root file with a tree. In this tree properties of jets are stored in std::vector format. To check what is inside do:
```
root -l
TBrowser b
```
Click on `jetTree` and play around.


## Contribute
* If you want to contribute to this code you need to have a github account. Go here to do so: https://github.com/join.
* Fork the original repository. Go to: https://github.com/JetQuenchingTools/JetToyHI and click 'Fork' in the upper right corner.
* Instead of cloning the original repository as shown above, clone your own.
* After committing your changes to your own branch, push them to your own fork. Don't know how to do this, ask your colleages or use google which might bring you here https://services.github.com/on-demand/downloads/github-git-cheat-sheet/
* Do a pull request once you have finished your developements.

## Samples
Event samples can be found in the jet quenching CERNBOX:
* From lxplus: /eos/project/j/jetquenching/JetWorkshop2017/samples
* Webbrowser: https://cernbox.cern.ch/index.php/s/JKvuTJMD0qCz6dJ
* Mount eos on a laptop or local desktop: https://cern.service-now.com/service-portal/article.do?n=KB0003493 

You will find samples from various event generators. For underlying event we have: 'thermal' which is independent particle production using a Boltzmann distribution with a fixed multiplicity and mean p<sub>T</sub> (indicated in the file names). For the hard signal we have PYTHIA8 and JEWEL events with various p<sub>T,hat</sub> settings.

More details about the available samples can be found here: https://twiki.cern.ch/twiki/bin/view/JetQuenchingTools/PU14Samples
