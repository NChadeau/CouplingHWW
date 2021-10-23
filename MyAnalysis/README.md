Command lines to compile the program:
$ mkdir build
$ cd build
$ cmake -C $ILCSOFT/ILCSoft.cmake .. 
$ make -j4 install

Then you have to go into the script/ folder. In the Python programs named: launchMyProcessor.py and launchMyProcessor_BG.py, you have to change the processID and the path where you find the .slcio files for the process you want to use. Each time it creat a ROOT file with the processID as name. You can then launch the program:
$ python3 launchMyProcessor.py or $ python3 launchMyProcessor_BG.py
/!\ /!\ /!\ Be Aware that only .slcio files with processID = 402007 or 402008 (n1n1h files) can be used for launchMyProcessor.py.

Once you have all .root files you need, you can use can go to the bin/ folder. Here you can do the study by running the executable files:
$ ./computeSignificance {electron polarisation} {positron polarisation} {integrated luminosity}
$ ./plots {electron polarisation} {positron polarisation} {integrated luminosity}
