
1. run
python get_unident.py 
to find UFOs in the PASIPHAE survey zone
it will create unidentPASIPHAE.dat and unidentALL.dat with names and RA DEC coordinates
the first one contains only sources within PASIPHAE field |b|>40


2. run
python get_ext.py
it will create unidentPASIPHAEextinct.dat where the 4th column is the maximum of expected PD in the field of every UFOs


3. run
python simulate.py

USNO_ident.csv - crossmatch of USNOB with unidentified sources
