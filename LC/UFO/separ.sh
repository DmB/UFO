#!/bin/bash


gawk '{if(($7=="Ignored")||($1=="#")) print $0}' RBPLJ0221+2518.dat > RBPLJ0221+2518_gamma_fl.dat
gawk '{if(($7!="Ignored")&&($7!="ERROR!")&&($7!="Error")) print $0}' RBPLJ0221+2518.dat > RBPLJ0221+2518_gamma_up.dat