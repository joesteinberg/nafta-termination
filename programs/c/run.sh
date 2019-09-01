#!/bin/bash

echo Baseline analysis...;
echo Baseline model...;
./bin/nafta > log.txt;

echo Dynamic ingredient sensitivity analyses...;
./bin/nafta -k > log_k.txt; # no capital adj cost
./bin/nafta -l > log_l.txt; # no labor adj cost
./bin/nafta -t > log_t.txt; # no trade adj cost
./bin/nafta -g > log_g.txt; # static exporting
./bin/nafta -j > log_j.txt; # no extensive margin
./bin/nafta -i > log_i.txt; # no IO
./bin/nafta -s > log_s.txt; # symmetric elasticities
./bin/nafta -e > log_e.txt; # cobb-douglas production
./bin/nafta -n > log_n.txt; # fixed trade balances

echo Alternative scenarios...;
./bin/nafta -d > log_d.txt; # USMCA
./bin/nafta -r > log_r.txt; # stricter domestic content requirements
./bin/nafta -q > log_q.txt; # us-canada FTA
./bin/nafta -c > log_c.txt; # can-mex FTA
./bin/nafta -u > log_u.txt; # higher US tariffs unilaterally
./bin/nafta -o > log_o.txt; # old MFN tariffs
./bin/nafta -e -o -n > log_o_e_n.txt # combo of fixed TB, old MFN, cobb-douglas prod
