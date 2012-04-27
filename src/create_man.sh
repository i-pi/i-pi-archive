help_dir=~/projects/NST_sim/wrap-pi/doc
src_dir=~/projects/NST_sim/wrap-pi/src
cd $src_dir
python help.py -x -o $help_dir/help
python help.py -l -n 1 -o $help_dir/simulation -r
python help.py -l -n 1 -o $help_dir/atoms -r -i atoms
python help.py -l -n 1 -o $help_dir/barostats -r -i barostats
python help.py -l -n 1 -o $help_dir/beads -r -i beads
python help.py -l -n 1 -o $help_dir/cell -r -i cell 
python help.py -l -n 1 -o $help_dir/forces -r -i forces
python help.py -l -n 1 -o $help_dir/ensembles -r -i ensembles
python help.py -l -n 1 -o $help_dir/interface -r -i interface
python help.py -l -n 1 -o $help_dir/thermostats -r -i thermostats
python help.py -l -n 1 -o $help_dir/random -r -i prng
