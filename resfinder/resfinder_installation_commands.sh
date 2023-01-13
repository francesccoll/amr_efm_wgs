
#####################################################################################################################################
###					Installing ResFinder 4.2.3 - latest version					  ###
#####################################################################################################################################

# Installed again on 22/12/2022 to avoid this error: https://bitbucket.org/genomicepidemiology/resfinder/issues/74/run_resfinderpy-aborts-on-git-error
# https://bitbucket.org/genomicepidemiology/resfinder/src/master/
# It is no longer recommended to clone the ResFinder bitbucket repository unless you plan to do development work on ResFinder.
# Instead we recommend installing ResFinder using pip as described below.

python3 -m venv resfinder_env
source resfinder_env/bin/activate
pip install resfinder
deactivate

cd ~/software/kma
git clone https://bitbucket.org/genomicepidemiology/kma.git
cd kma
make
cp kma kma_index ~/bin/

base_dir=""; # local directory to run and save ResFinder results: to be edited

cd $base_dir
git clone https://bitbucket.org/genomicepidemiology/resfinder_db/
git clone https://bitbucket.org/genomicepidemiology/pointfinder_db/
git clone https://bitbucket.org/genomicepidemiology/disinfinder_db/

cd resfinder_db
python3 INSTALL.py ~/bin/kma_index
cd ../pointfinder_db
python3 INSTALL.py ~/bin/kma_index
cd ../disinfinder_db
python3 INSTALL.py ~/bin/kma_index


