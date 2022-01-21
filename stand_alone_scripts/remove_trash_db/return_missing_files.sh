#!/bin/bash

################################################################################
#                              Identify missing db files                       #
#                                                                              #
# Identifies xyz files you wanted to simulate                                  #
# Identifies yaml files that have been generated                               #
# Compares them                                                                #
# Writes names of xyz files that were not simulated                            #
#                                                                              #
################################################################################
################################################################################
################################################################################
#                                                                              #
#  Copyright (C) 2022, 2022 Artem Fediai                                       #
#  artem.fediai@both.org                                                       #
#                                                                              #
#  This program is free software; you can redistribute it and/or modify        #
#  it under the terms of the GNU General Public License as published by        #
#  the Free Software Foundation; either version 2 of the License, or           #
#  (at your option) any later version.                                         #
#                                                                              #
#  This program is distributed in the hope that it will be usefult 
#  WITHOUT ANY WARRANTY; without even the implied warranty of                  #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               #
#  GNU General Public License for more details.                                #
#                                                                              #
#  You should have received a copy of the GNU General Public License           #
#  along with this program; if not, write to the Free Software                 #
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA   #
#  or call you mom                                                             #
################################################################################
################################################################################
################################################################################



 Help()
{
   # Display Help
   echo "."
   echo
   echo "Syntax: scriptTemplate [-h|V]"
   echo "options:"
   echo "h     Print this Help."
   echo "V     Print software version and exit."
   echo
}


while getopts ":hV" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      V) # display Help
         echo "Script writes names of xyz files that were not simulated"
         echo "version 0.0"
         exit;;
#     \?) # incorrect option
#         echo "Error: Invalid option"
#         exit;;
   esac
done





###############################################################################
###############################################################################
###############################################################################




printf '\nI will try to return missing files, i.e. the files which did not produced the simulation results\n\n'

DB_FOLDER_PATH='.'  #<-- specify
DATASET_PATH='../../dataset_qm17' #<-- specify

DATASET_FILES='db_folder_content_no_fixes.txt'   # optionally specify
DB_FILES='dataset_folder_content_no_fixes.txt'   # --||--
NOT_PRESENTING_FILES='not_presenting_files.txt'  # --||--

echo The path to the dataset: `realpath $DATASET_PATH`
echo The path to db folder: `realpath $DB_FOLDER_PATH`
#find $DATASET_PATH/  -printf "%f\n"

# 1. db folder (*yaml files)
# returns only the name of the molecule -->
ls $DB_FOLDER_PATH/DB_*.yaml |  # valid files start with DB_ and end with *.yaml
sed 's/'"$DB_FOLDER_PATH"'//' | # exclude also the path
sed 's/\///' |  # and the slash that will remain
sed 's/DB_//' | # exclude DB_
sed 's/.yaml//' > db_folder_content_no_fixes.txt # db_folder_content_no_fixes.txt

printf '\nsaved db files\n'


# 2. dataset folder (*.xyz files)
# returns only the name of the molecule -->
ls $DATASET_PATH/*.xyz |  # valid files are all ending with *.xyz
sed 's~'"$DATASET_PATH"'~~'|  # exclude also the path
sed 's~\/~~' |   # and the slash that will remain
sed 's~.xyz~~' > dataset_folder_content_no_fixes.txt # db_folder_content_no_fixes.txt

printf '\nsaved dataset files\n'


printf '\nSaving not presenting filenames into '$NOT_PRESENTING_FILES' \n'
# get the numebers that are not in the database (db)
comm -13  <(cat $DATASET_FILES | sort) <(cat $DB_FILES | sort) | sed 's/$/.xyz/'  > $NOT_PRESENTING_FILES



