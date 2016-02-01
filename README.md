Description
------------------
Software suite that applies hidden Markov models (HMM) to assembly and NGS read assessment and repair.

Installation
------------------

1. Clone the repository:

        git clone git@github.com:jordiabante/ngs_hmm.git 

2. Add the bin folder in the repository to your path variable by adding the following line to the `.bashrc` (`.bash_profile` in Mac) file:

        export PATH="path/to/bin:$PATH"

3. If the machine running the code is Machintosh install coreutils (sudo port install coreutils) as well as: sed, gzip, grep, getopt, mysql, openssl, readline, boost, zlib. Add the following line to the `.bashrc` (`.bash_profile` in Mac) file:

        export PATH="/opt/local/bin:$PATH"
        export PATH="/opt/local/libexec/gnubin:$PATH"

Dependancies
-----------------
1. perl{Algorithm::Combinatorics}

Usage
-----------------

1. Every folder (`script_name`) will contain the following files and directories:
  1. `sript_name.sh`: bash script
  2. `main/`
    * `main.sh`: bash script to try the command with toy example.
    * Files required for the toy example.

