#!/usr/bin/bash
#installer script for MG-TK
#stay in helpers/install/ dir while executing

#some basic housekeeping
#set -e
ulimit -c 0;
#set

echo "MG-TK installer script"
echo "This script will install several conda environments with most dependencies for MG-TK. It will use micromamba to run the installations, please ensure micromamba is installed natively for your account.";
echo "You can (and probably should) rerun this script every time you pull a major MF update. Changes to dependencies will be automatically updated and rerunning the script will be significantly faster than running it the first time.";

if ! command -v micromamba &> /dev/null
then
    echo "micromamba could not be found"
	echo "Make sure micromamba is in your \$PATH"
	echo "Aborting"
    exit
fi

if [ -z "${MAMBA_EXE}" ] ; then
MAMBA_E=$MAMBA_EXE
else
MAMBA_E=micromamba
fi
#which micromamba


eval "$($MAMBA_E shell hook --shell=bash)"

find_in_mamba_env(){
	$MAMBA_E env list | grep "${@}" >/dev/null 2>/dev/null
}

find_in_bashrc(){
	grep "${@}" ~/.bashrc >/dev/null 2>/dev/null
}


echo "Using micromamba version:"
$MAMBA_E --version
#exit


SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
MFdir=$(realpath -s $SCRIPT_DIR/../..)

#remove stone that declares programs are checked by MG-TK
rm -f $MFdir/helpers/install/progsChecked.sto
touch $MFdir/helpers/install/runningInstall.sto

mkdir -p $MFdir/gits/ 
INSTdir=$MFdir/helpers/install/
DBdir=$MFdir/data/DBs/

#MFLRDir	$MF3DIR
if [ ! -f $MFdir/config.txt ] ; then
	cp -f $MFdir/Mods/config.old $MFdir/Mods/MATAFILERcfg.txt
	#should be defaulted to $MF3DIR now
	#sed -i "s+MFLRDir.*+$\t$MFdir+" $MFdir/Mods/MATAFILERcfg.txt
	ln -s $MFdir/Mods/MATAFILERcfg.txt $MFdir/config.txt
	echo "Rewrote config.txt. Please modify as needed to local paths"
fi


if ! find_in_bashrc "##------------> MG-TK ADDED" ; then
	printf "\n\n##------------> MG-TK ADDED <----------##\nexport MGTKDIR=$MFdir/\nexport PERL5LIB=\"\$PERL5LIB:$MFdir/\"\n##------------> MG-TK ADDED <----------##\n\n" >> ~/.bashrc
	echo "Added MG-TK modules to .bashrc"
fi


#exit;


# For all micromamba installs, we use --channel-priority 1. This sets to 
# "flexible" priority - from conda docs "the solver may reach into lower 
# priority channels to fulfill dependencies, rather than raising an 
# unsatisfiable error". This should help avoid a lot of dependency problems
# during install, as micromamba defaults to "strict".
if ! find_in_mamba_env "MGTK\s" ; then
	echo "Creating base MGTK conda environment.. This might take awhile"
	#first install spades that seems to require a lot of mem..
	#$MAMBA_E create --channel-priority 1 -q -y -n MFF spades
#just in case it crashes, this often recovers it..
	$MAMBA_E create --channel-priority 1 -q -y -f $INSTdir/MG-TK.yml #-q -y
	$MAMBA_E activate MGTK
	#echo "Installing R packages in MGTK environment";	Rscript $INSTdir/reqPackages.R
	#pip install biopython
else 
	echo "Updating base MGTK conda environment.. Please be patient"
					#	$MAMBA_E activate MGTK
	$MAMBA_E update --channel-priority 1 -q -y -f $INSTdir/MG-TK.yml #-q -y
	
	#echo "Updating R packages in MGTK environment"
	#{ Rscript $INSTdir/reqPackages.R
	#} || { 		echo "Rscript install failed.. trying direct excecution";		$INSTdir/./reqPackages.R;	}
fi

#exit

$MAMBA_E activate MGTK

if command -v foldseek &> /dev/null ; then
	#prepare foldseek search
	if [ ! -d $DBdir/PtostT5_W ];then
		foldseek databases ProstT5 $DBdir/PtostT5_W $DBdir/tmp;
	fi

fi

echo ""
echo "Installing/updating further dependencies in additional conda environments.."
echo ""
echo "" 


#if ! find_in_mamba_env "checkm2" ; then
#	#git clone https://github.com/chklovski/CheckM2.git $MFdir/gits/checkm2/
#	echo "Installing checkm2"
#	$MAMBA_E create -y -q -f $INSTdir/checkm2.yml -n checkm2
#	$MAMBA_E activate checkm2
#	pip3 install --upgrade pip

#	pip install CheckM2 packaging
#	$MAMBA_E activate MFF
#	echo "checkm2 installed. you can verify this environment by running \"micromamba activate checkm2\" and \"checkm2\""
#else
#	$MAMBA_E activate checkm2
#	if [ ! command -v checkm2 > /dev/null 2>&1 ]; then
#		echo "Could not find checkm2. Please install via \"micromamba activate checkm2;pip install CheckM2 packaging\" and restart installer"
#		exit 5
#	fi
#	$MAMBA_E activate MFF
#fi

if [ ! -f "$SCRIPT_DIR/../../gits/XGTDB/extract_gtdb_mg.py" ]; then
	echo "Installing extractGTDB into $MFdir/gits/XGTDB/"
	git clone https://github.com/4less/extract_gtdb_mg.git $MFdir/gits/XGTDB/
fi

#additional dependencies not in the main yml..
if ! find_in_mamba_env "MGTKgtdbtk" ; then
	echo "Installing MGTKgtdbtk environment"
	$MAMBA_E create --channel-priority 1 -q -y -f $INSTdir/GTDBTK.yml 
else 
	echo "Updating MGTKgtdbtk environment"
	$MAMBA_E update --channel-priority 1 -q -y -f $INSTdir/GTDBTK.yml 
fi

if ! find_in_mamba_env "MGTKbinners" ; then
	echo "Installing MGTKbinners environment"
	$MAMBA_E create --channel-priority 1 -q -y -f $INSTdir/Binners.yml
else 
	echo "Updating MGTKbinners environment"
	$MAMBA_E update --channel-priority 1 -q -y -f $INSTdir/Binners.yml
fi





if ! find_in_mamba_env "MGTKcheckm2" ; then
	echo "Installing MGTKcheckm2 environment"
	$MAMBA_E create --channel-priority 1 -q -y -f $INSTdir/checkm2.yml 
	$MAMBA_E activate MGTKcheckm2
	checkm2 database --download --path $MFdir/DBs/
	$MAMBA_E deactivate
else 
	echo "Updating checkm2 environment"
	$MAMBA_E update -y -q -f $INSTdir/checkm2.yml 
	if [ ! -d $DBdir/CM2/ ]; then
		$MAMBA_E activate MGTKcheckm2
		checkm2 database --download --path $DBdir/CM2/
		$MAMBA_E deactivate
	fi
	
fi


#if ! find_in_mamba_env "comeBin" ; then
#	echo "Installing comeBin environment"
#	$MAMBA_E create -y -q -f $INSTdir/comeBin.yml
#else 
#	echo "Updating gtdbtk environment"
#	$MAMBA_E update -y -q -f $INSTdir/comeBin.yml 
#fi

#additional dependencies not in the main yml..
#if ! find_in_mamba_env "MGTKmetaMDBG" ; then
#	echo "Installing MGTKmetaMDBG environment"
#	$MAMBA_E create --channel-priority 0 -q -y -f $INSTdir/metaMDBG.yml
#else 
#	echo "Updating MGTKmetaMDBG environment"
#	$MAMBA_E activate MGTKmetaMDBG
#	$MAMBA_E update --channel-priority 0 -q -y -f $INSTdir/metaMDBG.yml 
#	$MAMBA_E deactivate 
#fi

#if ! find_in_mamba_env "motus" ; then
#	echo "Installing motus environment"
#	$MAMBA_E create -y -q -f $INSTdir/motus.yml
#else 
#	echo "Updating motus environment"
#	$MAMBA_E update -y -q -f $INSTdir/motus.yml 
#fi

if ! find_in_mamba_env "MGTKphylo" ; then
	echo "Installing MGTKphylo environment"
	$MAMBA_E create --channel-priority 1 -q -y -f $INSTdir/phylo.yml 
else 
	echo "Updating MGTKphylo environment"
	$MAMBA_E update --channel-priority 1 -q -y -f $INSTdir/phylo.yml 
fi

if ! find_in_mamba_env "MGTKwhokar" ; then
	echo "Installing MGTKwhokar environment"
	$MAMBA_E create --channel-priority 1 -q -y -f $INSTdir/whokaryote.yml 
else 
	echo "Updating MGTKwhokar environment"
	$MAMBA_E update --channel-priority 1 -q -y -f $INSTdir/whokaryote.yml 
fi




#if ! find_in_mamba_env "Rbase" ; then
#	$MAMBA_E create -y -f $INSTdir/Rbase.yml
#	$MAMBA_E activate Rbase
#	Rscript reqPackages.R
#else 
#	$MAMBA_E update -y -f $INSTdir/Rbase.yml --prune --allow-uninstall
#fi

#later toadd..
#git clone https://github.com/GaetanBenoitDev/metaMDBG.git;mima create -y -f conda_env.yml;activate metaMDBG

rm -f $MFdir/helpers/install/runningInstall.sto


echo "Finished MG-TK install"
echo ""
echo "To run MG-TK, make sure you are in the MGTK environment (micromamba activate MGTK)."
echo "You can rerun the installer.sh anytime, to ensure package were installed or are being updated."
echo "Run \"MG-TK.pl -checkInstall\" to ensure that the installation was successful."
exit 
