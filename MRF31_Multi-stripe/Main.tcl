
##############################################################################################################
# Reagan Chandramohan                                                                                        #
# John A. Blume Earthquake Engineering Center                                                                #
# Stanford University                                                                                        #
# Last edited: 02-Jun-2015

# Edited by Kuanshi Zhong
# Last edited: 03-Feb-2017

# Edited by wyy
# Last edited: 04-April-2018
##############################################################################################################

# Run a stripe analysis using multiple cores.
# For best results, run with as many processors as there are analyses to be run.

##############################################################################################################

# Set up
wipe all;
setMaxOpenFiles 2048;
# Output main directory
set outmain [file join [pwd] Result]; # change here
# set outmain trial;

# Some constants
set g 386.1;
set col_drift 0.10;
set num_stories 35; # change here
# set num_piers 4; # change here
set num_bays 5; # change here
set TmaxAnalysis [expr 5*60]; # five minute maximum
set elemConvTol 1e-6; # element convergence tolerance

# Source procedures
set userdirectory [pwd];

source $userdirectory/PanelZone2D2.tcl;
source $userdirectory/preNRBeam.tcl;
source $userdirectory/Wsection.tcl;
source $userdirectory/max_drift_check.tcl;

# Initialize the total number of processors and the id of this processor
set numprocs [getNP]
set this_procid [getPID]


# Initialize the list of ground motion folders to be run
# Define the sub-folder where ground motion records are stored
# Input directory nesting: CurrentDirectory/GM/S688/HZL1/GMData 
set gmset NR
set inpath GM/$gmset
set inpath_length [string length $inpath]
# Find and sort the interested hazard levels
# indirlist = {GM/S688/HZL1 GM/S688/HZL2 GM/S688/HZL3}
set indirlist [lsort [glob -directory $inpath -type d *]]

# Define the output path
# Output directory nesting: MSA_CS4630/S688/HZL1/RSN447841_Seismogram_s688_68_7/IDR and PFA results
set outpath $outmain/$gmset

# Create the list of ground motions to be run by looping over all the ground motion folders
set serial 0
foreach indir $indirlist {

    # Parse the input directory from the input path
    # ex. indir = GM/S688/HZL1
    #     indir_length = 13
    #     inpath_length = 7
    #     gmdir = HZL1
    set indir_length [string length $indir]
    set gmdir [string range $indir [expr {$inpath_length + 1}] [expr {$indir_length - 1}]]

    # Import information about each ground motion in the GMInfo.txt file and add the information to the
    # "gminfo_dict" dictionary
    # One example line in GMInfo: 1 RSN447841_Seismogram_s688_68_7.txt  0.025000    5764    1.051648
    set gminfofile [open $indir/GMInfo.txt r]
    while {[gets $gminfofile line] >= 0} {
        
        # Read the filename and dt
        set filename [lindex $line 1]
        set dt [lindex $line 2]
	    set scalor [lindex $line 4]

        # Count the number of points
        set numpts 0
        set gmfile [open $indir/$filename r]
        while {[gets $gmfile line] >= 0} {
            incr numpts
        }
        close $gmfile

        # Add the ground motion information to "gminfo_dict"
        dict set gminfo_dict $serial indir $indir
        dict set gminfo_dict $serial gmdir $gmdir
        dict set gminfo_dict $serial filename $filename
        dict set gminfo_dict $serial dt $dt
        dict set gminfo_dict $serial numpts $numpts
        dict set gminfo_dict $serial scalor $scalor

        incr serial
    }
    close $gminfofile
}

# Loop over all ground motion files and run the analyses for this processor
# For 3 stripes each with 100 GMs, serial = 0 ~ 299
dict for {serial gminfo} $gminfo_dict {
    if {[expr {$serial % $numprocs}] == $this_procid} {
        dict with gminfo {
            # gminfo = {$indir $gmdir $filename $dt $numpts $scalor}
	
            # Define and create the output directory
            set filename_length [string length $filename]
            # "-5" to get rid of the ".txt"
            set outfilename [string range $filename 0 [expr {$filename_length - 5}]]
            
            # outdir = MSA_CS4630/S688/HZL1/RSN447841_Seismogram_s688_68_7
            set outdir $outpath/$gmdir/$outfilename
            file mkdir $outdir

            # Run the analysis
            # stripe.txt stores $col_drift or $max_drift
            set stripe_file [open $outdir/stripe.txt w]
            # set converge_file [open $outdir/convergence.txt w]
            wipe;

            model BasicBuilder -ndm 2 -ndf 3;

            source $userdirectory/MRF31.tcl;
    
            # source EigenvalueAnalysis.tcl;
    
            source $userdirectory/DynamicAnalysis.tcl;
            # source $userdirectory/recorders_analysis_stripe.tcl;

            #source $userdirectory/create_model.tcl
	        #source $userdirectory/run_eigen.tcl
            #source $userdirectory/recorders_analysis_stripe.tcl
            close $stripe_file
        }

        # Display the status of the analysis
        set time [clock seconds]
        # puts "Processor [expr {$this_procid + 1}]: Ground Motion #[expr {$serial + 1}] complete  -  \
        #         [clock format $time -format {%D %H:%M:%S}]"
    }
}
