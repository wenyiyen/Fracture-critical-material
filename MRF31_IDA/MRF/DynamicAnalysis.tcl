##############################################################################################################
# Reagan Chandramohan                                                                                        #
# John A. Blume Earthquake Engineering Center                                                                #
# Stanford University                                                                                        #
# Last edited: 02-Jun-2015
##############################################################################################################

# Script called from "run_ida_mp_hpc.tcl".
# Define recorders and run one ground motion at one Sa(T1) level.

##############################################################################################################

# Define the recorder directory
set filename_length [string length $filename]
set outfilename [string range $filename 0 [expr {$filename_length - 4}]]
set recorderdir $outpath/$indir/$outfilename/Scale_[format "%.2f" $sat1]
file mkdir $recorderdir

# Define drift recorders for all stories
for {set story 1} {$story <= $num_stories} {incr story} {
    recorder EnvelopeDrift -file $recorderdir/story${story}_drift_env.out -iNode [lindex $ctrl_nodes \
            [expr {$story - 1}]] -jNode [lindex $ctrl_nodes $story] -dof 1 -perpDirn 2
    # recorder Drift -file $recorderdir/story${story}_drift.out -time -iNode [lindex $ctrl_nodes \
    #         [expr {$story - 1}]] -jNode [lindex $ctrl_nodes $story] -dof 1 -perpDirn 2
}

# # Define the acceleration recorders
# for {set floor 1} {$floor <= [expr $num_stories+1]} {incr floor} {
#     # recorder Node -file $recorderdir/floor${floor}_acc.out -time\
#     #  -node [lindex $ctrl_nodes $floor-1] -dof 1 accel; # relative acc
#     recorder Node -file $recorderdir/floor${floor}_acc.out -timeSeries [expr {10 + $serial}]\
#      -time -node [lindex $ctrl_nodes $floor-1] -dof 1 accel; # absolute acc
# }


# Define the ground motion time series
timeSeries Path [expr {10 + $i}] -dt $dt -filePath $inpath/$indir/$filename -factor [expr {$g*$sat1/$sat1_gm}]
set eq_load_pattern 3
pattern UniformExcitation $eq_load_pattern 1 -accel [expr {10 + $i}]

# ------ changed below --------

# initialize collapse / convergence variables
set collapse_flag false;
set tLast 0.0;
set max_drift_last 0.0;
set state AtStart;
# Build recorders for beams / columns
source buildRecorderBeam.tcl; # change here
source buildRecorderColumn.tcl;

# save log file
logFile $recorderdir/log.out;

# run analysis (implicit)
source SolverGeneral.tcl;

# append result to tolerance note file
set tol_note_file [open $outpath/$indir/$outfilename/tolerance_note.txt a];
puts $tol_note_file "[format "%.3f" $sat1]\t[format "%.5f" $max_drift]\t$tier\t[format "%.5f" $tLast]\t[format "%.5f" $max_drift_last]\t[format "%s" $state]\t[format "%.0e" $elemConvTol]"; # changed here
close $tol_note_file