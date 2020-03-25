## 
# SolverGeneral.tcl
# Kuanshi Zhong, kuanshi@stanford.edu
# 09/2019
##

# PDE System 
constraints Transformation;
numberer RCM;
# LinearSOE type for sparse matrix systems (using SuperLU)
# Other options can be "BandGeneral", "BandSPD", "FullGeneral" (no cut-down in memory)
system SparseGEN;
# Start with classical Newton-Raphson for nonlinear solution
algorithm Newton;
# No numerical (pseudo) damping
integrator Newmark 0.5 0.25;
# Analysis type
analysis Transient;
# Test type
set testPack {NormDispIncr EnergyIncr NormUnbalance RelativeNormUnbalance RelativeEnergyIncr};
set testtag 0;
# Start with Norm incremental displacement
set testType [lindex $testPack $testtag];
# Initialize test tolerance
set testTol 1.0e-6;
# Define the maximum tolerance
set testTolMax 1.0e-3;
# Iteration steps
set testIter 500;
# Construct the default convergence test
test $testType $testTol $testIter;

# Time step
# Initial trial step (0.005 is used for higher modes if the recording dt is large)
set dtAna [expr min(0.005,$dt)];
# Minimum time step (before report Newmark convergence issue)
set dtMin 1.0e-8;
# Maximum time step to recover (same as the initial)
set dtMax $dtAna;
# Reduction factor of time step
set dtReF 1.5; # 2.0 is common to cut a half, but can be smaller

# Initiate analysis states
set ok 0;
set tFinal [expr $numpts*$dt];
set tCurrent [getTime];

# Switch recorders on
record

# Start time-history analysis
while {$ok == 0 && $tCurrent < $tFinal  && !$collapse_flag} {
    # March one time step
    set ok [analyze 1 $dtAna];
    # Loop-1: loose tolerance
    while {$ok != 0 && $testTol <= $testTolMax} {
	    # Loop-2: different test methods
	    while {$ok != 0 && $testtag < [llength $testPack]} {
		    set testType [lindex $testPack $testtag];
			puts "Use $testType";
		    test $testType $testTol $testIter;
            # Loop-3: decrease time step
            while {$ok != 0 && [expr $dtAna/$dtReF] >= $dtMin} {
                # Decrease time step
                set dtAna [expr $dtAna/$dtReF];
                puts [format "Reducing time step size (dtNew = %1.6e)" $dtAna];
                # Make a trial
                set ok [analyze 1 $dtAna];
                # At this step, try other improved Newton methods
                # Try NewtonLineSearch
                if {$ok != 0} {
                    puts "Try NewtonLineSearch";
                    algorithm NewtonLineSearch 0.6;
                    set ok [analyze 1 $dtAna];
                    # switch back Newton
                    algorithm Newton
                }
                # Try ModifiedNewton
                if {$ok != 0} {
                    puts "Try ModifiedNewton";
                    algorithm ModifiedNewton -initial;
                    set ok [analyze 1 $dtAna];
                    # switch back Newton
                    algorithm Newton
                }
                # Try Broyden
                if {$ok != 0} {
                    puts "Try Broyden";
                    algorithm Broyden 10
                    set ok [analyze 1 $dtAna];
                    # switch back Newton
                    algorithm Newton
                }
                # Try KrylovNewton
                if {$ok != 0} {
                    puts "Try KrylovNewton";
                    algorithm KrylovNewton;
                    set ok [analyze 1 $dtAna];
                    # switch back Newton
                    algorithm Newton;
                }
                # Try BFGS
                if {$ok != 0} {
                    puts "Try BFGS";
                    algorithm BFGS;
                    set ok [analyze 1 $dtAna];
                    # switch back Newton
                    algorithm Newton;
                }
            }
			if {$ok != 0} {
			    # If these algorithm fails, try to use a different test method
                set testtag [expr $testtag+1];
				# resume dt
				set dtAna $dtMax;
				puts "Resume dt = $dtAna";
			} else {
			    # recover the initial test method
			    set testtag 0;
                set testType [lindex $testPack $testtag];
                test $testType $testTol $testIter;
				puts "Use $testType";
				# resume dt
				set dtAna $dtMax;
				puts "Resume dt = $dtAna";
			}
		}
        if {$ok != 0} {
			# Increase tolerance and resume max time step
			set testTol [expr $testTol*10.0];
			puts "Increase tolerance to $testTol";
			set dtAna $dtMax;
			puts "Resume dt = $dtAna";
			test $testType $testTol $testIter;
			set ok [analyze 1 $dtAna];
		}
	}
    if {$ok == 0 && [expr $dtAna*$dtReF] <= $dtMax} {
	    # Recover the dtAna if current step converged but no more than dtMax
	    set dtAna [expr $dtAna*$dtReF];
	    # resume the initial testTol
		set testtag 0;
		set testType [lindex $testPack $testtag];
	    set testTol 1.0e-6;
		puts "Resume tolerance $testTol";
	    test $testType $testTol $testIter;
		puts "Use $testType";
        puts [format "Increasing time step size (dtNew = %1.6e)" $dtAna]
	}  
	set tCurrent [getTime];
	# "in-analysis collapse check" before the GM input ends, changed here
	if {$ok == 0} {
        set max_drift [max_drift_model $ctrl_nodes];
        if {$max_drift >= $col_drift} {
            set collapse_flag true;
            set tLast $tCurrent;
            set max_drift_last $max_drift;
            set state "exceedCollapseDrift";
        }
    }
}

if {$ok != 0} {
    puts "Wipe analysis"
    wipeAnalysis;
	puts "Try GeneralizedAlpha integrator";
	constraints Transformation;
	numberer RCM;
	# LinearSOE type for sparse matrix systems (using SuperLU)
	# Other options can be "BandGeneral", "BandSPD", "FullGeneral" (no cut-down in memory)
	system SparseGEN;
    # Introduce numerical (pseudo) damping
    integrator GeneralizedAlpha 0.8 0.5;
	analysis Transient;
	# Test type
	set testtag 0;
	# Start with Norm incremental displacement
	set testType [lindex $testPack $testtag];
	# Initialize test tolerance
	set testTol 1.0e-6;
	# Define the maximum tolerance
	set testTolMax 1.0e-3;
	# Construct the default convergence test
	test $testType $testTol $testIter;

	# Time step
	# Initial trial step (0.005 is used for higher modes if the recording dt is large)
	set dtAna [expr min(0.005,$dt)];
	# Reduction factor of time step
	set dtReF 1.5; # 2.0 is common to cut a half, but can be smaller

	# Initiate analysis states
	set ok 0;
	set tFinal [expr $numpts*$dt];
	set tCurrent [getTime];

	# Switch recorders on
	record

	# Start time-history analysis
	while {$ok == 0 && $tCurrent < $tFinal && !$collapse_flag} { # changed here
		# March one time step
		set ok [analyze 1 $dtAna];
		
		# Loop-1: loose tolerance
		while {$ok != 0 && $testTol <= $testTolMax} {
			# Loop-2: different test methods
			while {$ok != 0 && $testtag < [llength $testPack]} {
				set testType [lindex $testPack $testtag];
				puts "Use $testType";
				test $testType $testTol $testIter;
				# Loop-3: decrease time step
				while {$ok != 0 && [expr $dtAna/$dtReF] >= $dtMin} {
					# Decrease time step
					set dtAna [expr $dtAna/$dtReF];
					puts [format "Reducing time step size (dtNew = %1.6e)" $dtAna];
					# Make a trial
					set ok [analyze 1 $dtAna];
					# At this step, try other improved Newton methods
					# Try NewtonLineSearch
					if {$ok != 0} {
						puts "Try NewtonLineSearch";
						algorithm NewtonLineSearch 0.6;
						set ok [analyze 1 $dtAna];
						# switch back Newton
						algorithm Newton
					}
					# Try ModifiedNewton
					if {$ok != 0} {
						puts "Try ModifiedNewton";
						algorithm ModifiedNewton -initial;
						set ok [analyze 1 $dtAna];
						# switch back Newton
						algorithm Newton
					}
					# Try Broyden
					if {$ok != 0} {
						puts "Try Broyden";
						algorithm Broyden 10
						set ok [analyze 1 $dtAna];
						# switch back Newton
						algorithm Newton
					}
					# Try KrylovNewton
					if {$ok != 0} {
						puts "Try KrylovNewton";
						algorithm KrylovNewton;
						set ok [analyze 1 $dtAna];
						# switch back Newton
						algorithm Newton;
					}
					# Try BFGS
					if {$ok != 0} {
						puts "Try BFGS";
						algorithm BFGS;
						set ok [analyze 1 $dtAna];
						# switch back Newton
						algorithm Newton;
					}
				}
				if {$ok != 0} {
					# If these algorithm fails, try to use a different test method
					set testtag [expr $testtag+1];
					# resume dt
					set dtAna $dtMax;
					puts "Resume dt = $dtAna";
				} else {
					# recover the initial test method
					set testtag 0;
					set testType [lindex $testPack $testtag];
					test $testType $testTol $testIter;
					puts "Use $testType";
					# resume dt
					set dtAna $dtMax;
					puts "Resume dt = $dtAna";
				}
			}
			if {$ok != 0} {
				# Increase tolerance and resume max time step
				set testTol [expr $testTol*10.0];
				puts "Increase tolerance to $testTol";
				set dtAna $dtMax;
				puts "Resume dt = $dtAna";
				test $testType $testTol $testIter;
				set ok [analyze 1 $dtAna];
			}
		}
		if {$ok == 0 && [expr $dtAna*$dtReF] <= $dtMax} {
			# Recover the dtAna if current step converged but no more than dtMax
			set dtAna [expr $dtAna*$dtReF];
			# resume the initial testTol
			set testtag 0;
			set testType [lindex $testPack $testtag];
			set testTol 1.0e-6;
			puts "Resume tolerance $testTol";
			test $testType $testTol $testIter;
			puts "Use $testType";
			puts [format "Increasing time step size (dtNew = %1.6e)" $dtAna]
		}  
		set tCurrent [getTime];
		# "in-analysis collapse check" before the GM input ends, changed here
		if {$ok == 0} {
	        set max_drift [max_drift_model $ctrl_nodes];
	        if {$max_drift >= $col_drift} {
	            set collapse_flag true;
	            set tLast $tCurrent;
	            set max_drift_last $max_drift;
	            set state "exceedCollapseDrift";
	        }
	    }
	}
}

if {$ok != 0} {
	# changed here
	set max_drift inf; # max_drift == inf, then there is convergence issue
	set collapse_flag true;
	set tLast $tCurrent;
	set max_drift_last [max_drift_outfile $recorderdir $num_stories];
	set state "inconvergence";
    puts "Convergence issue at $tCurrent";
} elseif {$ok == 0 && !$collapse_flag} {
	# changed here
	set max_drift [max_drift_outfile $recorderdir $num_stories];
	set tLast $tCurrent;
	set max_drift_last $max_drift;
	set state "completeNotCollapse";
    puts "Complete analysis";
}