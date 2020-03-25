############################################################################
# This script defines the fracture materials for flanges.
# created by: Wen-Yi Yen
# last edit: 18/05/25
############################################################################

proc flangeMat {matTag elBilinTag minMaxTag elPPGapTag maxStrain E Fy eta} {

# tension side
set EP1 $E;
set EP2 0;
set epsP2 [expr $Fy/$E];
set EN1 0;
set EN2 0;
set epsN2 [expr $Fy/$E];
uniaxialMaterial ElasticBilin $elBilinTag $EP1 $EP2 $epsP2 $EN1 $EN2 $epsN2;
uniaxialMaterial MinMax $minMaxTag $elBilinTag -max $maxStrain;

# compression side
set gap 0.0;
uniaxialMaterial ElasticPPGap $elPPGapTag $E $Fy $gap $eta damage;

# parallel
uniaxialMaterial Parallel $matTag $minMaxTag $elPPGapTag;

}

proc flangeMat2 {matTag matTag2 elTag minMaxTag zeroTensionTag steel01Tag maxStress E Fy eta} {

# elastic and minmax
uniaxialMaterial Elastic $elTag [expr $E*100.0];
uniaxialMaterial MinMax $minMaxTag $elTag -min -3.00000E+7 -max [expr $maxStress/$E/100.0];

# zero tension
uniaxialMaterial ENT $zeroTensionTag [expr $E*100];

# steel01
uniaxialMaterial Steel01 $steel01Tag $Fy $E $eta;

# series and parallel
uniaxialMaterial Parallel $matTag2 $minMaxTag $zeroTensionTag;
# uniaxialMaterial Parallel $matTag $minMaxTag $zeroTensionTag;
uniaxialMaterial Series $matTag $steel01Tag $matTag2;

}

############################################################################
# This script builds the section for beam hinges with "fracture" properties
# created by: Wen-Yi Yen
# last edit: 18/05/25
############################################################################

proc fracHingeSection {secTag topFMatTag botFMatTag webMatTag webMatTag1 d bf tf tw tabThickness boltLocation boltDiameter FuvBolt FuTab FyTab Lc} {

set widthDivNum 5; # need to be a variable
set effectiveWidth 3.0; # need to be a variable

# ----------- web material -------------------
set E 29000.0;
set as 0.03; # change here?
set FyP $FyTab;
set FyN -$FyTab;
set Lamda 10000;
set c 1.0;
set th_pP 20;
set th_pN 20;
set th_pcP 20;
set th_pcN 20;
set ResP 1.0;
set ResN 1.0;
set th_uP 100;
set th_uN 100;
set FprP 0.3;
set FprN 0.3;
set Apinch 0.8;

uniaxialMaterial ModIMKPinching $webMatTag1 $E $as $as $FyP $FyN $FprP $FprN $Apinch $Lamda $Lamda $Lamda $Lamda\
 $c $c $c $c $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN 1.0 1.0;

set AWebFiber [expr ($effectiveWidth-$boltDiameter)*$tabThickness];
set pi [expr acos(-1)];
set boltStressLimit [expr $FuvBolt*$pi*$boltDiameter*$boltDiameter/4/$AWebFiber]; # changed here
set tabStressLimit [expr ($effectiveWidth-$boltDiameter)*$tabThickness*$FuTab/$AWebFiber];
set bearingStressLimit [expr min(1.2*$Lc*$tabThickness*$FuTab, 2.4*$boltDiameter*$tabThickness*$FuTab)];
set webStressLimit [expr min($boltStressLimit, $tabStressLimit, $bearingStressLimit)];
set webStrainLimit [expr $FyP/$E + max(0, ($webStressLimit-$FyP)/($E*$as))];
# set webStrainLimit [expr $FyP/$E + max(0, (110.0-$FyP)/($E*$as))];
set webStrainLimit 10000.0; # debug, change here
# set webStrainLimit 0.02; # debug, change here

# puts $boltStressLimit;
# puts $tabStressLimit;
# puts $bearingStressLimit;
# puts $webStressLimit;
# puts $webStrainLimit;
# puts [expr $FyP/$E];
# puts [expr max(0, ($webStressLimit-$FyP)/($E*$as))];

# uniaxialMaterial MinMax $webMatTag $webMatTag1 -min [expr -$webStrainLimit] -max $webStrainLimit;
# set webMatTag $webMatTag1;

# uniaxialMaterial Steel01 $webMatTag $FyTab 29000.0 0.01;

# changed here
set matTag 100;
set E 2119.0;
set as 0.35; # change here?
set FyP 8.9;
set FyN -8.9;
set Lamda 10000.0;
set c 1.0;
set th_pP 0.045;
set th_pN 0.045;
set th_pcP 0.01;
set th_pcN 0.01;
set ResP 0.5;
set ResN 0.5;
set th_uP 100;
set th_uN 100;
set FprP 0.2;
set FprN 0.2;
set Apinch 0.8;
uniaxialMaterial ModIMKPinching $webMatTag $E $as $as $FyP $FyN $FprP $FprN $Apinch $Lamda $Lamda $Lamda $Lamda\
 $c $c $c $c $th_pP $th_pN $th_pcP $th_pcN $ResP $ResN $th_uP $th_uN 1.0 1.0;

# -------------------------------------------------------------

# fiber area and location (flange)
set ybot [expr ($d-$tf)/2]; # think about the local axis - bottom flange in tension at start
set ytop [expr -$ybot];
set z 0.0;
set AFlange [expr $bf*$tf];

# define fiber section
section Fiber $secTag {
	
	# top and bottom flanges # changed here
	fiber $ytop $z $AFlange $topFMatTag;
	fiber $ybot $z $AFlange $botFMatTag;

	# patch rect $topFMatTag 4 1 $yI $zI $yJ $zJ;
	# patch rect $botFMatTag 4 1 $yI $zI $yJ $zJ;

	# fiber [expr $ytop+3.0/8*$tf] $z [expr $AFlange/4] $topFMatTag;
	# fiber [expr $ytop+1.0/8*$tf] $z [expr $AFlange/4] $topFMatTag;
	# fiber [expr $ytop-3.0/8*$tf] $z [expr $AFlange/4] $topFMatTag;
	# fiber [expr $ytop-1.0/8*$tf] $z [expr $AFlange/4] $topFMatTag;

	# fiber [expr $ybot+3.0/8*$tf] $z [expr $AFlange/4] $botFMatTag;
	# fiber [expr $ybot+1.0/8*$tf] $z [expr $AFlange/4] $botFMatTag;
	# fiber [expr $ybot-3.0/8*$tf] $z [expr $AFlange/4] $botFMatTag;
	# fiber [expr $ybot-1.0/8*$tf] $z [expr $AFlange/4] $botFMatTag;

	# web springs
	for {set i 0} {$i < [llength $boltLocation]} {incr i} {
		set yCenter [lindex $boltLocation $i];
		set numSubdivY $widthDivNum;
		set numSubdivZ 1;
		set yI [expr $yCenter-$effectiveWidth/2];
		set zI [expr -$tabThickness/2];
		set yJ [expr $yCenter+$effectiveWidth/2];
		set zJ [expr $tabThickness/2];
		patch rect $webMatTag $numSubdivY $numSubdivZ $yI $zI $yJ $zJ;
	}

	# # web springs
	# for {set i 0} {$i < [llength $boltLocation]} {incr i} {
	# 	set y [lindex $boltLocation $i];
	# 	fiber $y $z $AWebFiber $webMatTag;
	# }

}


}

############################################################################
# This script builds the section for beam hinges with "fracture" properties in both top and bottom flange
# created by: Wen-Yi Yen
# last edit: 18/06/30
############################################################################

proc preNRBeam1 { eleTag node1 node2 hingeLgth transfTag elSecList flangeMatLftTopList flangeMatLftBotList flangeMatRgtTopList flangeMatRgtBotList fracHingeSecLftList fracHingeSecRgtList} {

# elastic section for beam
section Elastic [lindex $elSecList 0] [lindex $elSecList 1] [lindex $elSecList 2] [lindex $elSecList 3];

# left hinge section for beam
flangeMat [lindex $flangeMatLftTopList 0] [lindex $flangeMatLftTopList 1] [lindex $flangeMatLftTopList 2] [lindex $flangeMatLftTopList 3]\
		[lindex $flangeMatLftTopList 4] [lindex $flangeMatLftTopList 5] [lindex $flangeMatLftTopList 6] [lindex $flangeMatLftTopList 7];

flangeMat [lindex $flangeMatLftBotList 0] [lindex $flangeMatLftBotList 1] [lindex $flangeMatLftBotList 2] [lindex $flangeMatLftBotList 3]\
		[lindex $flangeMatLftBotList 4] [lindex $flangeMatLftBotList 5] [lindex $flangeMatLftBotList 6] [lindex $flangeMatLftBotList 7];

fracHingeSection [lindex $fracHingeSecLftList 0] [lindex $fracHingeSecLftList 1] [lindex $fracHingeSecLftList 2] [lindex $fracHingeSecLftList 3]\
		[lindex $fracHingeSecLftList 4] [lindex $fracHingeSecLftList 5] [lindex $fracHingeSecLftList 6] [lindex $fracHingeSecLftList 7]\
		[lindex $fracHingeSecLftList 8] [lindex $fracHingeSecLftList 9] [lindex $fracHingeSecLftList 10] [lindex $fracHingeSecLftList 11]\
		[lindex $fracHingeSecLftList 12] [lindex $fracHingeSecLftList 13] [lindex $fracHingeSecLftList 14] [lindex $fracHingeSecLftList 15];

# right hinge section for beam
flangeMat [lindex $flangeMatRgtTopList 0] [lindex $flangeMatRgtTopList 1] [lindex $flangeMatRgtTopList 2] [lindex $flangeMatRgtTopList 3]\
		[lindex $flangeMatRgtTopList 4] [lindex $flangeMatRgtTopList 5] [lindex $flangeMatRgtTopList 6] [lindex $flangeMatRgtTopList 7];

flangeMat [lindex $flangeMatRgtBotList 0] [lindex $flangeMatRgtBotList 1] [lindex $flangeMatRgtBotList 2] [lindex $flangeMatRgtBotList 3]\
		[lindex $flangeMatRgtBotList 4] [lindex $flangeMatRgtBotList 5] [lindex $flangeMatRgtBotList 6] [lindex $flangeMatRgtBotList 7];

fracHingeSection [lindex $fracHingeSecRgtList 0] [lindex $fracHingeSecRgtList 1] [lindex $fracHingeSecRgtList 2] [lindex $fracHingeSecRgtList 3]\
		[lindex $fracHingeSecRgtList 4] [lindex $fracHingeSecRgtList 5] [lindex $fracHingeSecRgtList 6] [lindex $fracHingeSecRgtList 7]\
		[lindex $fracHingeSecRgtList 8] [lindex $fracHingeSecRgtList 9] [lindex $fracHingeSecRgtList 10] [lindex $fracHingeSecRgtList 11]\
		[lindex $fracHingeSecRgtList 12] [lindex $fracHingeSecRgtList 13] [lindex $fracHingeSecRgtList 14] [lindex $fracHingeSecRgtList 15];

# create element
set elSecID [lindex $elSecList 0];
set hingeSecLftID [lindex $fracHingeSecLftList 0];
set hingeSecRgtID [lindex $fracHingeSecRgtList 0];

# set integration "HingeMidpoint $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
set integration "HingeRadau $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
element forceBeamColumn $eleTag $node1 $node2 $transfTag $integration;

}

############################################################################
# This script builds the section for beam hinges with "fracture" properties in bottom flanges only.
# The top flange material is defined as bilinear.
# created by: Wen-Yi Yen
# last edit: 18/06/30
############################################################################

proc preNRBeam2 { eleTag node1 node2 hingeLgth transfTag elSecList flangeMatLftTopList flangeMatLftBotList flangeMatRgtTopList flangeMatRgtBotList fracHingeSecLftList fracHingeSecRgtList} {

# elastic section for beam
# set elSecList [list $elSecTag $E $Abeam $Ibeam];
section Elastic [lindex $elSecList 0] [lindex $elSecList 1] [lindex $elSecList 2] [lindex $elSecList 3];

# left hinge section for beam
# set flangeMatLftTopList [list $topFMatTag $Fy $E $eta]; # recommend using eta = 0.05
uniaxialMaterial Steel01 [lindex $flangeMatLftTopList 0] [lindex $flangeMatLftTopList 1] [lindex $flangeMatLftTopList 2] [lindex $flangeMatLftTopList 3];

# set flangeMatLftBotList [list $botFMatTag $elBilinTag $minMaxTag $elPPGapTag $maxStrain $E $Fy $eta]; # recommend using eta = 0.02
flangeMat [lindex $flangeMatLftBotList 0] [lindex $flangeMatLftBotList 1] [lindex $flangeMatLftBotList 2] [lindex $flangeMatLftBotList 3]\
		[lindex $flangeMatLftBotList 4] [lindex $flangeMatLftBotList 5] [lindex $flangeMatLftBotList 6] [lindex $flangeMatLftBotList 7];

# set fracHingeSecLftList [list $secTag $topFMatTag $botFMatTag $webMatTag $webMatTag1 $d $bf $tf $tw $tabThickness $boltLocation $boltDiameter $FuBolt $FuTab $FyTab $Lc];
fracHingeSection [lindex $fracHingeSecLftList 0] [lindex $fracHingeSecLftList 1] [lindex $fracHingeSecLftList 2] [lindex $fracHingeSecLftList 3]\
		[lindex $fracHingeSecLftList 4] [lindex $fracHingeSecLftList 5] [lindex $fracHingeSecLftList 6] [lindex $fracHingeSecLftList 7]\
		[lindex $fracHingeSecLftList 8] [lindex $fracHingeSecLftList 9] [lindex $fracHingeSecLftList 10] [lindex $fracHingeSecLftList 11]\
		[lindex $fracHingeSecLftList 12] [lindex $fracHingeSecLftList 13] [lindex $fracHingeSecLftList 14] [lindex $fracHingeSecLftList 15];

# right hinge section for beam
uniaxialMaterial Steel01 [lindex $flangeMatRgtTopList 0] [lindex $flangeMatRgtTopList 1] [lindex $flangeMatRgtTopList 2] [lindex $flangeMatRgtTopList 3];

flangeMat [lindex $flangeMatRgtBotList 0] [lindex $flangeMatRgtBotList 1] [lindex $flangeMatRgtBotList 2] [lindex $flangeMatRgtBotList 3]\
		[lindex $flangeMatRgtBotList 4] [lindex $flangeMatRgtBotList 5] [lindex $flangeMatRgtBotList 6] [lindex $flangeMatRgtBotList 7];

fracHingeSection [lindex $fracHingeSecRgtList 0] [lindex $fracHingeSecRgtList 1] [lindex $fracHingeSecRgtList 2] [lindex $fracHingeSecRgtList 3]\
		[lindex $fracHingeSecRgtList 4] [lindex $fracHingeSecRgtList 5] [lindex $fracHingeSecRgtList 6] [lindex $fracHingeSecRgtList 7]\
		[lindex $fracHingeSecRgtList 8] [lindex $fracHingeSecRgtList 9] [lindex $fracHingeSecRgtList 10] [lindex $fracHingeSecRgtList 11]\
		[lindex $fracHingeSecRgtList 12] [lindex $fracHingeSecRgtList 13] [lindex $fracHingeSecRgtList 14] [lindex $fracHingeSecRgtList 15];

# create element
set elSecID [lindex $elSecList 0];
set hingeSecLftID [lindex $fracHingeSecLftList 0];
set hingeSecRgtID [lindex $fracHingeSecRgtList 0];

set integration "HingeMidpoint $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
# set integration "HingeRadau $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
element forceBeamColumn $eleTag $node1 $node2 $transfTag $integration;

}

#######################
proc preNRBeam3 { eleTag node1 node2 hingeLgth transfTag elSecList flangeMatLftTopList flangeMatLftBotList flangeMatRgtTopList flangeMatRgtBotList fracHingeSecLftList fracHingeSecRgtList} {

# elastic section for beam
# set elSecList [list $elSecTag $E $Abeam $Ibeam];
section Elastic [lindex $elSecList 0] [lindex $elSecList 1] [lindex $elSecList 2] [lindex $elSecList 3];

# left hinge section for beam
# set flangeMatLftTopList [list $topFMatTag $Fy $E $eta]; # recommend using eta = 0.05
uniaxialMaterial Steel01 [lindex $flangeMatLftTopList 0] [lindex $flangeMatLftTopList 1] [lindex $flangeMatLftTopList 2] [lindex $flangeMatLftTopList 3];

# # set flangeMatLftBotList [list $botFMatTag $elBilinTag $minMaxTag $elPPGapTag $maxStrain $E $Fy $eta]; # recommend using eta = 0.02
# flangeMat [lindex $flangeMatLftBotList 0] [lindex $flangeMatLftBotList 1] [lindex $flangeMatLftBotList 2] [lindex $flangeMatLftBotList 3]\
# 		[lindex $flangeMatLftBotList 4] [lindex $flangeMatLftBotList 5] [lindex $flangeMatLftBotList 6] [lindex $flangeMatLftBotList 7];

uniaxialMaterial Steel01 [lindex $flangeMatLftBotList 0] [lindex $flangeMatLftBotList 1] [lindex $flangeMatLftBotList 2] [lindex $flangeMatLftBotList 3];

# set fracHingeSecLftList [list $secTag $topFMatTag $botFMatTag $webMatTag $webMatTag1 $d $bf $tf $tw $tabThickness $boltLocation $boltDiameter $FuBolt $FuTab $FyTab $Lc];
fracHingeSection [lindex $fracHingeSecLftList 0] [lindex $fracHingeSecLftList 1] [lindex $fracHingeSecLftList 2] [lindex $fracHingeSecLftList 3]\
		[lindex $fracHingeSecLftList 4] [lindex $fracHingeSecLftList 5] [lindex $fracHingeSecLftList 6] [lindex $fracHingeSecLftList 7]\
		[lindex $fracHingeSecLftList 8] [lindex $fracHingeSecLftList 9] [lindex $fracHingeSecLftList 10] [lindex $fracHingeSecLftList 11]\
		[lindex $fracHingeSecLftList 12] [lindex $fracHingeSecLftList 13] [lindex $fracHingeSecLftList 14] [lindex $fracHingeSecLftList 15];

# right hinge section for beam
uniaxialMaterial Steel01 [lindex $flangeMatRgtTopList 0] [lindex $flangeMatRgtTopList 1] [lindex $flangeMatRgtTopList 2] [lindex $flangeMatRgtTopList 3];

# flangeMat [lindex $flangeMatRgtBotList 0] [lindex $flangeMatRgtBotList 1] [lindex $flangeMatRgtBotList 2] [lindex $flangeMatRgtBotList 3]\
# 		[lindex $flangeMatRgtBotList 4] [lindex $flangeMatRgtBotList 5] [lindex $flangeMatRgtBotList 6] [lindex $flangeMatRgtBotList 7];

uniaxialMaterial Steel01 [lindex $flangeMatRgtBotList 0] [lindex $flangeMatRgtBotList 1] [lindex $flangeMatRgtBotList 2] [lindex $flangeMatRgtBotList 3];

fracHingeSection [lindex $fracHingeSecRgtList 0] [lindex $fracHingeSecRgtList 1] [lindex $fracHingeSecRgtList 2] [lindex $fracHingeSecRgtList 3]\
		[lindex $fracHingeSecRgtList 4] [lindex $fracHingeSecRgtList 5] [lindex $fracHingeSecRgtList 6] [lindex $fracHingeSecRgtList 7]\
		[lindex $fracHingeSecRgtList 8] [lindex $fracHingeSecRgtList 9] [lindex $fracHingeSecRgtList 10] [lindex $fracHingeSecRgtList 11]\
		[lindex $fracHingeSecRgtList 12] [lindex $fracHingeSecRgtList 13] [lindex $fracHingeSecRgtList 14] [lindex $fracHingeSecRgtList 15];

# create element
set elSecID [lindex $elSecList 0];
set hingeSecLftID [lindex $fracHingeSecLftList 0];
set hingeSecRgtID [lindex $fracHingeSecRgtList 0];

set integration "HingeMidpoint $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
# set integration "HingeRadau $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
element forceBeamColumn $eleTag $node1 $node2 $transfTag $integration;

}

################ with SteelFractureDI ###########################
proc preNRBeam4 { eleTag node1 node2 hingeLgth transfTag elSecList flangeMatLftTopList flangeMatLftBotList flangeMatRgtTopList flangeMatRgtBotList fracHingeSecLftList fracHingeSecRgtList} {

# elastic section for beam
# set elSecList [list $elSecTag $E $Abeam $Ibeam];
section Elastic [lindex $elSecList 0] [lindex $elSecList 1] [lindex $elSecList 2] [lindex $elSecList 3];

# left hinge section for beam
# set flangeMatLftTopList [list $topFMatTag $Fy $E $eta]; # recommend using eta = 0.05
# SteelFractureDI $matTag $Fy $E0 $b $R0 $cR1 $cR2 $a1 $a2 $a3 $a4 $sigcr $m
uniaxialMaterial SteelFractureDI [lindex $flangeMatLftTopList 0] [lindex $flangeMatLftTopList 1] [lindex $flangeMatLftTopList 2] [lindex $flangeMatLftTopList 3] \
		[lindex $flangeMatLftTopList 4] [lindex $flangeMatLftTopList 5] [lindex $flangeMatLftTopList 6] [lindex $flangeMatLftTopList 7] \
		[lindex $flangeMatLftTopList 8] [lindex $flangeMatLftTopList 9] [lindex $flangeMatLftTopList 10] [lindex $flangeMatLftTopList 11] \
		[lindex $flangeMatLftTopList 12];

# set flangeMatLftBotList [list $botFMatTag $elBilinTag $minMaxTag $elPPGapTag $maxStrain $E $Fy $eta]; # recommend using eta = 0.02
uniaxialMaterial SteelFractureDI [lindex $flangeMatLftBotList 0] [lindex $flangeMatLftBotList 1] [lindex $flangeMatLftBotList 2] [lindex $flangeMatLftBotList 3] \
		[lindex $flangeMatLftBotList 4] [lindex $flangeMatLftBotList 5] [lindex $flangeMatLftBotList 6] [lindex $flangeMatLftBotList 7] \
		[lindex $flangeMatLftBotList 8] [lindex $flangeMatLftBotList 9] [lindex $flangeMatLftBotList 10] [lindex $flangeMatLftBotList 11] \
		[lindex $flangeMatLftBotList 12];

# set fracHingeSecLftList [list $secTag $topFMatTag $botFMatTag $webMatTag $webMatTag1 $d $bf $tf $tw $tabThickness $boltLocation $boltDiameter $FuBolt $FuTab $FyTab $Lc];
fracHingeSection [lindex $fracHingeSecLftList 0] [lindex $fracHingeSecLftList 1] [lindex $fracHingeSecLftList 2] [lindex $fracHingeSecLftList 3]\
		[lindex $fracHingeSecLftList 4] [lindex $fracHingeSecLftList 5] [lindex $fracHingeSecLftList 6] [lindex $fracHingeSecLftList 7]\
		[lindex $fracHingeSecLftList 8] [lindex $fracHingeSecLftList 9] [lindex $fracHingeSecLftList 10] [lindex $fracHingeSecLftList 11]\
		[lindex $fracHingeSecLftList 12] [lindex $fracHingeSecLftList 13] [lindex $fracHingeSecLftList 14] [lindex $fracHingeSecLftList 15];

# right hinge section for beam
# uniaxialMaterial Steel01 [lindex $flangeMatRgtTopList 0] [lindex $flangeMatRgtTopList 1] [lindex $flangeMatRgtTopList 2] [lindex $flangeMatRgtTopList 3];

uniaxialMaterial SteelFractureDI [lindex $flangeMatRgtTopList 0] [lindex $flangeMatRgtTopList 1] [lindex $flangeMatRgtTopList 2] [lindex $flangeMatRgtTopList 3] \
		[lindex $flangeMatRgtTopList 4] [lindex $flangeMatRgtTopList 5] [lindex $flangeMatRgtTopList 6] [lindex $flangeMatRgtTopList 7] \
		[lindex $flangeMatRgtTopList 8] [lindex $flangeMatRgtTopList 9] [lindex $flangeMatRgtTopList 10] [lindex $flangeMatRgtTopList 11] \
		[lindex $flangeMatRgtTopList 12];

uniaxialMaterial SteelFractureDI [lindex $flangeMatRgtBotList 0] [lindex $flangeMatRgtBotList 1] [lindex $flangeMatRgtBotList 2] [lindex $flangeMatRgtBotList 3] \
		[lindex $flangeMatRgtBotList 4] [lindex $flangeMatRgtBotList 5] [lindex $flangeMatRgtBotList 6] [lindex $flangeMatRgtBotList 7] \
		[lindex $flangeMatRgtBotList 8] [lindex $flangeMatRgtBotList 9] [lindex $flangeMatRgtBotList 10] [lindex $flangeMatRgtBotList 11] \
		[lindex $flangeMatRgtBotList 12];
		
fracHingeSection [lindex $fracHingeSecRgtList 0] [lindex $fracHingeSecRgtList 1] [lindex $fracHingeSecRgtList 2] [lindex $fracHingeSecRgtList 3]\
		[lindex $fracHingeSecRgtList 4] [lindex $fracHingeSecRgtList 5] [lindex $fracHingeSecRgtList 6] [lindex $fracHingeSecRgtList 7]\
		[lindex $fracHingeSecRgtList 8] [lindex $fracHingeSecRgtList 9] [lindex $fracHingeSecRgtList 10] [lindex $fracHingeSecRgtList 11]\
		[lindex $fracHingeSecRgtList 12] [lindex $fracHingeSecRgtList 13] [lindex $fracHingeSecRgtList 14] [lindex $fracHingeSecRgtList 15];

# create element
set elSecID [lindex $elSecList 0];
set hingeSecLftID [lindex $fracHingeSecLftList 0];
set hingeSecRgtID [lindex $fracHingeSecRgtList 0];

# set integration "HingeMidpoint $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
set integration "HingeRadau $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
# element forceBeamColumn $eleTag $node1 $node2 $transfTag $integration;
element forceBeamColumn $eleTag $node1 $node2 $transfTag $integration -iter 50 1e-6;

}

################ with SteelFractureDI, displacement-based beam ###########################
proc preNRBeam5 { eleTag node1 node2 hingeLgth transfTag elSecList flangeMatLftTopList flangeMatLftBotList flangeMatRgtTopList flangeMatRgtBotList fracHingeSecLftList fracHingeSecRgtList} {

# elastic section for beam
# set elSecList [list $elSecTag $E $Abeam $Ibeam];
section Elastic [lindex $elSecList 0] [lindex $elSecList 1] [lindex $elSecList 2] [lindex $elSecList 3];

# left hinge section for beam
# set flangeMatLftTopList [list $topFMatTag $Fy $E $eta]; # recommend using eta = 0.05
# SteelFractureDI $matTag $Fy $E0 $b $R0 $cR1 $cR2 $a1 $a2 $a3 $a4 $sigcr $m
uniaxialMaterial SteelFractureDI [lindex $flangeMatLftTopList 0] [lindex $flangeMatLftTopList 1] [lindex $flangeMatLftTopList 2] [lindex $flangeMatLftTopList 3] \
		[lindex $flangeMatLftTopList 4] [lindex $flangeMatLftTopList 5] [lindex $flangeMatLftTopList 6] [lindex $flangeMatLftTopList 7] \
		[lindex $flangeMatLftTopList 8] [lindex $flangeMatLftTopList 9] [lindex $flangeMatLftTopList 10] [lindex $flangeMatLftTopList 11] \
		[lindex $flangeMatLftTopList 12];

# set flangeMatLftBotList [list $botFMatTag $elBilinTag $minMaxTag $elPPGapTag $maxStrain $E $Fy $eta]; # recommend using eta = 0.02
uniaxialMaterial SteelFractureDI [lindex $flangeMatLftBotList 0] [lindex $flangeMatLftBotList 1] [lindex $flangeMatLftBotList 2] [lindex $flangeMatLftBotList 3] \
		[lindex $flangeMatLftBotList 4] [lindex $flangeMatLftBotList 5] [lindex $flangeMatLftBotList 6] [lindex $flangeMatLftBotList 7] \
		[lindex $flangeMatLftBotList 8] [lindex $flangeMatLftBotList 9] [lindex $flangeMatLftBotList 10] [lindex $flangeMatLftBotList 11] \
		[lindex $flangeMatLftBotList 12];

# set fracHingeSecLftList [list $secTag $topFMatTag $botFMatTag $webMatTag $webMatTag1 $d $bf $tf $tw $tabThickness $boltLocation $boltDiameter $FuBolt $FuTab $FyTab $Lc];
fracHingeSection [lindex $fracHingeSecLftList 0] [lindex $fracHingeSecLftList 1] [lindex $fracHingeSecLftList 2] [lindex $fracHingeSecLftList 3]\
		[lindex $fracHingeSecLftList 4] [lindex $fracHingeSecLftList 5] [lindex $fracHingeSecLftList 6] [lindex $fracHingeSecLftList 7]\
		[lindex $fracHingeSecLftList 8] [lindex $fracHingeSecLftList 9] [lindex $fracHingeSecLftList 10] [lindex $fracHingeSecLftList 11]\
		[lindex $fracHingeSecLftList 12] [lindex $fracHingeSecLftList 13] [lindex $fracHingeSecLftList 14] [lindex $fracHingeSecLftList 15];

# right hinge section for beam
# uniaxialMaterial Steel01 [lindex $flangeMatRgtTopList 0] [lindex $flangeMatRgtTopList 1] [lindex $flangeMatRgtTopList 2] [lindex $flangeMatRgtTopList 3];

uniaxialMaterial SteelFractureDI [lindex $flangeMatRgtTopList 0] [lindex $flangeMatRgtTopList 1] [lindex $flangeMatRgtTopList 2] [lindex $flangeMatRgtTopList 3] \
		[lindex $flangeMatRgtTopList 4] [lindex $flangeMatRgtTopList 5] [lindex $flangeMatRgtTopList 6] [lindex $flangeMatRgtTopList 7] \
		[lindex $flangeMatRgtTopList 8] [lindex $flangeMatRgtTopList 9] [lindex $flangeMatRgtTopList 10] [lindex $flangeMatRgtTopList 11] \
		[lindex $flangeMatRgtTopList 12];

uniaxialMaterial SteelFractureDI [lindex $flangeMatRgtBotList 0] [lindex $flangeMatRgtBotList 1] [lindex $flangeMatRgtBotList 2] [lindex $flangeMatRgtBotList 3] \
		[lindex $flangeMatRgtBotList 4] [lindex $flangeMatRgtBotList 5] [lindex $flangeMatRgtBotList 6] [lindex $flangeMatRgtBotList 7] \
		[lindex $flangeMatRgtBotList 8] [lindex $flangeMatRgtBotList 9] [lindex $flangeMatRgtBotList 10] [lindex $flangeMatRgtBotList 11] \
		[lindex $flangeMatRgtBotList 12];
		
fracHingeSection [lindex $fracHingeSecRgtList 0] [lindex $fracHingeSecRgtList 1] [lindex $fracHingeSecRgtList 2] [lindex $fracHingeSecRgtList 3]\
		[lindex $fracHingeSecRgtList 4] [lindex $fracHingeSecRgtList 5] [lindex $fracHingeSecRgtList 6] [lindex $fracHingeSecRgtList 7]\
		[lindex $fracHingeSecRgtList 8] [lindex $fracHingeSecRgtList 9] [lindex $fracHingeSecRgtList 10] [lindex $fracHingeSecRgtList 11]\
		[lindex $fracHingeSecRgtList 12] [lindex $fracHingeSecRgtList 13] [lindex $fracHingeSecRgtList 14] [lindex $fracHingeSecRgtList 15];

# create element
set elSecID [lindex $elSecList 0];
set hingeSecLftID [lindex $fracHingeSecLftList 0];
set hingeSecRgtID [lindex $fracHingeSecRgtList 0];

# set integration "HingeMidpoint $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
# set integration "HingeRadau $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
# element dispBeamColumn $eleTag $node1 $node2 5 -sections $hingeSecLftID $elSecID $elSecID $elSecID $hingeSecRgtID $transfTag;

# build dispBeamColumns
set discritize 8.0; # change here
if {$discritize != 1} {
	# build extra nodes
	set n1Coord [nodeCoord $node1];
	set n2Coord [nodeCoord $node2];
	set xGap [expr [lindex $n2Coord 0] - [lindex $n1Coord 0]];
	set x1 [lindex $n1Coord 0];
	set y [lindex $n2Coord 1];
	for {set interN 1} {$interN <= [expr $discritize-1]} {incr interN} {
		node [expr $node1+$interN+20] [expr $x1 + $xGap / $discritize * $interN] $y;
	}
	# build dispBeamColumn elements
	for {set segment 1} {$segment <= $discritize} {incr segment} {
		set eleTagTmp [expr $eleTag+$segment-1];
		if {$segment == 1} {
			element dispBeamColumn $eleTagTmp $node1 [expr $node1+$segment+20] 5 -sections $hingeSecLftID $elSecID $elSecID $elSecID $elSecID $transfTag;
		} elseif {$segment == $discritize} {
			element dispBeamColumn $eleTagTmp [expr $node1+$segment-1+20] $node2 5 -sections $elSecID $elSecID $elSecID $elSecID $hingeSecRgtID $transfTag;
		} else {
			element dispBeamColumn $eleTagTmp [expr $node1+$segment-1+20] [expr $node1+$segment+20] 5 -sections $elSecID $elSecID $elSecID $elSecID $elSecID $transfTag;
		}	
	}
} else {
	element dispBeamColumn $eleTag $node1 $node2 5 -sections $hingeSecLftID $elSecID $elSecID $elSecID $hingeSecRgtID $transfTag;
}

}

################ Using Stillmaker's paralleling method ###########################
proc preNRBeam6 { eleTag node1 node2 hingeLgth transfTag elSecList flangeMatLftTopList flangeMatLftBotList flangeMatRgtTopList flangeMatRgtBotList fracHingeSecLftList fracHingeSecRgtList} {

# elastic section for beam
# set elSecList [list $elSecTag $E $Abeam $Ibeam];
section Elastic [lindex $elSecList 0] [lindex $elSecList 1] [lindex $elSecList 2] [lindex $elSecList 3];

# left hinge section for beam
# set flangeMatLftTopList [list $topFMatTag $Fy $E $eta]; # recommend using eta = 0.05
# SteelFractureDI $matTag $Fy $E0 $b $R0 $cR1 $cR2 $a1 $a2 $a3 $a4 $sigcr $m
flangeMat2 [lindex $flangeMatLftTopList 0] [lindex $flangeMatLftTopList 1] [lindex $flangeMatLftTopList 2] [lindex $flangeMatLftTopList 3]\
		[lindex $flangeMatLftTopList 4] [lindex $flangeMatLftTopList 5] [lindex $flangeMatLftTopList 6] [lindex $flangeMatLftTopList 7]\
		[lindex $flangeMatLftTopList 8] [lindex $flangeMatLftTopList 9];

# set flangeMatLftBotList [list $botFMatTag $elBilinTag $minMaxTag $elPPGapTag $maxStrain $E $Fy $eta]; # recommend using eta = 0.02
flangeMat2 [lindex $flangeMatLftBotList 0] [lindex $flangeMatLftBotList 1] [lindex $flangeMatLftBotList 2] [lindex $flangeMatLftBotList 3]\
		[lindex $flangeMatLftBotList 4] [lindex $flangeMatLftBotList 5] [lindex $flangeMatLftBotList 6] [lindex $flangeMatLftBotList 7]\
		[lindex $flangeMatLftBotList 8] [lindex $flangeMatLftBotList 9];

# set fracHingeSecLftList [list $secTag $topFMatTag $botFMatTag $webMatTag $webMatTag1 $d $bf $tf $tw $tabThickness $boltLocation $boltDiameter $FuBolt $FuTab $FyTab $Lc];
fracHingeSection [lindex $fracHingeSecLftList 0] [lindex $fracHingeSecLftList 1] [lindex $fracHingeSecLftList 2] [lindex $fracHingeSecLftList 3]\
		[lindex $fracHingeSecLftList 4] [lindex $fracHingeSecLftList 5] [lindex $fracHingeSecLftList 6] [lindex $fracHingeSecLftList 7]\
		[lindex $fracHingeSecLftList 8] [lindex $fracHingeSecLftList 9] [lindex $fracHingeSecLftList 10] [lindex $fracHingeSecLftList 11]\
		[lindex $fracHingeSecLftList 12] [lindex $fracHingeSecLftList 13] [lindex $fracHingeSecLftList 14] [lindex $fracHingeSecLftList 15];

# right hinge section for beam
# uniaxialMaterial Steel01 [lindex $flangeMatRgtTopList 0] [lindex $flangeMatRgtTopList 1] [lindex $flangeMatRgtTopList 2] [lindex $flangeMatRgtTopList 3];

flangeMat2 [lindex $flangeMatRgtTopList 0] [lindex $flangeMatRgtTopList 1] [lindex $flangeMatRgtTopList 2] [lindex $flangeMatRgtTopList 3]\
		[lindex $flangeMatRgtTopList 4] [lindex $flangeMatRgtTopList 5] [lindex $flangeMatRgtTopList 6] [lindex $flangeMatRgtTopList 7]\
		[lindex $flangeMatRgtTopList 8] [lindex $flangeMatRgtTopList 9];

flangeMat2 [lindex $flangeMatRgtBotList 0] [lindex $flangeMatRgtBotList 1] [lindex $flangeMatRgtBotList 2] [lindex $flangeMatRgtBotList 3]\
		[lindex $flangeMatRgtBotList 4] [lindex $flangeMatRgtBotList 5] [lindex $flangeMatRgtBotList 6] [lindex $flangeMatRgtBotList 7]\
		[lindex $flangeMatRgtBotList 8] [lindex $flangeMatRgtBotList 9];
		
fracHingeSection [lindex $fracHingeSecRgtList 0] [lindex $fracHingeSecRgtList 1] [lindex $fracHingeSecRgtList 2] [lindex $fracHingeSecRgtList 3]\
		[lindex $fracHingeSecRgtList 4] [lindex $fracHingeSecRgtList 5] [lindex $fracHingeSecRgtList 6] [lindex $fracHingeSecRgtList 7]\
		[lindex $fracHingeSecRgtList 8] [lindex $fracHingeSecRgtList 9] [lindex $fracHingeSecRgtList 10] [lindex $fracHingeSecRgtList 11]\
		[lindex $fracHingeSecRgtList 12] [lindex $fracHingeSecRgtList 13] [lindex $fracHingeSecRgtList 14] [lindex $fracHingeSecRgtList 15];

# create element
set elSecID [lindex $elSecList 0];
set hingeSecLftID [lindex $fracHingeSecLftList 0];
set hingeSecRgtID [lindex $fracHingeSecRgtList 0];

set integration "HingeMidpoint $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
# set integration "HingeRadau $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
element forceBeamColumn $eleTag $node1 $node2 $transfTag $integration;

}

# add the tolerance limit (working...)
################ with SteelFractureDI ###########################
proc preNRBeam7 { eleTag node1 node2 hingeLgth transfTag elSecList flangeMatLftTopList flangeMatLftBotList flangeMatRgtTopList flangeMatRgtBotList fracHingeSecLftList fracHingeSecRgtList elemConvTol} {

# elastic section for beam
# set elSecList [list $elSecTag $E $Abeam $Ibeam];
section Elastic [lindex $elSecList 0] [lindex $elSecList 1] [lindex $elSecList 2] [lindex $elSecList 3];

# left hinge section for beam
# set flangeMatLftTopList [list $topFMatTag $Fy $E $eta]; # recommend using eta = 0.05
# SteelFractureDI $matTag $Fy $E0 $b $R0 $cR1 $cR2 $a1 $a2 $a3 $a4 $sigcr $m
uniaxialMaterial SteelFractureDI [lindex $flangeMatLftTopList 0] [lindex $flangeMatLftTopList 1] [lindex $flangeMatLftTopList 2] [lindex $flangeMatLftTopList 3] \
		[lindex $flangeMatLftTopList 4] [lindex $flangeMatLftTopList 5] [lindex $flangeMatLftTopList 6] [lindex $flangeMatLftTopList 7] \
		[lindex $flangeMatLftTopList 8] [lindex $flangeMatLftTopList 9] [lindex $flangeMatLftTopList 10] [lindex $flangeMatLftTopList 11] \
		[lindex $flangeMatLftTopList 12];

# set flangeMatLftBotList [list $botFMatTag $elBilinTag $minMaxTag $elPPGapTag $maxStrain $E $Fy $eta]; # recommend using eta = 0.02
uniaxialMaterial SteelFractureDI [lindex $flangeMatLftBotList 0] [lindex $flangeMatLftBotList 1] [lindex $flangeMatLftBotList 2] [lindex $flangeMatLftBotList 3] \
		[lindex $flangeMatLftBotList 4] [lindex $flangeMatLftBotList 5] [lindex $flangeMatLftBotList 6] [lindex $flangeMatLftBotList 7] \
		[lindex $flangeMatLftBotList 8] [lindex $flangeMatLftBotList 9] [lindex $flangeMatLftBotList 10] [lindex $flangeMatLftBotList 11] \
		[lindex $flangeMatLftBotList 12];

# set fracHingeSecLftList [list $secTag $topFMatTag $botFMatTag $webMatTag $webMatTag1 $d $bf $tf $tw $tabThickness $boltLocation $boltDiameter $FuBolt $FuTab $FyTab $Lc];
fracHingeSection [lindex $fracHingeSecLftList 0] [lindex $fracHingeSecLftList 1] [lindex $fracHingeSecLftList 2] [lindex $fracHingeSecLftList 3]\
		[lindex $fracHingeSecLftList 4] [lindex $fracHingeSecLftList 5] [lindex $fracHingeSecLftList 6] [lindex $fracHingeSecLftList 7]\
		[lindex $fracHingeSecLftList 8] [lindex $fracHingeSecLftList 9] [lindex $fracHingeSecLftList 10] [lindex $fracHingeSecLftList 11]\
		[lindex $fracHingeSecLftList 12] [lindex $fracHingeSecLftList 13] [lindex $fracHingeSecLftList 14] [lindex $fracHingeSecLftList 15];

# right hinge section for beam
# uniaxialMaterial Steel01 [lindex $flangeMatRgtTopList 0] [lindex $flangeMatRgtTopList 1] [lindex $flangeMatRgtTopList 2] [lindex $flangeMatRgtTopList 3];

uniaxialMaterial SteelFractureDI [lindex $flangeMatRgtTopList 0] [lindex $flangeMatRgtTopList 1] [lindex $flangeMatRgtTopList 2] [lindex $flangeMatRgtTopList 3] \
		[lindex $flangeMatRgtTopList 4] [lindex $flangeMatRgtTopList 5] [lindex $flangeMatRgtTopList 6] [lindex $flangeMatRgtTopList 7] \
		[lindex $flangeMatRgtTopList 8] [lindex $flangeMatRgtTopList 9] [lindex $flangeMatRgtTopList 10] [lindex $flangeMatRgtTopList 11] \
		[lindex $flangeMatRgtTopList 12];

uniaxialMaterial SteelFractureDI [lindex $flangeMatRgtBotList 0] [lindex $flangeMatRgtBotList 1] [lindex $flangeMatRgtBotList 2] [lindex $flangeMatRgtBotList 3] \
		[lindex $flangeMatRgtBotList 4] [lindex $flangeMatRgtBotList 5] [lindex $flangeMatRgtBotList 6] [lindex $flangeMatRgtBotList 7] \
		[lindex $flangeMatRgtBotList 8] [lindex $flangeMatRgtBotList 9] [lindex $flangeMatRgtBotList 10] [lindex $flangeMatRgtBotList 11] \
		[lindex $flangeMatRgtBotList 12];
		
fracHingeSection [lindex $fracHingeSecRgtList 0] [lindex $fracHingeSecRgtList 1] [lindex $fracHingeSecRgtList 2] [lindex $fracHingeSecRgtList 3]\
		[lindex $fracHingeSecRgtList 4] [lindex $fracHingeSecRgtList 5] [lindex $fracHingeSecRgtList 6] [lindex $fracHingeSecRgtList 7]\
		[lindex $fracHingeSecRgtList 8] [lindex $fracHingeSecRgtList 9] [lindex $fracHingeSecRgtList 10] [lindex $fracHingeSecRgtList 11]\
		[lindex $fracHingeSecRgtList 12] [lindex $fracHingeSecRgtList 13] [lindex $fracHingeSecRgtList 14] [lindex $fracHingeSecRgtList 15];

# create element
set elSecID [lindex $elSecList 0];
set hingeSecLftID [lindex $fracHingeSecLftList 0];
set hingeSecRgtID [lindex $fracHingeSecRgtList 0];

# # set integration "HingeMidpoint $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
# set integration "HingeRadau $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
# element forceBeamColumn $eleTag $node1 $node2 $transfTag $integration;

# if {[lsearch -exact $inConvEleList $eleTag] >= 0} {
# 	set integration "HingeRadau $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
# 	element forceBeamColumn $eleTag $node1 $node2 $transfTag $integration -iter 50 $convTol;
# } else {
# 	set integration "HingeRadau $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
# 	element forceBeamColumn $eleTag $node1 $node2 $transfTag $integration;
# }

set integration "HingeRadau $hingeSecLftID $hingeLgth $hingeSecRgtID $hingeLgth $elSecID";
element forceBeamColumn $eleTag $node1 $node2 $transfTag $integration -iter 50 $elemConvTol;

}

