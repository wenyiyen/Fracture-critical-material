# 

proc PanelZone2D {eleID Es A_pz I_pz transfTag Fy dc bc tcf tp db alpha} {


################################################
# 
# Build rigid links
#
################################################

set floor [expr $eleID / 10000]; # change here accordingly
set pier [expr ($eleID-$floor*10000) / 100]; # change here accordingly
set tag_tmp [expr $floor * 10000 + $pier * 100]; 
set nid_tmp [expr $floor * 10000 + $pier * 100];

# 8 rigid links with large A and I
# eleID nodei nodej A E I transfTag
element elasticBeamColumn [expr $tag_tmp + 1] [expr $nid_tmp + 2] [expr $nid_tmp + 11] $A_pz $Es $I_pz $transfTag;
element elasticBeamColumn [expr $tag_tmp + 2] [expr $nid_tmp + 11] [expr $nid_tmp + 3] $A_pz $Es $I_pz $transfTag;
element elasticBeamColumn [expr $tag_tmp + 3] [expr $nid_tmp + 4] [expr $nid_tmp + 9] $A_pz $Es $I_pz $transfTag;
element elasticBeamColumn [expr $tag_tmp + 4] [expr $nid_tmp + 9] [expr $nid_tmp + 5] $A_pz $Es $I_pz $transfTag;
element elasticBeamColumn [expr $tag_tmp + 5] [expr $nid_tmp + 6] [expr $nid_tmp + 12] $A_pz $Es $I_pz $transfTag;
element elasticBeamColumn [expr $tag_tmp + 6] [expr $nid_tmp + 12] [expr $nid_tmp + 7] $A_pz $Es $I_pz $transfTag;
element elasticBeamColumn [expr $tag_tmp + 7] [expr $nid_tmp + 8] [expr $nid_tmp + 10] $A_pz $Es $I_pz $transfTag;
element elasticBeamColumn [expr $tag_tmp + 8] [expr $nid_tmp + 10] [expr $nid_tmp + 1] $A_pz $Es $I_pz $transfTag;

# pin 3 corners of panel zone (leave one corner for spring element)
equalDOF [expr $nid_tmp + 1] [expr $nid_tmp + 2] 1 2;
equalDOF [expr $nid_tmp + 3] [expr $nid_tmp + 4] 1 2;
equalDOF [expr $nid_tmp + 7] [expr $nid_tmp + 8] 1 2;


################################################
# 
# Build panel zone spring
#
################################################

# Trilinear shear-strain curve (Gupta and Krawinkler, 1999)

set G [expr $Es / (2.0 * (1.0 + 0.3))];
set Vy [expr 0.55 * $Fy * $dc * $tp];
set Ke [expr 0.95 * $dc * $tp * $G];
set Kp [expr 0.95 * $bc * pow($tcf, 2) * $G / $db];

set gamma_y [expr $Vy / $Ke];
set gamma_p [expr 4*$gamma_y];
set gamma_3 [expr 100*$gamma_y]; # third point of curve

# Transition to M-gamma relationship
# K_moment = db * K_shear

set My [expr $gamma_y * ($Ke * $db)];
set Mp [expr $My + ($Kp * $db) * ($gamma_p - $gamma_y)];
set M3 [expr $Mp + ($Ke * $alpha * $db) * ($gamma_3 - $gamma_p)];

# Nodes at the spring element
set node_R [expr $nid_tmp + 6];
set node_C [expr $nid_tmp + 5];


# material with no pinching and damage, matID same as eleID

uniaxialMaterial Hysteretic $eleID $My $gamma_y $Mp $gamma_p $M3 $gamma_3 \
		[expr -$My] [expr -$gamma_y] [expr -$Mp] [expr -$gamma_p] [expr -$M3] [expr -$gamma_3] \
		1 1 0.0 0.0 0.0;

# Build spring element
element zeroLength $eleID $node_R $node_C -mat $eleID -dir 6

# constrain translation movement
equalDOF $node_R $node_C 1 2



}