file mkdir $recorderdir/beams;

 # fracture recorder
recorder Element -file $recorderdir/beams/frac_LB.out -time -ele 1020100 1020200 1020300 1020400 1020500 1030100 1030200 1030300 1030400 1030500 1040100 1040200 1040300 1040400 1040500 1050100 1050200 1050300 1050400 1050500 1060100 1060200 1060300 1060400 1060500 1070100 1070200 1070300 1070400 1070500 1080100 1080200 1080300 1080400 1080500 1090100 1090200 1090300 1090400 1090500 1100100 1100200 1100300 1100400 1100500 1110100 1110200 1110300 1110400 1110500 1120100 1120200 1120300 1120400 1120500 1130100 1130200 1130300 1130400 1130500 1140100 1140200 1140300 1140400 1140500 1150100 1150200 1150300 1150400 1150500 1160100 1160200 1160300 1160400 1160500 1170100 1170200 1170300 1170400 1170500 1180100 1180200 1180300 1180400 1180500 1190100 1190200 1190300 1190400 1190500 1200100 1200200 1200300 1200400 1200500 1210100 1210200 1210300 1210400 1210500 1220100 1220200 1220300 1220400 1220500 1230100 1230200 1230300 1230400 1230500 1240100 1240200 1240300 1240400 1240500 1250100 1250200 1250300 1250400 1250500 1260100 1260200 1260300 1260400 1260500 1270100 1270200 1270300 1270400 1270500 1280100 1280200 1280300 1280400 1280500 1290100 1290200 1290300 1290400 1290500 1300100 1300200 1300300 1300400 1300500 1310100 1310200 1310300 1310400 1310500 1320100 1320200 1320300 1320400 1320500 1330100 1330200 1330300 1330400 1330500 1340100 1340200 1340300 1340400 1340500 1350100 1350200 1350300 1350400 1350500 1360100 1360200 1360300 1360400 1360500 section 1 fiber -30.0 0.0 failure;
recorder Element -file $recorderdir/beams/frac_LT.out -time -ele 1020100 1020200 1020300 1020400 1020500 1030100 1030200 1030300 1030400 1030500 1040100 1040200 1040300 1040400 1040500 1050100 1050200 1050300 1050400 1050500 1060100 1060200 1060300 1060400 1060500 1070100 1070200 1070300 1070400 1070500 1080100 1080200 1080300 1080400 1080500 1090100 1090200 1090300 1090400 1090500 1100100 1100200 1100300 1100400 1100500 1110100 1110200 1110300 1110400 1110500 1120100 1120200 1120300 1120400 1120500 1130100 1130200 1130300 1130400 1130500 1140100 1140200 1140300 1140400 1140500 1150100 1150200 1150300 1150400 1150500 1160100 1160200 1160300 1160400 1160500 1170100 1170200 1170300 1170400 1170500 1180100 1180200 1180300 1180400 1180500 1190100 1190200 1190300 1190400 1190500 1200100 1200200 1200300 1200400 1200500 1210100 1210200 1210300 1210400 1210500 1220100 1220200 1220300 1220400 1220500 1230100 1230200 1230300 1230400 1230500 1240100 1240200 1240300 1240400 1240500 1250100 1250200 1250300 1250400 1250500 1260100 1260200 1260300 1260400 1260500 1270100 1270200 1270300 1270400 1270500 1280100 1280200 1280300 1280400 1280500 1290100 1290200 1290300 1290400 1290500 1300100 1300200 1300300 1300400 1300500 1310100 1310200 1310300 1310400 1310500 1320100 1320200 1320300 1320400 1320500 1330100 1330200 1330300 1330400 1330500 1340100 1340200 1340300 1340400 1340500 1350100 1350200 1350300 1350400 1350500 1360100 1360200 1360300 1360400 1360500 section 1 fiber 30.0 0.0 failure;
recorder Element -file $recorderdir/beams/frac_RB.out -time -ele 1020100 1020200 1020300 1020400 1020500 1030100 1030200 1030300 1030400 1030500 1040100 1040200 1040300 1040400 1040500 1050100 1050200 1050300 1050400 1050500 1060100 1060200 1060300 1060400 1060500 1070100 1070200 1070300 1070400 1070500 1080100 1080200 1080300 1080400 1080500 1090100 1090200 1090300 1090400 1090500 1100100 1100200 1100300 1100400 1100500 1110100 1110200 1110300 1110400 1110500 1120100 1120200 1120300 1120400 1120500 1130100 1130200 1130300 1130400 1130500 1140100 1140200 1140300 1140400 1140500 1150100 1150200 1150300 1150400 1150500 1160100 1160200 1160300 1160400 1160500 1170100 1170200 1170300 1170400 1170500 1180100 1180200 1180300 1180400 1180500 1190100 1190200 1190300 1190400 1190500 1200100 1200200 1200300 1200400 1200500 1210100 1210200 1210300 1210400 1210500 1220100 1220200 1220300 1220400 1220500 1230100 1230200 1230300 1230400 1230500 1240100 1240200 1240300 1240400 1240500 1250100 1250200 1250300 1250400 1250500 1260100 1260200 1260300 1260400 1260500 1270100 1270200 1270300 1270400 1270500 1280100 1280200 1280300 1280400 1280500 1290100 1290200 1290300 1290400 1290500 1300100 1300200 1300300 1300400 1300500 1310100 1310200 1310300 1310400 1310500 1320100 1320200 1320300 1320400 1320500 1330100 1330200 1330300 1330400 1330500 1340100 1340200 1340300 1340400 1340500 1350100 1350200 1350300 1350400 1350500 1360100 1360200 1360300 1360400 1360500 section 6 fiber -30.0 0.0 failure;
recorder Element -file $recorderdir/beams/frac_RT.out -time -ele 1020100 1020200 1020300 1020400 1020500 1030100 1030200 1030300 1030400 1030500 1040100 1040200 1040300 1040400 1040500 1050100 1050200 1050300 1050400 1050500 1060100 1060200 1060300 1060400 1060500 1070100 1070200 1070300 1070400 1070500 1080100 1080200 1080300 1080400 1080500 1090100 1090200 1090300 1090400 1090500 1100100 1100200 1100300 1100400 1100500 1110100 1110200 1110300 1110400 1110500 1120100 1120200 1120300 1120400 1120500 1130100 1130200 1130300 1130400 1130500 1140100 1140200 1140300 1140400 1140500 1150100 1150200 1150300 1150400 1150500 1160100 1160200 1160300 1160400 1160500 1170100 1170200 1170300 1170400 1170500 1180100 1180200 1180300 1180400 1180500 1190100 1190200 1190300 1190400 1190500 1200100 1200200 1200300 1200400 1200500 1210100 1210200 1210300 1210400 1210500 1220100 1220200 1220300 1220400 1220500 1230100 1230200 1230300 1230400 1230500 1240100 1240200 1240300 1240400 1240500 1250100 1250200 1250300 1250400 1250500 1260100 1260200 1260300 1260400 1260500 1270100 1270200 1270300 1270400 1270500 1280100 1280200 1280300 1280400 1280500 1290100 1290200 1290300 1290400 1290500 1300100 1300200 1300300 1300400 1300500 1310100 1310200 1310300 1310400 1310500 1320100 1320200 1320300 1320400 1320500 1330100 1330200 1330300 1330400 1330500 1340100 1340200 1340300 1340400 1340500 1350100 1350200 1350300 1350400 1350500 1360100 1360200 1360300 1360400 1360500 section 6 fiber 30.0 0.0 failure;

 # DI recorder
recorder Element -file $recorderdir/beams/DI_LB.out -time -ele 1020100 1020200 1020300 1020400 1020500 1030100 1030200 1030300 1030400 1030500 1040100 1040200 1040300 1040400 1040500 1050100 1050200 1050300 1050400 1050500 1060100 1060200 1060300 1060400 1060500 1070100 1070200 1070300 1070400 1070500 1080100 1080200 1080300 1080400 1080500 1090100 1090200 1090300 1090400 1090500 1100100 1100200 1100300 1100400 1100500 1110100 1110200 1110300 1110400 1110500 1120100 1120200 1120300 1120400 1120500 1130100 1130200 1130300 1130400 1130500 1140100 1140200 1140300 1140400 1140500 1150100 1150200 1150300 1150400 1150500 1160100 1160200 1160300 1160400 1160500 1170100 1170200 1170300 1170400 1170500 1180100 1180200 1180300 1180400 1180500 1190100 1190200 1190300 1190400 1190500 1200100 1200200 1200300 1200400 1200500 1210100 1210200 1210300 1210400 1210500 1220100 1220200 1220300 1220400 1220500 1230100 1230200 1230300 1230400 1230500 1240100 1240200 1240300 1240400 1240500 1250100 1250200 1250300 1250400 1250500 1260100 1260200 1260300 1260400 1260500 1270100 1270200 1270300 1270400 1270500 1280100 1280200 1280300 1280400 1280500 1290100 1290200 1290300 1290400 1290500 1300100 1300200 1300300 1300400 1300500 1310100 1310200 1310300 1310400 1310500 1320100 1320200 1320300 1320400 1320500 1330100 1330200 1330300 1330400 1330500 1340100 1340200 1340300 1340400 1340500 1350100 1350200 1350300 1350400 1350500 1360100 1360200 1360300 1360400 1360500 section 1 fiber -30.0 0.0 damage;
recorder Element -file $recorderdir/beams/DI_LT.out -time -ele 1020100 1020200 1020300 1020400 1020500 1030100 1030200 1030300 1030400 1030500 1040100 1040200 1040300 1040400 1040500 1050100 1050200 1050300 1050400 1050500 1060100 1060200 1060300 1060400 1060500 1070100 1070200 1070300 1070400 1070500 1080100 1080200 1080300 1080400 1080500 1090100 1090200 1090300 1090400 1090500 1100100 1100200 1100300 1100400 1100500 1110100 1110200 1110300 1110400 1110500 1120100 1120200 1120300 1120400 1120500 1130100 1130200 1130300 1130400 1130500 1140100 1140200 1140300 1140400 1140500 1150100 1150200 1150300 1150400 1150500 1160100 1160200 1160300 1160400 1160500 1170100 1170200 1170300 1170400 1170500 1180100 1180200 1180300 1180400 1180500 1190100 1190200 1190300 1190400 1190500 1200100 1200200 1200300 1200400 1200500 1210100 1210200 1210300 1210400 1210500 1220100 1220200 1220300 1220400 1220500 1230100 1230200 1230300 1230400 1230500 1240100 1240200 1240300 1240400 1240500 1250100 1250200 1250300 1250400 1250500 1260100 1260200 1260300 1260400 1260500 1270100 1270200 1270300 1270400 1270500 1280100 1280200 1280300 1280400 1280500 1290100 1290200 1290300 1290400 1290500 1300100 1300200 1300300 1300400 1300500 1310100 1310200 1310300 1310400 1310500 1320100 1320200 1320300 1320400 1320500 1330100 1330200 1330300 1330400 1330500 1340100 1340200 1340300 1340400 1340500 1350100 1350200 1350300 1350400 1350500 1360100 1360200 1360300 1360400 1360500 section 1 fiber 30.0 0.0 damage;
recorder Element -file $recorderdir/beams/DI_RB.out -time -ele 1020100 1020200 1020300 1020400 1020500 1030100 1030200 1030300 1030400 1030500 1040100 1040200 1040300 1040400 1040500 1050100 1050200 1050300 1050400 1050500 1060100 1060200 1060300 1060400 1060500 1070100 1070200 1070300 1070400 1070500 1080100 1080200 1080300 1080400 1080500 1090100 1090200 1090300 1090400 1090500 1100100 1100200 1100300 1100400 1100500 1110100 1110200 1110300 1110400 1110500 1120100 1120200 1120300 1120400 1120500 1130100 1130200 1130300 1130400 1130500 1140100 1140200 1140300 1140400 1140500 1150100 1150200 1150300 1150400 1150500 1160100 1160200 1160300 1160400 1160500 1170100 1170200 1170300 1170400 1170500 1180100 1180200 1180300 1180400 1180500 1190100 1190200 1190300 1190400 1190500 1200100 1200200 1200300 1200400 1200500 1210100 1210200 1210300 1210400 1210500 1220100 1220200 1220300 1220400 1220500 1230100 1230200 1230300 1230400 1230500 1240100 1240200 1240300 1240400 1240500 1250100 1250200 1250300 1250400 1250500 1260100 1260200 1260300 1260400 1260500 1270100 1270200 1270300 1270400 1270500 1280100 1280200 1280300 1280400 1280500 1290100 1290200 1290300 1290400 1290500 1300100 1300200 1300300 1300400 1300500 1310100 1310200 1310300 1310400 1310500 1320100 1320200 1320300 1320400 1320500 1330100 1330200 1330300 1330400 1330500 1340100 1340200 1340300 1340400 1340500 1350100 1350200 1350300 1350400 1350500 1360100 1360200 1360300 1360400 1360500 section 6 fiber -30.0 0.0 damage;
recorder Element -file $recorderdir/beams/DI_RT.out -time -ele 1020100 1020200 1020300 1020400 1020500 1030100 1030200 1030300 1030400 1030500 1040100 1040200 1040300 1040400 1040500 1050100 1050200 1050300 1050400 1050500 1060100 1060200 1060300 1060400 1060500 1070100 1070200 1070300 1070400 1070500 1080100 1080200 1080300 1080400 1080500 1090100 1090200 1090300 1090400 1090500 1100100 1100200 1100300 1100400 1100500 1110100 1110200 1110300 1110400 1110500 1120100 1120200 1120300 1120400 1120500 1130100 1130200 1130300 1130400 1130500 1140100 1140200 1140300 1140400 1140500 1150100 1150200 1150300 1150400 1150500 1160100 1160200 1160300 1160400 1160500 1170100 1170200 1170300 1170400 1170500 1180100 1180200 1180300 1180400 1180500 1190100 1190200 1190300 1190400 1190500 1200100 1200200 1200300 1200400 1200500 1210100 1210200 1210300 1210400 1210500 1220100 1220200 1220300 1220400 1220500 1230100 1230200 1230300 1230400 1230500 1240100 1240200 1240300 1240400 1240500 1250100 1250200 1250300 1250400 1250500 1260100 1260200 1260300 1260400 1260500 1270100 1270200 1270300 1270400 1270500 1280100 1280200 1280300 1280400 1280500 1290100 1290200 1290300 1290400 1290500 1300100 1300200 1300300 1300400 1300500 1310100 1310200 1310300 1310400 1310500 1320100 1320200 1320300 1320400 1320500 1330100 1330200 1330300 1330400 1330500 1340100 1340200 1340300 1340400 1340500 1350100 1350200 1350300 1350400 1350500 1360100 1360200 1360300 1360400 1360500 section 6 fiber 30.0 0.0 damage;

 # F-D recorder
recorder Element -file $recorderdir/beams/D_L.out -time -ele 1020100 1020200 1020300 1020400 1020500 1030100 1030200 1030300 1030400 1030500 1040100 1040200 1040300 1040400 1040500 1050100 1050200 1050300 1050400 1050500 1060100 1060200 1060300 1060400 1060500 1070100 1070200 1070300 1070400 1070500 1080100 1080200 1080300 1080400 1080500 1090100 1090200 1090300 1090400 1090500 1100100 1100200 1100300 1100400 1100500 1110100 1110200 1110300 1110400 1110500 1120100 1120200 1120300 1120400 1120500 1130100 1130200 1130300 1130400 1130500 1140100 1140200 1140300 1140400 1140500 1150100 1150200 1150300 1150400 1150500 1160100 1160200 1160300 1160400 1160500 1170100 1170200 1170300 1170400 1170500 1180100 1180200 1180300 1180400 1180500 1190100 1190200 1190300 1190400 1190500 1200100 1200200 1200300 1200400 1200500 1210100 1210200 1210300 1210400 1210500 1220100 1220200 1220300 1220400 1220500 1230100 1230200 1230300 1230400 1230500 1240100 1240200 1240300 1240400 1240500 1250100 1250200 1250300 1250400 1250500 1260100 1260200 1260300 1260400 1260500 1270100 1270200 1270300 1270400 1270500 1280100 1280200 1280300 1280400 1280500 1290100 1290200 1290300 1290400 1290500 1300100 1300200 1300300 1300400 1300500 1310100 1310200 1310300 1310400 1310500 1320100 1320200 1320300 1320400 1320500 1330100 1330200 1330300 1330400 1330500 1340100 1340200 1340300 1340400 1340500 1350100 1350200 1350300 1350400 1350500 1360100 1360200 1360300 1360400 1360500 section 1 deformation;
recorder Element -file $recorderdir/beams/F_L.out -time -ele 1020100 1020200 1020300 1020400 1020500 1030100 1030200 1030300 1030400 1030500 1040100 1040200 1040300 1040400 1040500 1050100 1050200 1050300 1050400 1050500 1060100 1060200 1060300 1060400 1060500 1070100 1070200 1070300 1070400 1070500 1080100 1080200 1080300 1080400 1080500 1090100 1090200 1090300 1090400 1090500 1100100 1100200 1100300 1100400 1100500 1110100 1110200 1110300 1110400 1110500 1120100 1120200 1120300 1120400 1120500 1130100 1130200 1130300 1130400 1130500 1140100 1140200 1140300 1140400 1140500 1150100 1150200 1150300 1150400 1150500 1160100 1160200 1160300 1160400 1160500 1170100 1170200 1170300 1170400 1170500 1180100 1180200 1180300 1180400 1180500 1190100 1190200 1190300 1190400 1190500 1200100 1200200 1200300 1200400 1200500 1210100 1210200 1210300 1210400 1210500 1220100 1220200 1220300 1220400 1220500 1230100 1230200 1230300 1230400 1230500 1240100 1240200 1240300 1240400 1240500 1250100 1250200 1250300 1250400 1250500 1260100 1260200 1260300 1260400 1260500 1270100 1270200 1270300 1270400 1270500 1280100 1280200 1280300 1280400 1280500 1290100 1290200 1290300 1290400 1290500 1300100 1300200 1300300 1300400 1300500 1310100 1310200 1310300 1310400 1310500 1320100 1320200 1320300 1320400 1320500 1330100 1330200 1330300 1330400 1330500 1340100 1340200 1340300 1340400 1340500 1350100 1350200 1350300 1350400 1350500 1360100 1360200 1360300 1360400 1360500 section 1 force;
recorder Element -file $recorderdir/beams/D_R.out -time -ele 1020100 1020200 1020300 1020400 1020500 1030100 1030200 1030300 1030400 1030500 1040100 1040200 1040300 1040400 1040500 1050100 1050200 1050300 1050400 1050500 1060100 1060200 1060300 1060400 1060500 1070100 1070200 1070300 1070400 1070500 1080100 1080200 1080300 1080400 1080500 1090100 1090200 1090300 1090400 1090500 1100100 1100200 1100300 1100400 1100500 1110100 1110200 1110300 1110400 1110500 1120100 1120200 1120300 1120400 1120500 1130100 1130200 1130300 1130400 1130500 1140100 1140200 1140300 1140400 1140500 1150100 1150200 1150300 1150400 1150500 1160100 1160200 1160300 1160400 1160500 1170100 1170200 1170300 1170400 1170500 1180100 1180200 1180300 1180400 1180500 1190100 1190200 1190300 1190400 1190500 1200100 1200200 1200300 1200400 1200500 1210100 1210200 1210300 1210400 1210500 1220100 1220200 1220300 1220400 1220500 1230100 1230200 1230300 1230400 1230500 1240100 1240200 1240300 1240400 1240500 1250100 1250200 1250300 1250400 1250500 1260100 1260200 1260300 1260400 1260500 1270100 1270200 1270300 1270400 1270500 1280100 1280200 1280300 1280400 1280500 1290100 1290200 1290300 1290400 1290500 1300100 1300200 1300300 1300400 1300500 1310100 1310200 1310300 1310400 1310500 1320100 1320200 1320300 1320400 1320500 1330100 1330200 1330300 1330400 1330500 1340100 1340200 1340300 1340400 1340500 1350100 1350200 1350300 1350400 1350500 1360100 1360200 1360300 1360400 1360500 section 6 deformation;
recorder Element -file $recorderdir/beams/F_R.out -time -ele 1020100 1020200 1020300 1020400 1020500 1030100 1030200 1030300 1030400 1030500 1040100 1040200 1040300 1040400 1040500 1050100 1050200 1050300 1050400 1050500 1060100 1060200 1060300 1060400 1060500 1070100 1070200 1070300 1070400 1070500 1080100 1080200 1080300 1080400 1080500 1090100 1090200 1090300 1090400 1090500 1100100 1100200 1100300 1100400 1100500 1110100 1110200 1110300 1110400 1110500 1120100 1120200 1120300 1120400 1120500 1130100 1130200 1130300 1130400 1130500 1140100 1140200 1140300 1140400 1140500 1150100 1150200 1150300 1150400 1150500 1160100 1160200 1160300 1160400 1160500 1170100 1170200 1170300 1170400 1170500 1180100 1180200 1180300 1180400 1180500 1190100 1190200 1190300 1190400 1190500 1200100 1200200 1200300 1200400 1200500 1210100 1210200 1210300 1210400 1210500 1220100 1220200 1220300 1220400 1220500 1230100 1230200 1230300 1230400 1230500 1240100 1240200 1240300 1240400 1240500 1250100 1250200 1250300 1250400 1250500 1260100 1260200 1260300 1260400 1260500 1270100 1270200 1270300 1270400 1270500 1280100 1280200 1280300 1280400 1280500 1290100 1290200 1290300 1290400 1290500 1300100 1300200 1300300 1300400 1300500 1310100 1310200 1310300 1310400 1310500 1320100 1320200 1320300 1320400 1320500 1330100 1330200 1330300 1330400 1330500 1340100 1340200 1340300 1340400 1340500 1350100 1350200 1350300 1350400 1350500 1360100 1360200 1360300 1360400 1360500 section 6 force;

 # stress-strain recorder
recorder Element -file $recorderdir/beams/SS_LB.out -time -ele 1020100 1020200 1020300 1020400 1020500 1030100 1030200 1030300 1030400 1030500 1040100 1040200 1040300 1040400 1040500 1050100 1050200 1050300 1050400 1050500 1060100 1060200 1060300 1060400 1060500 1070100 1070200 1070300 1070400 1070500 1080100 1080200 1080300 1080400 1080500 1090100 1090200 1090300 1090400 1090500 1100100 1100200 1100300 1100400 1100500 1110100 1110200 1110300 1110400 1110500 1120100 1120200 1120300 1120400 1120500 1130100 1130200 1130300 1130400 1130500 1140100 1140200 1140300 1140400 1140500 1150100 1150200 1150300 1150400 1150500 1160100 1160200 1160300 1160400 1160500 1170100 1170200 1170300 1170400 1170500 1180100 1180200 1180300 1180400 1180500 1190100 1190200 1190300 1190400 1190500 1200100 1200200 1200300 1200400 1200500 1210100 1210200 1210300 1210400 1210500 1220100 1220200 1220300 1220400 1220500 1230100 1230200 1230300 1230400 1230500 1240100 1240200 1240300 1240400 1240500 1250100 1250200 1250300 1250400 1250500 1260100 1260200 1260300 1260400 1260500 1270100 1270200 1270300 1270400 1270500 1280100 1280200 1280300 1280400 1280500 1290100 1290200 1290300 1290400 1290500 1300100 1300200 1300300 1300400 1300500 1310100 1310200 1310300 1310400 1310500 1320100 1320200 1320300 1320400 1320500 1330100 1330200 1330300 1330400 1330500 1340100 1340200 1340300 1340400 1340500 1350100 1350200 1350300 1350400 1350500 1360100 1360200 1360300 1360400 1360500 section 1 fiber -30.0 0.0 stressStrain;
recorder Element -file $recorderdir/beams/SS_LT.out -time -ele 1020100 1020200 1020300 1020400 1020500 1030100 1030200 1030300 1030400 1030500 1040100 1040200 1040300 1040400 1040500 1050100 1050200 1050300 1050400 1050500 1060100 1060200 1060300 1060400 1060500 1070100 1070200 1070300 1070400 1070500 1080100 1080200 1080300 1080400 1080500 1090100 1090200 1090300 1090400 1090500 1100100 1100200 1100300 1100400 1100500 1110100 1110200 1110300 1110400 1110500 1120100 1120200 1120300 1120400 1120500 1130100 1130200 1130300 1130400 1130500 1140100 1140200 1140300 1140400 1140500 1150100 1150200 1150300 1150400 1150500 1160100 1160200 1160300 1160400 1160500 1170100 1170200 1170300 1170400 1170500 1180100 1180200 1180300 1180400 1180500 1190100 1190200 1190300 1190400 1190500 1200100 1200200 1200300 1200400 1200500 1210100 1210200 1210300 1210400 1210500 1220100 1220200 1220300 1220400 1220500 1230100 1230200 1230300 1230400 1230500 1240100 1240200 1240300 1240400 1240500 1250100 1250200 1250300 1250400 1250500 1260100 1260200 1260300 1260400 1260500 1270100 1270200 1270300 1270400 1270500 1280100 1280200 1280300 1280400 1280500 1290100 1290200 1290300 1290400 1290500 1300100 1300200 1300300 1300400 1300500 1310100 1310200 1310300 1310400 1310500 1320100 1320200 1320300 1320400 1320500 1330100 1330200 1330300 1330400 1330500 1340100 1340200 1340300 1340400 1340500 1350100 1350200 1350300 1350400 1350500 1360100 1360200 1360300 1360400 1360500 section 1 fiber 30.0 0.0 stressStrain;
recorder Element -file $recorderdir/beams/SS_RB.out -time -ele 1020100 1020200 1020300 1020400 1020500 1030100 1030200 1030300 1030400 1030500 1040100 1040200 1040300 1040400 1040500 1050100 1050200 1050300 1050400 1050500 1060100 1060200 1060300 1060400 1060500 1070100 1070200 1070300 1070400 1070500 1080100 1080200 1080300 1080400 1080500 1090100 1090200 1090300 1090400 1090500 1100100 1100200 1100300 1100400 1100500 1110100 1110200 1110300 1110400 1110500 1120100 1120200 1120300 1120400 1120500 1130100 1130200 1130300 1130400 1130500 1140100 1140200 1140300 1140400 1140500 1150100 1150200 1150300 1150400 1150500 1160100 1160200 1160300 1160400 1160500 1170100 1170200 1170300 1170400 1170500 1180100 1180200 1180300 1180400 1180500 1190100 1190200 1190300 1190400 1190500 1200100 1200200 1200300 1200400 1200500 1210100 1210200 1210300 1210400 1210500 1220100 1220200 1220300 1220400 1220500 1230100 1230200 1230300 1230400 1230500 1240100 1240200 1240300 1240400 1240500 1250100 1250200 1250300 1250400 1250500 1260100 1260200 1260300 1260400 1260500 1270100 1270200 1270300 1270400 1270500 1280100 1280200 1280300 1280400 1280500 1290100 1290200 1290300 1290400 1290500 1300100 1300200 1300300 1300400 1300500 1310100 1310200 1310300 1310400 1310500 1320100 1320200 1320300 1320400 1320500 1330100 1330200 1330300 1330400 1330500 1340100 1340200 1340300 1340400 1340500 1350100 1350200 1350300 1350400 1350500 1360100 1360200 1360300 1360400 1360500 section 6 fiber -30.0 0.0 stressStrain;
recorder Element -file $recorderdir/beams/SS_RT.out -time -ele 1020100 1020200 1020300 1020400 1020500 1030100 1030200 1030300 1030400 1030500 1040100 1040200 1040300 1040400 1040500 1050100 1050200 1050300 1050400 1050500 1060100 1060200 1060300 1060400 1060500 1070100 1070200 1070300 1070400 1070500 1080100 1080200 1080300 1080400 1080500 1090100 1090200 1090300 1090400 1090500 1100100 1100200 1100300 1100400 1100500 1110100 1110200 1110300 1110400 1110500 1120100 1120200 1120300 1120400 1120500 1130100 1130200 1130300 1130400 1130500 1140100 1140200 1140300 1140400 1140500 1150100 1150200 1150300 1150400 1150500 1160100 1160200 1160300 1160400 1160500 1170100 1170200 1170300 1170400 1170500 1180100 1180200 1180300 1180400 1180500 1190100 1190200 1190300 1190400 1190500 1200100 1200200 1200300 1200400 1200500 1210100 1210200 1210300 1210400 1210500 1220100 1220200 1220300 1220400 1220500 1230100 1230200 1230300 1230400 1230500 1240100 1240200 1240300 1240400 1240500 1250100 1250200 1250300 1250400 1250500 1260100 1260200 1260300 1260400 1260500 1270100 1270200 1270300 1270400 1270500 1280100 1280200 1280300 1280400 1280500 1290100 1290200 1290300 1290400 1290500 1300100 1300200 1300300 1300400 1300500 1310100 1310200 1310300 1310400 1310500 1320100 1320200 1320300 1320400 1320500 1330100 1330200 1330300 1330400 1330500 1340100 1340200 1340300 1340400 1340500 1350100 1350200 1350300 1350400 1350500 1360100 1360200 1360300 1360400 1360500 section 6 fiber 30.0 0.0 stressStrain;
