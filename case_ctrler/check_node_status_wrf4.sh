#/bin/sh 

NODE_LST=("hqlx47" "hqlx100" "hqlx111" "hqlx133" "hqlx132" "hqlx62" "hqlx65" "hqlx69" "hqlx129" "hqlx130")
QUE_LST=("P1" "P2")
for NODE in ${NODE_LST[@]}
do
    #FULL_NAME=${NODE//hq/hqlx}
    echo ">>>>>>>>>>>>>>>>>>>"$NODE" START<<<<<<<<<<<<<<<<<<<<<<<<<"
    for QUE in ${QUE_LST[@]}
    do
        WRF_DIR=/home/lzhenn/array${NODE:4}/WRF-4.3_${QUE}/run
        ls -l $WRF_DIR/wrfbdy* | tail -1
        ls -l $WRF_DIR/wrfout_d04* | tail -1
        ls -l $WRF_DIR/rsl.out.0000 
        tail -5 $WRF_DIR/rsl.out.0000
    done
    ssh $NODE "quota -uvs lzhenn"
    echo ">>>>>>>>>>>>>>>>>>>"$NODE" END<<<<<<<<<<<<<<<<<<<<<<<<<"
    echo ""
    echo ""
done

