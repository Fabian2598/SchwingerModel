#!/bin/bash
set -euo pipefail

for N in 64 128 256; do
    for BETA in 2; do
        NCONF=0

        if [ "$BETA" -eq 2 ]; then
            if [ "$N" -eq 128 ]; then
                NCONF=3
            elif [ "$N" -eq 256 ]; then
                NCONF=20
            fi
            MASS="-018"
            MASS2="-01884"
            M0=-0.18840579710144945

        elif [ "$BETA" -eq 4 ]; then
            MASS="-01023"
            MASS2="-01023"
            M0=-0.1023

        elif [ "$BETA" -eq 6 ]; then
            MASS="-0709"
            MASS2="-00709"
            M0=-0.0709
        fi

        CONF_PATH="confs/b${BETA}_${N}x${N}/m${MASS}/2D_U1_Ns${N}_Nt${N}_b${BETA}0000_m${MASS2}_${NCONF}.ctxt"
        printf "%s\n%s\n%s\n%s\n" "${M0}" "${BETA}" "${CONF_PATH}" "${NCONF}" >> inputs
        echo "${CONF_PATH}"
        ./dirac${N} < inputs
        rm inputs
    done
done