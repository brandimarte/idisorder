#!/bin/bash

red1="./reduceBiasPts.sh 500 504"
red3="./reduceEnergyPts.sh 250 504"

${red1} ExVxIel.tot > ExVxIel.red1
${red1} ExVxItot.tot > ExVxItot.red1
${red1} ExVxIsy.tot > ExVxIsy.red1
${red1} ExVxIasy.tot > ExVxIasy.red1
${red3} ExVxIel.red1 > ExVxIel.red3
${red3} ExVxItot.red1 > ExVxItot.red3
${red3} ExVxIsy.red1 > ExVxIsy.red3
${red3} ExVxIasy.red1 > ExVxIasy.red3

${red1} ExVxdIel.tot > ExVxdIel.red1
${red1} ExVxdItot.tot > ExVxdItot.red1
${red1} ExVxdIsy.tot > ExVxdIsy.red1
${red1} ExVxdIasy.tot > ExVxdIasy.red1
${red3} ExVxdIel.red1 > ExVxdIel.red3
${red3} ExVxdItot.red1 > ExVxdItot.red3
${red3} ExVxdIsy.red1 > ExVxdIsy.red3
${red3} ExVxdIasy.red1 > ExVxdIasy.red3

${red1} ExVxd2Iel.tot > ExVxd2Iel.red1
${red1} ExVxd2Itot.tot > ExVxd2Itot.red1
${red1} ExVxd2Isy.tot > ExVxd2Isy.red1
${red1} ExVxd2Iasy.tot > ExVxd2Iasy.red1
${red3} ExVxd2Iel.red1 > ExVxd2Iel.red3
${red3} ExVxd2Itot.red1 > ExVxd2Itot.red3
${red3} ExVxd2Isy.red1 > ExVxd2Isy.red3
${red3} ExVxd2Iasy.red1 > ExVxd2Iasy.red3

./IETSnoGap.sh

${red1} IETSel.tot > IETSel.red1
${red1} IETStot.tot > IETStot.red1
${red1} IETSsy.tot > IETSsy.red1
${red1} IETSasy.tot > IETSasy.red1
${red3} IETSel.red1 > IETSel.red3
${red3} IETStot.red1 > IETStot.red3
${red3} IETSsy.red1 > IETSsy.red3
${red3} IETSasy.red1 > IETSasy.red3

