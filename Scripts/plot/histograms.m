(* ******************************************************************* *)
(* I-Disorder Fortran Code 2007-2014                                   *)
(*                                                                     *)
(* Written by Pedro Brandimarte (brandimarte@gmail.com),               *)
(*            Alberto Torres (alberto.trj@gmail.com) and               *)
(*            Alexandre Reily Rocha (reilya@ift.unesp.br).             *)
(*                                                                     *)
(* Copyright (c), All Rights Reserved                                  *)
(*                                                                     *)
(* This program is free software. You can redistribute it and/or       *)
(* modify it under the terms of the GNU General Public License         *)
(* (version 3 or later) as published by the Free Software Foundation   *)
(* <http://fsf.org/>.                                                  *)
(*                                                                     *)
(* This program is distributed in the hope that it will be useful, but *)
(* WITHOUT ANY WARRANTY, without even the implied warranty of          *)
(* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU    *)
(* General Public License for more details (file 'LICENSE_GPL'         *)
(* distributed along with this program or at                           *)
(* <http://www.gnu.org/licenses/gpl.html>).                            *)
(* ******************************************************************* *)
(*                            histograms.m                             *)
(* ******************************************************************* *)
(* Description: given the '*.dat' files generated by 'genHist.sh'      *)
(* script, it generates the histograms and stores them as jpeg images. *)
(*                                                                     *)
(* Requirements: Wolfram Mathematica software                          *)
(*                                                                     *)
(* How to run: [linux]  > math -script histograms.m                    *)
(*             [Mac OS] > MathKernel -script histograms.m              *)
(*                                                                     *)
(* Written by Pedro Brandimarte, Feb 2014.                             *)
(* Instituto de Fisica                                                 *)
(* Universidade de Sao Paulo                                           *)
(* e-mail: brandimarte@gmail.com                                       *)
(* ***************************** HISTORY ***************************** *)
(* Original version:    February 2014                                  *)
(* ******************************************************************* *)

(* Working directory (containing the '.dat' files). *) 
(* SetDirectory["../../teste/disorder/SimpleDefect05/histo"] *)

files = FileNames["*el.dat"]
For [i = 1, i < Length[files], i = i + 5, 
 Export[StringJoin[StringDrop[files[[i]], -3], "jpg"], 
  Histogram[Import[files[[i]], "List"]], ImageResolution -> 200]]

files = FileNames["*sym.dat"]
For [i = 1, i < Length[files], i = i + 5, 
 Export[StringJoin[StringDrop[files[[i]], -3], "jpg"], 
  Histogram[Import[files[[i]], "List"]], ImageResolution -> 200]]

files = FileNames["*asy.dat"]
For [i = 1, i < Length[files], i = i + 5, 
 Export[StringJoin[StringDrop[files[[i]], -3], "jpg"], 
  Histogram[Import[files[[i]], "List"]], ImageResolution -> 200]]

files = FileNames["*tot.dat"]
For [i = 1, i < Length[files], i = i + 5, 
 Export[StringJoin[StringDrop[files[[i]], -3], "jpg"], 
  Histogram[Import[files[[i]], "List"]], ImageResolution -> 200]]

