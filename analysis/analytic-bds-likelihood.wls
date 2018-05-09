#!/usr/bin/wolframscript 
(* ::Package:: *)

(* ::Section:: *)
(*Check that PGF satisfies PDE*)


H=((\[Alpha]1 - \[Alpha]2 * ((s - \[Alpha]1) / (s - \[Alpha]2)) * Exp[-\[Lambda] * (\[Alpha]2 - \[Alpha]1) * u * t])/ (1 - (s - \[Alpha]1) / (s - \[Alpha]2) * Exp[-\[Lambda] * (\[Alpha]2 - \[Alpha]1) * u * t]))^i *((\[Alpha]1 - \[Alpha]2)/(s - \[Alpha]2 - (s - \[Alpha]1) * Exp[-\[Lambda] * (\[Alpha]2 - \[Alpha]1) * u * t]))^(\[Nu] / \[Lambda]) * Exp[-\[Nu] * (1 - u * \[Alpha]1) * t]
Hnoi=(\[Alpha]1 - \[Alpha]2 * ((s - \[Alpha]1) / (s - \[Alpha]2)) * Exp[-\[Lambda] * (\[Alpha]2 - \[Alpha]1) * u * t])/ (1 - (s - \[Alpha]1) / (s - \[Alpha]2) * Exp[-\[Lambda] * (\[Alpha]2 - \[Alpha]1) * u * t])
D[Hnoi, t] //FullSimplify


D[Hnoi, s] //FullSimplify


foo =D[Hnoi, t] / D[Hnoi, s] //FullSimplify //Expand


alphas = Table[(\[Lambda]+\[Mu]+w+(-1)^x * Sqrt[(\[Lambda]+\[Mu]+w)^ 2\[Minus]4*\[Lambda]*\[Mu]*u*v]) / (2 *\[Lambda]*u), {x, 2}]


repAlpha = {\[Alpha]1 -> alphas[[1]], \[Alpha]2 -> alphas[[2]]}
Print[foo /. repAlpha //FullSimplify //Expand]


(* ::Text:: *)
(*		Agrees  with S-1 in Doss et al. (http://dx.doi.org/10.1214%2F13-AOAS673) :*)


(* ::DisplayFormula:: *)
(*s^2 u \[Lambda] - (\[Lambda] + \[Mu] + w) s + v \[Mu]*)


(* ::Section:: *)
(*Calculation of expected area under population size curve*)
