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


-D[Hnoi/. repAlpha /.{v->1, s->1, u->1}, w]/.{w->0}//FullSimplify


Print[prob = Hnoi /. repAlpha /. {v->1,s->1,u->1} //FullSimplify]


prob /. {\[Lambda]->1, \[Mu]->.2, t->3}


D[prob, t] + (\[Lambda] + \[Mu] + w) * prob - \[Lambda] * prob^2 -\[Mu]// FullSimplify


(* agrees with differential equation (2) in Kuhnert et al. *)


(* ::Input:: *)
(*altAns = DSolve[{p'[t]== -(\[Lambda] + \[Mu] + w) * p[t] - \[Lambda] p[t]^2 + mu, p[0]==1}, p[t],t]*)


Print[altAns[[1]][[1]]//FullSimplify]


eqs = {p1'[t] == -(\[Lambda]11 + \[Lambda]12 + m12 + \[Mu]1 + w1) * p1[t] + \[Lambda]11 * p1[t] * p1[t] + \[Lambda]12 * p1[t] * p2[t] + m21 * p2[t] + \[Mu]1,p2'[t] == -(\[Lambda]21 + \[Lambda]22 + m21 + \[Mu]2 + w2) * p2[t] + \[Lambda]21 * p1[t] * p2[t] + \[Lambda]22 * p2[t] * p2[t] + m12 * p1[t] + \[Mu]2}


DSolve[eqs, {p1[t], p2[t]}, t]
