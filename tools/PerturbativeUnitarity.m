(* ::Package:: *)
BeginPackage["PerturbativeUnitarity`"]

scatterMatrix::usage = "Return the full scatter matrix of the states constructed from fields using the interactions from the scalar potential V. The second element contains the states corresponding to the rows, and the third element the rules that charge conjugate a state.";
scatterMatrixForm::usage ="nicely print a scattering matrix using the states as row and column labels";
splitBlockMatrix::usage = "Split the given scatter matrix into block matrices keeping track of the corresponding states";
uniqueEV::usage = "Computes the unique Eigenvalues of the given scatter block matrices";
maxUnitaryParRange::usage = "numerically computes the maximum possible range for each parameter without violating the unitarity constraints from the eigenvalues ev";

Begin["`Private`"]

chargeConjugationRules[fields_List]:=Flatten[#/.{f_List:>{Conjugate[f[[1]]]->f[[2]],Conjugate[f[[2]]]->f[[1]]},x_:>Conjugate[x]->x}&/@fields]
chargeConjugationRules::usage = "get the rules that charge conjugate the fields";

chargeRules[fields_List]:=Flatten[#/.{f_List:>{f[[1]]->1,f[[2]]->-1},x_:>x->0}&/@fields]
chargeRules::usage = "map each field to its electric charge";

generateStates[fields_List]:=Module[{cRules},
cRules = chargeRules[fields];
{Map[{If[#[[1]]===#[[2]],1/Sqrt[2],1],#}&,Select[DeleteDuplicates[Map[Sort,Tuples[Flatten[fields],2]]],Total[#/.cRules]>=0&]],
chargeConjugationRules[fields]}
]
generateStates::usage = "Returns a list of appropriately normalized scattering states and the rules that charge conjugate a state.";

combStates[s1_,s2_]:={s1[[1]]*s2[[1]],Flatten[{s1[[2]],s2[[2]]}]}
combStates::usage = "combine two states";

scatterMatrix[V_,fields_]:=Module[{states,conjrule},
{states ,conjrule}= generateStates[fields];
{Map[ Simplify[#[[1]]* D[V,Sequence@@#[[2]]]]&,Outer[combStates,states,Conjugate[states]/.conjrule,1],{2}],
states,
conjrule}
]

scatterMatrixForm[sMat_List]:=Module[{matrix,states,conjrule},
{matrix,states,conjrule}=sMat;
TableForm[matrix,TableHeadings->{Times@@@Flatten/@states,Times@@@Flatten/@(Conjugate[states]/.conjrule)}]
]

findBlock[x_List,matrix_]:=Union[Flatten[({#,SparseArray[matrix[[#]]]["NonzeroPositions"]})&/@x]]
findBlock::usage ="finds the indices of the smallest sub-block of the matrix that contains the xth rows and columns";

splitBlockMatrix[smat_]:=Module[{blocks,matrix,states,conjrule},
{matrix,states,conjrule}=smat;
blocks=Union[FixedPoint[findBlock[#,matrix]&,{#},4]&/@Range[Length[matrix]]];
{Part[matrix,#,#],Part[states,#],conjrule}&/@blocks
]

uniqueEV[blockSmat_]:=Union[Flatten[FullSimplify[Eigenvalues[#[[1]]]&/@blockSmat,Reals]]]



maxUnitaryParRange[ev_,params_,additionalConstraints_:True] := Module[{constraints},
constraints = And[And@@(Abs[#]<8 Pi &/@ev),additionalConstraints];
({#,
{NMinimize[{#,constraints},params,MaxIterations->1000][[1]],
NMaximize[{#,constraints},params,MaxIterations->1000][[1]]}}&/@params)]

End[]

EndPackage[]
