(* Content-type: application/vnd.wolfram.cdf.text *)

(*** Wolfram CDF File ***)
(* http://www.wolfram.com/cdf *)

(* CreatedBy='Mathematica 11.2' *)

(***************************************************************************)
(*                                                                         *)
(*                                                                         *)
(*  Under the Wolfram FreeCDF terms of use, this file and its content are  *)
(*  bound by the Creative Commons BY-SA Attribution-ShareAlike license.    *)
(*                                                                         *)
(*        For additional information concerning CDF licensing, see:        *)
(*                                                                         *)
(*         www.wolfram.com/cdf/adopting-cdf/licensing-options.html         *)
(*                                                                         *)
(*                                                                         *)
(***************************************************************************)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[      1088,         20]
NotebookDataLength[     38456,        866]
NotebookOptionsPosition[     38961,        863]
NotebookOutlinePosition[     39299,        878]
CellTagsIndexPosition[     39256,        875]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell[BoxData[
 RowBox[{"Manipulate", "[", "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"ContourPlot3D", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       RowBox[{
        RowBox[{
         RowBox[{
          SuperscriptBox["10", "c"], " ", 
          RowBox[{"(", 
           RowBox[{
            RowBox[{
             RowBox[{"-", 
              SuperscriptBox["10", "x"]}], " ", "KGR76"}], "+", 
            FractionBox[
             RowBox[{"a1R", " ", "aL", " ", "capacity", " ", "KGR76"}], 
             RowBox[{"dL", "+", "growth"}]]}], ")"}]}], "+", 
         RowBox[{
          RowBox[{"(", 
           RowBox[{
            RowBox[{"-", 
             SuperscriptBox["10", "x"]}], "+", 
            FractionBox[
             RowBox[{"a76", " ", "aL", " ", "capacity"}], 
             RowBox[{"dL", "+", "growth"}]]}], ")"}], " ", 
          SuperscriptBox[
           RowBox[{"(", 
            RowBox[{"1", "+", 
             FractionBox[
              RowBox[{
               SuperscriptBox["a1S", "4"], " ", 
               SuperscriptBox["aT", "4"], " ", 
               SuperscriptBox["capacity", "4"], " ", 
               SuperscriptBox[
                RowBox[{"(", 
                 RowBox[{
                  SuperscriptBox["10", "c\[Prime]"], "+", 
                  FractionBox[
                   RowBox[{
                    SuperscriptBox[
                    RowBox[{"(", 
                    RowBox[{"1", "+", 
                    SuperscriptBox["10", "x"]}], ")"}], "2"], " ", "a81"}], 
                   RowBox[{"a1S", " ", "KGS81"}]]}], ")"}], "4"], " ", 
               SuperscriptBox["KGS81", "4"]}], 
              RowBox[{
               SuperscriptBox[
                RowBox[{"(", 
                 RowBox[{"dT", "+", "growth"}], ")"}], "4"], " ", 
               SuperscriptBox[
                RowBox[{"(", 
                 RowBox[{
                  SuperscriptBox[
                   RowBox[{"(", 
                    RowBox[{"1", "+", 
                    SuperscriptBox["10", "x"]}], ")"}], "2"], "+", 
                  RowBox[{
                   SuperscriptBox["10", "c\[Prime]"], " ", "KGS81"}]}], ")"}],
                 "4"]}]]}], ")"}], "2"]}]}], "\[Equal]", "0"}], ",", 
       RowBox[{
        RowBox[{
         RowBox[{
          RowBox[{"(", 
           RowBox[{
            RowBox[{"-", 
             SuperscriptBox["10", "x"]}], "+", 
            FractionBox[
             RowBox[{"a81", " ", "aT", " ", "capacity"}], 
             RowBox[{"dT", "+", "growth"}]]}], ")"}], " ", 
          SuperscriptBox[
           RowBox[{"(", 
            RowBox[{"1", "+", 
             FractionBox[
              RowBox[{"a1R", " ", "aL", " ", "capacity", " ", 
               RowBox[{"(", 
                RowBox[{
                 SuperscriptBox["10", "c"], "+", 
                 FractionBox[
                  RowBox[{
                   SuperscriptBox[
                    RowBox[{"(", 
                    RowBox[{"1", "+", 
                    SuperscriptBox["10", 
                    RowBox[{"4", " ", "x"}]]}], ")"}], "2"], " ", "a76"}], 
                  RowBox[{"a1R", " ", "KGR76"}]]}], ")"}], " ", "KGR76"}], 
              RowBox[{
               RowBox[{"(", 
                RowBox[{"dL", "+", "growth"}], ")"}], " ", 
               RowBox[{"(", 
                RowBox[{
                 SuperscriptBox[
                  RowBox[{"(", 
                   RowBox[{"1", "+", 
                    SuperscriptBox["10", 
                    RowBox[{"4", " ", "x"}]]}], ")"}], "2"], "+", 
                 RowBox[{
                  SuperscriptBox["10", "c"], " ", "KGR76"}]}], ")"}]}]]}], 
            ")"}], "2"]}], "+", 
         RowBox[{
          SuperscriptBox["10", "c\[Prime]"], " ", 
          RowBox[{"(", 
           RowBox[{
            RowBox[{
             RowBox[{"-", 
              SuperscriptBox["10", "x"]}], " ", "KGS81"}], "+", 
            FractionBox[
             RowBox[{"a1S", " ", "aT", " ", "capacity", " ", "KGS81"}], 
             RowBox[{"dT", "+", "growth"}]]}], ")"}]}]}], "\[Equal]", "0"}]}],
       "}"}], ",", "\[IndentingNewLine]", "\[IndentingNewLine]", 
     RowBox[{"(*", " ", 
      RowBox[{
      "range", " ", "of", " ", "input", " ", "space", " ", "c", " ", 
       "c\[Prime]"}], " ", "*)"}], "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"c", ",", 
       RowBox[{"-", "2"}], ",", 
       FractionBox[
        RowBox[{"Log", "[", 
         FractionBox[
          RowBox[{
           SuperscriptBox["aR33", "2"], " ", 
           SuperscriptBox["capacity", "2"]}], 
          SuperscriptBox[
           RowBox[{"(", 
            RowBox[{"dR", "+", "growth"}], ")"}], "2"]], "]"}], 
        RowBox[{"Log", "[", "10", "]"}]]}], "}"}], ",", "\[IndentingNewLine]",
      "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"c\[Prime]", ",", "2", ",", 
       FractionBox[
        RowBox[{"Log", "[", 
         FractionBox[
          RowBox[{
           SuperscriptBox["aS175", "2"], " ", 
           SuperscriptBox["capacity", "2"]}], 
          SuperscriptBox[
           RowBox[{"(", 
            RowBox[{"dS", "+", "growth"}], ")"}], "2"]], "]"}], 
        RowBox[{"Log", "[", "10", "]"}]]}], "}"}], ",", "\[IndentingNewLine]",
      "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"x", ",", 
       RowBox[{"-", "2"}], ",", 
       RowBox[{
        FractionBox["1", 
         RowBox[{"Log", "[", "10", "]"}]], 
        RowBox[{"Log", "[", 
         RowBox[{"Max", "[", 
          RowBox[{
           FractionBox[
            RowBox[{"a1R", " ", "aL", " ", "capacity", " ", 
             RowBox[{"(", 
              RowBox[{
               FractionBox[
                RowBox[{
                 SuperscriptBox["aR33", "2"], " ", 
                 SuperscriptBox["capacity", "2"]}], 
                SuperscriptBox[
                 RowBox[{"(", 
                  RowBox[{"dR", "+", "growth"}], ")"}], "2"]], "+", 
               FractionBox[
                RowBox[{"a76", " ", 
                 SuperscriptBox[
                  RowBox[{"(", 
                   RowBox[{"1", "+", 
                    FractionBox[
                    RowBox[{
                    SuperscriptBox["a81", "4"], " ", 
                    SuperscriptBox["aT", "4"], " ", 
                    SuperscriptBox["capacity", "4"]}], 
                    SuperscriptBox[
                    RowBox[{"(", 
                    RowBox[{"dT", "+", "growth"}], ")"}], "4"]]}], ")"}], 
                  "2"]}], 
                RowBox[{"a1R", " ", "KGR76"}]]}], ")"}], " ", "KGR76"}], 
            RowBox[{
             RowBox[{"(", 
              RowBox[{"dL", "+", "growth"}], ")"}], " ", 
             RowBox[{"(", 
              RowBox[{
               SuperscriptBox[
                RowBox[{"(", 
                 RowBox[{"1", "+", 
                  FractionBox[
                   RowBox[{
                    SuperscriptBox["a81", "4"], " ", 
                    SuperscriptBox["aT", "4"], " ", 
                    SuperscriptBox["capacity", "4"]}], 
                   SuperscriptBox[
                    RowBox[{"(", 
                    RowBox[{"dT", "+", "growth"}], ")"}], "4"]]}], ")"}], 
                "2"], "+", 
               FractionBox[
                RowBox[{
                 SuperscriptBox["aR33", "2"], " ", 
                 SuperscriptBox["capacity", "2"], " ", "KGR76"}], 
                SuperscriptBox[
                 RowBox[{"(", 
                  RowBox[{"dR", "+", "growth"}], ")"}], "2"]]}], ")"}]}]], 
           ",", 
           FractionBox[
            RowBox[{"a1S", " ", "aT", " ", "capacity", " ", 
             RowBox[{"(", 
              RowBox[{
               FractionBox[
                RowBox[{
                 SuperscriptBox["aS175", "2"], " ", 
                 SuperscriptBox["capacity", "2"]}], 
                SuperscriptBox[
                 RowBox[{"(", 
                  RowBox[{"dS", "+", "growth"}], ")"}], "2"]], "+", 
               FractionBox[
                RowBox[{"a81", " ", 
                 SuperscriptBox[
                  RowBox[{"(", 
                   RowBox[{"1", "+", 
                    FractionBox[
                    RowBox[{"a76", " ", "aL", " ", "capacity"}], 
                    RowBox[{"dL", "+", "growth"}]]}], ")"}], "2"]}], 
                RowBox[{"a1S", " ", "KGS81"}]]}], ")"}], " ", "KGS81"}], 
            RowBox[{
             RowBox[{"(", 
              RowBox[{"dT", "+", "growth"}], ")"}], " ", 
             RowBox[{"(", 
              RowBox[{
               SuperscriptBox[
                RowBox[{"(", 
                 RowBox[{"1", "+", 
                  FractionBox[
                   RowBox[{"a76", " ", "aL", " ", "capacity"}], 
                   RowBox[{"dL", "+", "growth"}]]}], ")"}], "2"], "+", 
               FractionBox[
                RowBox[{
                 SuperscriptBox["aS175", "2"], " ", 
                 SuperscriptBox["capacity", "2"], " ", "KGS81"}], 
                SuperscriptBox[
                 RowBox[{"(", 
                  RowBox[{"dS", "+", "growth"}], ")"}], "2"]]}], ")"}]}]]}], 
          "]"}], "]"}]}]}], "}"}], ",", "\[IndentingNewLine]", 
     "\[IndentingNewLine]", 
     RowBox[{"AxesLabel", "\[Rule]", 
      RowBox[{"{", 
       RowBox[{
       "\"\<\!\(\*SubscriptBox[\(Log\), \(10\)]\)(c)\>\"", ",", 
        "\"\<\!\(\*SubscriptBox[\(Log\), \(10\)]\)(c\[Prime])\>\"", ",", 
        "\"\<\!\(\*SubscriptBox[\(Log\), \(10\)]\)(L), \
\!\(\*SubscriptBox[\(Log\), \(10\)]\)(T)\>\""}], "}"}]}], ",", 
     "\[IndentingNewLine]", "\[IndentingNewLine]", 
     RowBox[{"MeshFunctions", "\[Rule]", 
      RowBox[{"{", "\[IndentingNewLine]", "\[IndentingNewLine]", 
       RowBox[{"(*", " ", 
        RowBox[{"local", " ", "mass", " ", "conserving", " ", "planes"}], " ",
         "*)"}], "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{"Function", "[", 
         RowBox[{
          RowBox[{"{", 
           RowBox[{"c", ",", "c\[Prime]", ",", "x"}], "}"}], ",", 
          "\[IndentingNewLine]", 
          RowBox[{
           RowBox[{"-", 
            FractionBox[
             SuperscriptBox["10", 
              RowBox[{"c", "/", "nR"}]], 
             RowBox[{
              RowBox[{"(", 
               RowBox[{
                SuperscriptBox["10", 
                 RowBox[{"c", "/", "nR"}]], "-", 
                SuperscriptBox[
                 RowBox[{"(", 
                  FractionBox[
                   RowBox[{
                    SuperscriptBox["aR33", "2"], " ", 
                    SuperscriptBox["capacity", "2"]}], 
                   SuperscriptBox[
                    RowBox[{"(", 
                    RowBox[{"dR", "+", "growth"}], ")"}], "2"]], ")"}], 
                 FractionBox["1", "nR"]]}], ")"}], " ", "KR6"}]]}], "-", 
           FractionBox[
            SuperscriptBox["10", 
             RowBox[{"c\[Prime]", "/", "nS"}]], 
            RowBox[{
             RowBox[{"(", 
              RowBox[{
               SuperscriptBox["10", 
                RowBox[{"c\[Prime]", "/", "nS"}]], "-", 
               SuperscriptBox[
                RowBox[{"(", 
                 FractionBox[
                  RowBox[{
                   SuperscriptBox["aS175", "2"], " ", 
                   SuperscriptBox["capacity", "2"]}], 
                  SuperscriptBox[
                   RowBox[{"(", 
                    RowBox[{"dS", "+", "growth"}], ")"}], "2"]], ")"}], 
                FractionBox["1", "nS"]]}], ")"}], " ", "KS12"}]]}]}], "]"}], 
        ",", "\[IndentingNewLine]", 
        RowBox[{"Function", "[", 
         RowBox[{
          RowBox[{"{", 
           RowBox[{"c", ",", "c\[Prime]", ",", "x"}], "}"}], ",", 
          RowBox[{"c", "-", 
           RowBox[{"Log10", "[", " ", 
            RowBox[{
             SuperscriptBox[
              RowBox[{"(", 
               FractionBox[
                RowBox[{"capacity", "*", "aR33"}], 
                RowBox[{"dR", "+", "growth"}]], ")"}], "2"], 
             SuperscriptBox[
              RowBox[{"(", 
               FractionBox[
                SubscriptBox["c", "6"], 
                RowBox[{
                 SubscriptBox["c", "6"], "+", 
                 RowBox[{"1", "/", "KR6"}]}]], ")"}], "nR"]}], "]"}]}]}], 
         "]"}], ",", "\[IndentingNewLine]", 
        RowBox[{"Function", "[", 
         RowBox[{
          RowBox[{"{", 
           RowBox[{"c", ",", "c\[Prime]", ",", "x"}], "}"}], ",", 
          RowBox[{"c\[Prime]", "-", 
           RowBox[{"Log10", "[", 
            RowBox[{
             SuperscriptBox[
              RowBox[{"(", 
               FractionBox[
                RowBox[{"capacity", "*", "aS175"}], 
                RowBox[{"dS", "+", "growth"}]], ")"}], "2"], 
             SuperscriptBox[
              RowBox[{"(", 
               FractionBox[
                SubscriptBox["c", "12"], 
                RowBox[{
                 SubscriptBox["c", "12"], "+", 
                 RowBox[{"1", "/", "KS12"}]}]], ")"}], "nS"]}], "]"}]}]}], 
         "]"}]}], "\[IndentingNewLine]", "}"}]}], ",", "\[IndentingNewLine]", 
     "\[IndentingNewLine]", 
     RowBox[{"(*", " ", 
      RowBox[{"colouring", " ", "options"}], " ", "*)"}], 
     "\[IndentingNewLine]", 
     RowBox[{"MeshStyle", "\[Rule]", 
      RowBox[{"{", "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{"{", 
         RowBox[{
          RowBox[{"Thickness", "[", "0.004", "]"}], ",", 
          RowBox[{"Darker", "[", "Blue", "]"}]}], "}"}], ",", 
        "\[IndentingNewLine]", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{"Thickness", "[", "0.001", "]"}], ",", 
          RowBox[{"Darker", "[", "Gray", "]"}]}], "}"}], ",", 
        "\[IndentingNewLine]", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{"Thickness", "[", "0.001", "]"}], ",", 
          RowBox[{"Darker", "[", "Gray", "]"}]}], "}"}]}], 
       "\[IndentingNewLine]", "}"}]}], ",", "\[IndentingNewLine]", 
     RowBox[{"Mesh", "\[Rule]", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"{", 
         RowBox[{
          SubscriptBox["c", "6"], "+", 
          SubscriptBox["c", "12"]}], "}"}], ",", "\[IndentingNewLine]", 
        RowBox[{"{", "0", "}"}], ",", 
        RowBox[{"{", "0", "}"}]}], "}"}]}], ",", 
     RowBox[{"ImageSize", "\[Rule]", "Full"}], ",", "\[IndentingNewLine]", 
     RowBox[{"ContourStyle", "\[Rule]", 
      RowBox[{"{", "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{"Directive", "[", 
         RowBox[{
          RowBox[{"RGBColor", "[", 
           RowBox[{"0.0", ",", "0.72", ",", "0.92"}], "]"}], ",", 
          RowBox[{"Opacity", "[", "0.5", "]"}]}], "]"}], ",", 
        "\[IndentingNewLine]", 
        RowBox[{"Directive", "[", 
         RowBox[{
          RowBox[{"RGBColor", "[", 
           RowBox[{"1.0", ",", "0.87", ",", "0.0"}], "]"}], ",", 
          RowBox[{"Opacity", "[", "0.5", "]"}]}], "]"}]}], "}"}]}]}], 
    "\[IndentingNewLine]", "]"}], ",", "\[IndentingNewLine]", 
   "\[IndentingNewLine]", 
   RowBox[{"(*", " ", 
    RowBox[{"biophysical", " ", "parameters"}], " ", "*)"}], 
   "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"a76", ",", "1.92"}], "}"}], ",", "0", ",", "10"}], "}"}], ",", 
   "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"a81", ",", "0.822"}], "}"}], ",", "0", ",", "2"}], "}"}], ",", 
   "\[IndentingNewLine]", "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"aL", ",", "0.00117"}], "}"}], ",", "0", ",", "10"}], "}"}], 
   ",", "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"a1R", ",", 
       RowBox[{"1.86", " ", 
        SuperscriptBox["10", "3"]}]}], "}"}], ",", "0", ",", "2000"}], "}"}], 
   ",", "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"KGR76", ",", "0.00145301678863517"}], "}"}], ",", "0", ",", 
     "1"}], "}"}], ",", "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"dL", ",", "0.000117"}], "}"}], ",", "0", ",", "2"}], "}"}], 
   ",", "\[IndentingNewLine]", "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"aT", ",", "0.000951"}], "}"}], ",", "0", ",", "10"}], "}"}], 
   ",", "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"a1S", ",", "704"}], "}"}], ",", "0", ",", "2000"}], "}"}], ",",
    "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"KGS81", ",", 
       RowBox[{"1.00086836995962", " ", 
        SuperscriptBox["10", 
         RowBox[{"-", "5"}]]}]}], "}"}], ",", "0", ",", "1"}], "}"}], ",", 
   "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"dT", ",", "0.666"}], "}"}], ",", "0", ",", "2"}], "}"}], ",", 
   "\[IndentingNewLine]", "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"aR33", ",", "9.1"}], "}"}], ",", "0", ",", "20"}], "}"}], ",", 
   "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"dR", ",", "0.267"}], "}"}], ",", "0", ",", "2"}], "}"}], ",", 
   "\[IndentingNewLine]", "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"aS175", ",", "3.58"}], "}"}], ",", "0", ",", "20"}], "}"}], 
   ",", "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"dS", ",", "0.319"}], "}"}], ",", "0", ",", "2"}], "}"}], ",", 
   "\[IndentingNewLine]", "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"KR6", ",", 
       RowBox[{"3.50695213045775", " ", 
        SuperscriptBox["10", 
         RowBox[{"-", "8"}]]}]}], "}"}], ",", "0", ",", "1"}], "}"}], ",", 
   "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"nR", ",", "0.699"}], "}"}], ",", "0", ",", "2"}], "}"}], ",", 
   "\[IndentingNewLine]", "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"KS12", ",", "0.0095893786931253"}], "}"}], ",", "0", ",", 
     "2"}], "}"}], ",", "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"nS", ",", "1.25"}], "}"}], ",", "0", ",", "2"}], "}"}], ",", 
   "\[IndentingNewLine]", "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"growth", ",", "1"}], "}"}], ",", "0", ",", "3"}], "}"}], ",", 
   "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"capacity", ",", "67.5"}], "}"}], ",", "0", ",", "70"}], "}"}], 
   ",", "\[IndentingNewLine]", "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       SubscriptBox["c", "6"], ",", "120"}], "}"}], ",", "0", ",", "120"}], 
    "}"}], ",", "\[IndentingNewLine]", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       SubscriptBox["c", "12"], ",", "120"}], "}"}], ",", "0", ",", "120"}], 
    "}"}], ",", "\[IndentingNewLine]", 
   RowBox[{"ContentSize", "\[Rule]", 
    RowBox[{"{", 
     RowBox[{"600", ",", "600"}], "}"}]}]}], "\[IndentingNewLine]", 
  "\[IndentingNewLine]", "]"}]], "Input",
 CellChangeTimes->{{3.747389031386036*^9, 3.747389032184869*^9}, {
   3.747389089118073*^9, 3.747389096129272*^9}, {3.747389358615514*^9, 
   3.747389392455998*^9}, {3.747389474974658*^9, 3.747389485177602*^9}, {
   3.747389528426054*^9, 3.747389543806888*^9}, {3.747389579534487*^9, 
   3.747389587093443*^9}, {3.74738963806575*^9, 3.747389703549094*^9}, {
   3.747389733599334*^9, 3.747389808824527*^9}, {3.74738995410583*^9, 
   3.7473900478525*^9}, {3.7473901021766853`*^9, 3.7473901653308773`*^9}, {
   3.747390312166627*^9, 3.74739037517773*^9}, {3.747390420980495*^9, 
   3.747390466868335*^9}, {3.7473905363315783`*^9, 3.7473905822407227`*^9}, {
   3.747390700790504*^9, 3.747390705970769*^9}, {3.747390743665019*^9, 
   3.747390746030817*^9}, {3.7473908568204527`*^9, 3.747390897327035*^9}, {
   3.7473921715429564`*^9, 3.747392187049419*^9}, {3.747392220336143*^9, 
   3.747392295087961*^9}, {3.7473923732117777`*^9, 3.7473923914785843`*^9}, {
   3.747392452402422*^9, 3.747392464897769*^9}, {3.747392613181941*^9, 
   3.747392726639402*^9}, {3.747392792176017*^9, 3.747392831741208*^9}, {
   3.747392962698248*^9, 3.747393099981729*^9}, {3.747393152889573*^9, 
   3.7473932415065613`*^9}, {3.747393803247024*^9, 3.747393811784584*^9}, {
   3.747394009086836*^9, 3.7473940222519407`*^9}, {3.7473940729118834`*^9, 
   3.747394092776291*^9}, {3.7473942177733097`*^9, 3.747394439759956*^9}, {
   3.747394494746484*^9, 3.747394513876671*^9}, 3.747394878960573*^9, {
   3.7473967053684683`*^9, 3.747396710094656*^9}, {3.747396747379023*^9, 
   3.747396777368492*^9}, {3.7473978250936337`*^9, 3.747397844035946*^9}, {
   3.747398247232318*^9, 3.747398313983431*^9}, {3.747398344284812*^9, 
   3.747398400996332*^9}, {3.7473984331100082`*^9, 3.7473984334096203`*^9}, {
   3.747398511845694*^9, 3.7473985127434072`*^9}, 3.747398543582432*^9, {
   3.747398576000943*^9, 3.747398599410327*^9}, {3.747398630863188*^9, 
   3.747398637578126*^9}, {3.747398687371264*^9, 3.747398687723837*^9}, {
   3.74739892109085*^9, 3.747398963217066*^9}, {3.747399023619132*^9, 
   3.747399048800289*^9}, {3.747399116596814*^9, 3.747399238863666*^9}, {
   3.747399271714319*^9, 3.747399272367382*^9}, {3.747399345056735*^9, 
   3.747399351477324*^9}, {3.747399413688953*^9, 3.7473994375974627`*^9}, {
   3.7473995274612303`*^9, 3.747399535520075*^9}, {3.7473998572520514`*^9, 
   3.7473999209532137`*^9}, {3.747400024056685*^9, 3.7474000244795647`*^9}, {
   3.747400223366687*^9, 3.747400278092937*^9}, {3.747401082717723*^9, 
   3.747401164798724*^9}, {3.7474012007470818`*^9, 3.7474012390977287`*^9}, {
   3.747401291524198*^9, 3.747401340086585*^9}, {3.7474014850487947`*^9, 
   3.747401489814784*^9}, 3.747401550983315*^9, {3.747401661171061*^9, 
   3.747401668694418*^9}, {3.747401789406528*^9, 3.74740179342679*^9}, {
   3.747401965149969*^9, 3.747401967094095*^9}, {3.747402156939337*^9, 
   3.747402158777237*^9}, {3.7474025641125813`*^9, 3.7474026042372847`*^9}, {
   3.7474026467352324`*^9, 3.747402768968134*^9}, {3.747402812375635*^9, 
   3.747402814305605*^9}, {3.747403114516549*^9, 3.747403169020748*^9}, 
   3.74740334830191*^9, {3.747725813209064*^9, 3.7477258458540487`*^9}, {
   3.747727544941637*^9, 3.747727586880343*^9}, {3.747728225086191*^9, 
   3.747728231576223*^9}, {3.747728275659957*^9, 3.7477282832962503`*^9}, {
   3.747730183228279*^9, 3.747730185549885*^9}, {3.747730404166831*^9, 
   3.747730404950549*^9}, {3.747730439291437*^9, 3.7477305300090237`*^9}, {
   3.7477306118940887`*^9, 
   3.747730612452181*^9}},ExpressionUUID->"52ed98e4-1cb8-421b-99d5-\
e525bc5730ce"],

Cell[BoxData[
 TagBox[
  StyleBox[
   DynamicModuleBox[{$CellContext`a1R$$ = 1860., $CellContext`a1S$$ = 
    704, $CellContext`a76$$ = 1.92, $CellContext`a81$$ = 
    0.822, $CellContext`aL$$ = 0.00117, $CellContext`aR33$$ = 
    9.1, $CellContext`aS175$$ = 3.58, $CellContext`aT$$ = 
    0.000951, $CellContext`capacity$$ = 67.5, $CellContext`dL$$ = 
    0.000117, $CellContext`dR$$ = 0.267, $CellContext`dS$$ = 
    0.319, $CellContext`dT$$ = 0.666, $CellContext`growth$$ = 
    1, $CellContext`KGR76$$ = 0.00145301678863517, $CellContext`KGS81$$ = 
    0.0000100086836995962, $CellContext`KR6$$ = 
    3.50695213045775*^-8, $CellContext`KS12$$ = 
    0.0095893786931253, $CellContext`nR$$ = 0.699, $CellContext`nS$$ = 
    1.25, $CellContext`$1155$$ = 120, $CellContext`$1156$$ = 78., 
    Typeset`show$$ = True, Typeset`bookmarkList$$ = {}, 
    Typeset`bookmarkMode$$ = "Menu", Typeset`animator$$, Typeset`animvar$$ = 
    1, Typeset`name$$ = "\"untitled\"", Typeset`specs$$ = {{{
       Hold[$CellContext`a76$$], 1.92}, 0, 10}, {{
       Hold[$CellContext`a81$$], 0.822}, 0, 2}, {{
       Hold[$CellContext`aL$$], 0.00117}, 0, 10}, {{
       Hold[$CellContext`a1R$$], 1860.}, 0, 2000}, {{
       Hold[$CellContext`KGR76$$], 0.00145301678863517}, 0, 1}, {{
       Hold[$CellContext`dL$$], 0.000117}, 0, 2}, {{
       Hold[$CellContext`aT$$], 0.000951}, 0, 10}, {{
       Hold[$CellContext`a1S$$], 704}, 0, 2000}, {{
       Hold[$CellContext`KGS81$$], 0.0000100086836995962}, 0, 1}, {{
       Hold[$CellContext`dT$$], 0.666}, 0, 2}, {{
       Hold[$CellContext`aR33$$], 9.1}, 0, 20}, {{
       Hold[$CellContext`dR$$], 0.267}, 0, 2}, {{
       Hold[$CellContext`aS175$$], 3.58}, 0, 20}, {{
       Hold[$CellContext`dS$$], 0.319}, 0, 2}, {{
       Hold[$CellContext`KR6$$], 3.50695213045775*^-8}, 0, 1}, {{
       Hold[$CellContext`nR$$], 0.699}, 0, 2}, {{
       Hold[$CellContext`KS12$$], 0.0095893786931253}, 0, 2}, {{
       Hold[$CellContext`nS$$], 1.25}, 0, 2}, {{
       Hold[$CellContext`growth$$], 1}, 0, 3}, {{
       Hold[$CellContext`capacity$$], 67.5}, 0, 70}, {{
       Hold[$CellContext`$1155$$], 120, 
       RawBoxes[
        SubscriptBox["c", "6"]]}, 0, 120}, {{
       Hold[$CellContext`$1156$$], 120, 
       RawBoxes[
        SubscriptBox["c", "12"]]}, 0, 120}}, Typeset`size$$ = {
    580., {244., 248.}}, Typeset`update$$ = 0, Typeset`initDone$$, 
    Typeset`skipInitDone$$ = True, $CellContext`a76$19038$$ = 
    0, $CellContext`a81$19039$$ = 0, $CellContext`aL$19040$$ = 
    0, $CellContext`a1R$19041$$ = 0, $CellContext`KGR76$19042$$ = 
    0, $CellContext`dL$19043$$ = 0, $CellContext`aT$19044$$ = 
    0, $CellContext`a1S$19045$$ = 0, $CellContext`KGS81$19046$$ = 
    0, $CellContext`dT$19047$$ = 0}, 
    DynamicBox[Manipulate`ManipulateBoxes[
     1, StandardForm, 
      "Variables" :> {$CellContext`a1R$$ = 1860., $CellContext`a1S$$ = 
        704, $CellContext`a76$$ = 1.92, $CellContext`a81$$ = 
        0.822, $CellContext`aL$$ = 0.00117, $CellContext`aR33$$ = 
        9.1, $CellContext`aS175$$ = 3.58, $CellContext`aT$$ = 
        0.000951, $CellContext`capacity$$ = 67.5, $CellContext`dL$$ = 
        0.000117, $CellContext`dR$$ = 0.267, $CellContext`dS$$ = 
        0.319, $CellContext`dT$$ = 0.666, $CellContext`growth$$ = 
        1, $CellContext`KGR76$$ = 0.00145301678863517, $CellContext`KGS81$$ = 
        0.0000100086836995962, $CellContext`KR6$$ = 
        3.50695213045775*^-8, $CellContext`KS12$$ = 
        0.0095893786931253, $CellContext`nR$$ = 0.699, $CellContext`nS$$ = 
        1.25, $CellContext`$1155$$ = 120, $CellContext`$1156$$ = 120}, 
      "ControllerVariables" :> {
        Hold[$CellContext`a76$$, $CellContext`a76$19038$$, 0], 
        Hold[$CellContext`a81$$, $CellContext`a81$19039$$, 0], 
        Hold[$CellContext`aL$$, $CellContext`aL$19040$$, 0], 
        Hold[$CellContext`a1R$$, $CellContext`a1R$19041$$, 0], 
        Hold[$CellContext`KGR76$$, $CellContext`KGR76$19042$$, 0], 
        Hold[$CellContext`dL$$, $CellContext`dL$19043$$, 0], 
        Hold[$CellContext`aT$$, $CellContext`aT$19044$$, 0], 
        Hold[$CellContext`a1S$$, $CellContext`a1S$19045$$, 0], 
        Hold[$CellContext`KGS81$$, $CellContext`KGS81$19046$$, 0], 
        Hold[$CellContext`dT$$, $CellContext`dT$19047$$, 0]}, 
      "OtherVariables" :> {
       Typeset`show$$, Typeset`bookmarkList$$, Typeset`bookmarkMode$$, 
        Typeset`animator$$, Typeset`animvar$$, Typeset`name$$, 
        Typeset`specs$$, Typeset`size$$, Typeset`update$$, Typeset`initDone$$,
         Typeset`skipInitDone$$}, "Body" :> 
      ContourPlot3D[{
        10^$CellContext`c ((-10^$CellContext`x) $CellContext`KGR76$$ + \
$CellContext`a1R$$ $CellContext`aL$$ $CellContext`capacity$$ \
$CellContext`KGR76$$/($CellContext`dL$$ + $CellContext`growth$$)) + \
(-10^$CellContext`x + $CellContext`a76$$ $CellContext`aL$$ \
$CellContext`capacity$$/($CellContext`dL$$ + $CellContext`growth$$)) (
             1 + $CellContext`a1S$$^4 $CellContext`aT$$^4 \
$CellContext`capacity$$^4 (
                10^$CellContext`c\[Prime] + (1 + 
                   10^$CellContext`x)^2 \
$CellContext`a81$$/($CellContext`a1S$$ $CellContext`KGS81$$))^4 \
$CellContext`KGS81$$^4/(($CellContext`dT$$ + $CellContext`growth$$)^4 ((1 + 
                  10^$CellContext`x)^2 + 
                10^$CellContext`c\[Prime] $CellContext`KGS81$$)^4))^2 == 
         0, (-10^$CellContext`x + $CellContext`a81$$ $CellContext`aT$$ \
$CellContext`capacity$$/($CellContext`dT$$ + $CellContext`growth$$)) (
             1 + $CellContext`a1R$$ $CellContext`aL$$ $CellContext`capacity$$ \
(10^$CellContext`c + (1 + 
                  10^(4 $CellContext`x))^2 \
$CellContext`a76$$/($CellContext`a1R$$ $CellContext`KGR76$$)) \
$CellContext`KGR76$$/(($CellContext`dL$$ + $CellContext`growth$$) ((1 + 
                 10^(4 $CellContext`x))^2 + 
               10^$CellContext`c $CellContext`KGR76$$)))^2 + 
          10^$CellContext`c\[Prime] ((-10^$CellContext`x) \
$CellContext`KGS81$$ + $CellContext`a1S$$ $CellContext`aT$$ \
$CellContext`capacity$$ $CellContext`KGS81$$/($CellContext`dT$$ + \
$CellContext`growth$$)) == 0}, {$CellContext`c, -2, 
         Log[$CellContext`aR33$$^2 \
$CellContext`capacity$$^2/($CellContext`dR$$ + $CellContext`growth$$)^2]/Log[
         10]}, {$CellContext`c\[Prime], 2, 
         Log[$CellContext`aS175$$^2 \
$CellContext`capacity$$^2/($CellContext`dS$$ + $CellContext`growth$$)^2]/Log[
         10]}, {$CellContext`x, -2, (1/Log[10]) Log[
           
           Max[$CellContext`a1R$$ $CellContext`aL$$ $CellContext`capacity$$ \
($CellContext`aR33$$^2 $CellContext`capacity$$^2/($CellContext`dR$$ + \
$CellContext`growth$$)^2 + $CellContext`a76$$ (
                1 + $CellContext`a81$$^4 $CellContext`aT$$^4 \
$CellContext`capacity$$^4/($CellContext`dT$$ + \
$CellContext`growth$$)^4)^2/($CellContext`a1R$$ $CellContext`KGR76$$)) \
$CellContext`KGR76$$/(($CellContext`dL$$ + $CellContext`growth$$) ((
               1 + $CellContext`a81$$^4 $CellContext`aT$$^4 \
$CellContext`capacity$$^4/($CellContext`dT$$ + $CellContext`growth$$)^4)^2 + \
$CellContext`aR33$$^2 $CellContext`capacity$$^2 \
$CellContext`KGR76$$/($CellContext`dR$$ + $CellContext`growth$$)^2)), \
$CellContext`a1S$$ $CellContext`aT$$ $CellContext`capacity$$ \
($CellContext`aS175$$^2 $CellContext`capacity$$^2/($CellContext`dS$$ + \
$CellContext`growth$$)^2 + $CellContext`a81$$ (
                1 + $CellContext`a76$$ $CellContext`aL$$ \
$CellContext`capacity$$/($CellContext`dL$$ + \
$CellContext`growth$$))^2/($CellContext`a1S$$ $CellContext`KGS81$$)) \
$CellContext`KGS81$$/(($CellContext`dT$$ + $CellContext`growth$$) ((
               1 + $CellContext`a76$$ $CellContext`aL$$ \
$CellContext`capacity$$/($CellContext`dL$$ + $CellContext`growth$$))^2 + \
$CellContext`aS175$$^2 $CellContext`capacity$$^2 \
$CellContext`KGS81$$/($CellContext`dS$$ + $CellContext`growth$$)^2))]]}, 
        AxesLabel -> {
         "\!\(\*SubscriptBox[\(Log\), \(10\)]\)(c)", 
          "\!\(\*SubscriptBox[\(Log\), \(10\)]\)(c\[Prime])", 
          "\!\(\*SubscriptBox[\(Log\), \(10\)]\)(L), \!\(\*SubscriptBox[\(Log\
\), \(10\)]\)(T)"}, MeshFunctions -> {
          
          Function[{$CellContext`c$, $CellContext`c\[Prime]$, \
$CellContext`x$}, -(
             10^($CellContext`c$/$CellContext`nR$$)/((
              10^($CellContext`c$/$CellContext`nR$$) - ($CellContext`aR33$$^2 \
$CellContext`capacity$$^2/($CellContext`dR$$ + $CellContext`growth$$)^2)^(
               1/$CellContext`nR$$)) $CellContext`KR6$$)) - 
           10^($CellContext`c\[Prime]$/$CellContext`nS$$)/((
            10^($CellContext`c\[Prime]$/$CellContext`nS$$) - \
($CellContext`aS175$$^2 $CellContext`capacity$$^2/($CellContext`dS$$ + \
$CellContext`growth$$)^2)^(1/$CellContext`nS$$)) $CellContext`KS12$$)], 
          
          Function[{$CellContext`c$, $CellContext`c\[Prime]$, \
$CellContext`x$}, $CellContext`c$ - 
           Log10[($CellContext`capacity$$ \
$CellContext`aR33$$/($CellContext`dR$$ + $CellContext`growth$$))^2 \
($CellContext`$1155$$/($CellContext`$1155$$ + 
              1/$CellContext`KR6$$))^$CellContext`nR$$]], 
          
          Function[{$CellContext`c$, $CellContext`c\[Prime]$, \
$CellContext`x$}, $CellContext`c\[Prime]$ - 
           Log10[($CellContext`capacity$$ \
$CellContext`aS175$$/($CellContext`dS$$ + $CellContext`growth$$))^2 \
($CellContext`$1156$$/($CellContext`$1156$$ + 
              1/$CellContext`KS12$$))^$CellContext`nS$$]]}, MeshStyle -> {{
           Thickness[0.004], 
           Darker[Blue]}, {
           Thickness[0.001], 
           Darker[Gray]}, {
           Thickness[0.001], 
           Darker[Gray]}}, 
        Mesh -> {{$CellContext`$1155$$ + $CellContext`$1156$$}, {0}, {0}}, 
        ImageSize -> Full, ContourStyle -> {
          Directive[
           RGBColor[0., 0.72, 0.92], 
           Opacity[0.5]], 
          Directive[
           RGBColor[1., 0.87, 0.], 
           Opacity[0.5]]}], 
      "Specifications" :> {{{$CellContext`a76$$, 1.92}, 0, 
         10}, {{$CellContext`a81$$, 0.822}, 0, 
         2}, {{$CellContext`aL$$, 0.00117}, 0, 
         10}, {{$CellContext`a1R$$, 1860.}, 0, 
         2000}, {{$CellContext`KGR76$$, 0.00145301678863517}, 0, 
         1}, {{$CellContext`dL$$, 0.000117}, 0, 
         2}, {{$CellContext`aT$$, 0.000951}, 0, 
         10}, {{$CellContext`a1S$$, 704}, 0, 
         2000}, {{$CellContext`KGS81$$, 0.0000100086836995962}, 0, 
         1}, {{$CellContext`dT$$, 0.666}, 0, 2}, {{$CellContext`aR33$$, 9.1}, 
         0, 20}, {{$CellContext`dR$$, 0.267}, 0, 
         2}, {{$CellContext`aS175$$, 3.58}, 0, 
         20}, {{$CellContext`dS$$, 0.319}, 0, 
         2}, {{$CellContext`KR6$$, 3.50695213045775*^-8}, 0, 
         1}, {{$CellContext`nR$$, 0.699}, 0, 
         2}, {{$CellContext`KS12$$, 0.0095893786931253}, 0, 
         2}, {{$CellContext`nS$$, 1.25}, 0, 2}, {{$CellContext`growth$$, 1}, 
         0, 3}, {{$CellContext`capacity$$, 67.5}, 0, 
         70}, {{$CellContext`$1155$$, 120, 
          RawBoxes[
           SubscriptBox["c", "6"]]}, 0, 120}, {{$CellContext`$1156$$, 120, 
          RawBoxes[
           SubscriptBox["c", "12"]]}, 0, 120}}, 
      "Options" :> {ContentSize -> {600, 600}}, "DefaultOptions" :> {}],
     ImageSizeCache->{882., {315., 320.}},
     SingleEvaluation->True],
    Deinitialization:>None,
    DynamicModuleValues:>{},
    SynchronousInitialization->True,
    UndoTrackedVariables:>{Typeset`show$$, Typeset`bookmarkMode$$},
    UnsavedVariables:>{Typeset`initDone$$},
    UntrackedVariables:>{Typeset`size$$}], "Manipulate",
   Deployed->True,
   StripOnInput->False],
  Manipulate`InterpretManipulate[1]]], "Output",
 GeneratedCell->False,
 CellAutoOverwrite->False,
 CellChangeTimes->{
  3.747389033748398*^9, 3.747389103232154*^9, 3.74738939581768*^9, 
   3.747389587791507*^9, {3.747389771483369*^9, 3.747389810322426*^9}, {
   3.7473899962897587`*^9, 3.74739005048315*^9}, {3.747390106389085*^9, 
   3.74739011944261*^9}, {3.74739014961051*^9, 3.7473901664628*^9}, 
   3.747390376207191*^9, 3.7473905158512497`*^9, {3.7473905634053307`*^9, 
   3.7473905745708447`*^9}, 3.747390706868517*^9, 3.747390747086545*^9, 
   3.747390910963441*^9, 3.74739178773517*^9, 3.747391972830592*^9, 
   3.747392188770481*^9, {3.747392225696656*^9, 3.747392239468734*^9}, {
   3.747392283528406*^9, 3.747392296005834*^9}, {3.7473923742543373`*^9, 
   3.747392392782117*^9}, {3.747392454493553*^9, 3.7473924656868467`*^9}, {
   3.747392715056752*^9, 3.7473927274942713`*^9}, {3.747392796123869*^9, 
   3.7473928325544033`*^9}, {3.747392977152059*^9, 3.7473931007143173`*^9}, {
   3.747393164130526*^9, 3.747393242268786*^9}, 3.747393275322805*^9, 
   3.7473940247559977`*^9, {3.747394241028718*^9, 3.7473944405544033`*^9}, {
   3.747394503094947*^9, 3.7473945153702707`*^9}, 3.7473946193185663`*^9, 
   3.747394880414385*^9, 3.747395193313669*^9, 3.747395352861328*^9, 
   3.7473967121870403`*^9, 3.74739736203403*^9, {3.7473978182397118`*^9, 
   3.747397845410696*^9}, {3.747398249782784*^9, 3.7473982972199507`*^9}, {
   3.747398385016014*^9, 3.747398402089199*^9}, 3.747398434927126*^9, 
   3.747398513752076*^9, 3.747398544834832*^9, 3.747398600570704*^9, 
   3.747398639829867*^9, 3.747398688828745*^9, {3.7473991861670103`*^9, 
   3.747399240072234*^9}, 3.747399274046216*^9, 3.7473993735285177`*^9, 
   3.747399441290389*^9, {3.7473995312540827`*^9, 3.7473995608647757`*^9}, {
   3.747399878435495*^9, 3.747399921944111*^9}, {3.747400232071891*^9, 
   3.747400279111176*^9}, 3.747401165892241*^9, {3.7474012182971163`*^9, 
   3.747401240740753*^9}, {3.747401296991517*^9, 3.747401343112605*^9}, 
   3.747401491428276*^9, 3.747401569164694*^9, 3.747401670110915*^9, 
   3.7474017944431334`*^9, {3.747401841759651*^9, 3.747401869298086*^9}, 
   3.74740194031612*^9, 3.747401971133642*^9, {3.747402150684732*^9, 
   3.7474021608536043`*^9}, 3.747402202090272*^9, 3.747402428375476*^9, 
   3.7474025763002987`*^9, 3.7474026639754553`*^9, {3.747402710643634*^9, 
   3.747402770784511*^9}, 3.747402815217278*^9, {3.747403101829006*^9, 
   3.7474031715417747`*^9}, 3.7474032233717403`*^9, 3.747403354054719*^9, 
   3.7474034257423973`*^9, 3.747403624948254*^9, 3.747725697034576*^9, {
   3.7477258352126207`*^9, 3.74772584699594*^9}, 3.747727219369602*^9, 
   3.747727590318657*^9, 3.747727641429296*^9, 3.747728237006011*^9, 
   3.74772828405066*^9, 3.747730190524681*^9, 3.747730355845407*^9, 
   3.7477304062015467`*^9, 3.7477304875773907`*^9, {3.747730522524756*^9, 
   3.747730535159205*^9}, 
   3.7477306164354067`*^9},ExpressionUUID->"f667eb85-8f03-4eec-a575-\
fda6e39e39c8"]
}, {2}]]
},
WindowSize->{960, 1030},
WindowMargins->{{0, Automatic}, {0, Automatic}},
FrontEndVersion->"11.2 for Linux x86 (64-bit) (September 10, 2017)",
StyleDefinitions->"Default.nb"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[CellGroupData[{
Cell[1510, 35, 22826, 564, 2374, "Input",ExpressionUUID->"52ed98e4-1cb8-421b-99d5-e525bc5730ce"],
Cell[24339, 601, 14609, 259, 654, "Output",ExpressionUUID->"f667eb85-8f03-4eec-a575-fda6e39e39c8"]
}, {2}]]
}
]
*)

(* NotebookSignature evpAI6WwEc#0uAwsNxvI#zCJ *)
