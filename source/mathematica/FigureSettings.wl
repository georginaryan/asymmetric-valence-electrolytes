(* ::Package:: *)

BeginPackage["FigureSettings`"];

InitializePlotSettings::usage = "InitializePlotSettings[] applies global SetOptions to Plot and defines the custom colors and linestyles variables.";
colors::usage = "A list of blue tones.";
DivergentColoursv2::usage = "A list of tones in a divergent colour scheme from purple -> grey -> blue.";
NumCompColours::usage= "Colours for Figure 2.";
dotting::usage = "A list of dashing patterns from solid to dots.";
dottingv2::usage = "A list of dashing patterns from dots to solid and back to dots.";
linestyles::usage = "A list containing Thick, Dashed, and Dotted.";
Thick1::usage = "A list containing 5 line thickness of 0.008.";
Thick2::usage = "A list containing 5 line thickness of 0.01.";


Begin["`Private`"];

InitializePlotSettings[] := (
  SetOptions[Plot,
    BaseStyle -> {FontSize -> 28, FontFamily -> "Latin Modern Roman"},
    PlotStyle -> Thickness[0.006],
    AxesStyle -> False,
    GridLines -> Automatic,
    GridLinesStyle->Directive[GrayLevel[0.7],Thick,Thickness[0.001]],
    Frame -> True,
    FrameStyle -> Directive[Black, Thickness[0.004]],
    FrameTicksStyle -> Directive[FontSize -> 28, FontFamily -> "Latin Modern Roman"],
    PlotRangeClipping -> False, 
    Background -> None
  ];
  
  SetOptions[ListLinePlot,
	BaseStyle->{FontSize->36,FontFamily->"Latin Modern Roman"},
	PlotStyle->Thickness[0.0060],
	Axes->False,
	GridLines->Automatic,
	GridLinesStyle->Directive[GrayLevel[0.7],Dashed],
	Frame->True,
	FrameStyle->Directive[Black,Thickness[0.004]],
	FrameTicksStyle->Directive[FontSize->36,FontFamily->"Latin Modern Roman"],
	PlotRangeClipping->True,
	ImageSize->500,Background->None];
  
  colors = Reverse[{
    RGBColor[0.00, 0.20, 0.40], (*very dark blue*)
    RGBColor[0.00, 0.32, 0.60], 
    RGBColor[0.18, 0.46, 0.70], 
    RGBColor[0.40, 0.62, 0.80], 
    RGBColor[0.62, 0.76, 0.88] (*light blue*)
  }];
  
  t = 3;n=0.008;
  dotting = {
    AbsoluteThickness[t], 
    Directive[AbsoluteThickness[t], Dashing[{0.03, 0.02}]], 
    Directive[AbsoluteThickness[t], Dashing[{0.02, 0.02}]], 
    Directive[AbsoluteThickness[t], Dashing[{0.01, 0.02}]], 
    Directive[AbsoluteThickness[t], Dashing[{0.005, 0.02}]]
  };

  linestyles = {Thick, Dashed, Dotted};


	DivergentColoursv2={
	RGBColor["#2C7BB6"],
	RGBColor["#ABD9E9"],
	RGBColor["#BDBDBD"],
	RGBColor["#AF8DC3"],
	RGBColor["#762A83"]   (*dark purple*)};
	
	NumCompColours={RGBColor["#90C3DD"],
	RGBColor["#D7B5D8"],
	RGBColor["#1F5A8A"],
	RGBColor["#762A83"]};
	
	dottingv2={
	Directive[AbsoluteThickness[t],Dashing[{0.01,0.02}]],
	Directive[AbsoluteThickness[t],Dashing[{0.04,0.02}]],
	AbsoluteThickness[t],
	Directive[AbsoluteThickness[t],Dashing[{0.04,0.02}]],
	Directive[AbsoluteThickness[t],Dashing[{0.01,0.02}]]};
	
	
	Thick1={Thickness[n],Thickness[n],Thickness[n],Thickness[n],Thickness[n]};
	Thick2={Thickness[0.01],Thickness[0.01],Thickness[0.01],Thickness[0.01],Thickness[0.01]};
);
End[];
EndPackage[];
