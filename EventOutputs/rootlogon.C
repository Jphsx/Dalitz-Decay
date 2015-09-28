{

//ild TStyle
TStyle* ildStyle = new  TStyle("ildStyle", "ILD Style (from ~/rootlogon.C)");

//set the background color to white
ildStyle->SetHistMinimumZero(true);               // Suppresses Histogram Zero Supression
ildStyle->SetFillColor(10);
ildStyle->SetFrameFillColor(10);
ildStyle->SetCanvasColor(10);
ildStyle->SetPadColor(10);
ildStyle->SetTitleFillColor(0);
ildStyle->SetStatColor(10);

//dont put a colored frame around the plots
ildStyle->SetFrameBorderMode(0);
ildStyle->SetCanvasBorderMode(0);
ildStyle->SetPadBorderMode(0);
ildStyle->SetLegendBorderSize(0);

//use the primary color palette
ildStyle->SetPalette(1,0);

//set the default line color for a histogram to be black
ildStyle->SetHistLineColor(kBlack);

//set the default line color for a fit function to be red
ildStyle->SetFuncColor(kRed);

//make the axis labels black
ildStyle->SetLabelColor(kBlack,"xyz");

//set the default title color to be white (no title) 
ildStyle->SetTitleColor(kBlack);
 
//set the margins
ildStyle->SetPadBottomMargin(0.12);
ildStyle->SetPadTopMargin(0.08);
ildStyle->SetPadRightMargin(0.08);
ildStyle->SetPadLeftMargin(0.17);

//set axis label and title text sizes
ildStyle->SetLabelFont(42,"xyz");
ildStyle->SetLabelSize(0.04,"xyz");
ildStyle->SetLabelOffset(0.005,"xyz");
ildStyle->SetTitleFont(42,"xyz");
ildStyle->SetTitleSize(0.05,"xyz");
//ildStyle->SetTitleOffset(1.6,"yz");  // This is important for the y-axis labels being resolved from the y-axis title
ildStyle->SetTitleOffset(1.1,"yz");  // This is important for the y-axis labels being resolved from the y-axis title
ildStyle->SetTitleOffset(0.8,"x");
ildStyle->SetStatFont(42);
ildStyle->SetStatFontSize(0.05);
ildStyle->SetTitleBorderSize(0);
ildStyle->SetStatBorderSize(2);
ildStyle->SetTextFont(42);
ildStyle->SetStatW(0.25);
ildStyle->SetStatH(0.20);

//GWW ??
//ildStyle->SetTitleX(0.5);
//ildStyle->SetTitleY(1.0);

//set line widths
ildStyle->SetFrameLineWidth(2);
ildStyle->SetFuncWidth(2);
ildStyle->SetHistLineWidth(2);

//set the number of divisions to show
ildStyle->SetNdivisions(506, "xy");

//turn off xy grids
ildStyle->SetPadGridX(0);
ildStyle->SetPadGridY(0);

//set the tick mark style
ildStyle->SetPadTickX(1);
ildStyle->SetPadTickY(1);

//turn off stats
ildStyle->SetOptStat("nemrou");
ildStyle->SetOptFit(0);

//marker settings
ildStyle->SetMarkerColor(kBlack);
ildStyle->SetMarkerStyle(20);
ildStyle->SetMarkerSize(0.8);
ildStyle->SetLineWidth(0); 

//ild TStyle
TStyle* ildStyle2 = new  TStyle("ildStyle2", "ILD Style 2 (from ~/rootlogon.C)");

//set the background color to white
ildStyle2->SetFillColor(10);
ildStyle2->SetFrameFillColor(10);
ildStyle2->SetCanvasColor(10);
ildStyle2->SetPadColor(10);
ildStyle2->SetTitleFillColor(0);
ildStyle2->SetStatColor(10);

//dont put a colored frame around the plots
ildStyle2->SetFrameBorderMode(0);
ildStyle2->SetCanvasBorderMode(0);
ildStyle2->SetPadBorderMode(0);
ildStyle2->SetLegendBorderSize(0);

//use the primary color palette
ildStyle2->SetPalette(1,0);

//set the default line color for a histogram to be black
ildStyle2->SetHistLineColor(kBlue);

//set the default line color for a fit function to be red
ildStyle2->SetFuncColor(kRed);

//make the axis labels black
ildStyle2->SetLabelColor(kBlack,"xyz");

//set the default title color to be white (no title) 
ildStyle2->SetTitleColor(kBlack);
 
//set the margins
ildStyle2->SetPadBottomMargin(0.18);
ildStyle2->SetPadTopMargin(0.08);
ildStyle2->SetPadRightMargin(0.08);
ildStyle2->SetPadLeftMargin(0.17);

//set axis label and title text sizes
ildStyle2->SetLabelFont(42,"xyz");
ildStyle2->SetLabelSize(0.04,"xyz");
ildStyle2->SetLabelOffset(0.015,"xyz");
ildStyle2->SetTitleFont(42,"xyz");
ildStyle2->SetTitleSize(0.04,"xyz");
ildStyle2->SetTitleOffset(1.1,"yz");
ildStyle2->SetTitleOffset(1.0,"x");
ildStyle2->SetStatFont(42);
ildStyle2->SetStatFontSize(0.05);
ildStyle2->SetTitleBorderSize(0);
ildStyle2->SetStatBorderSize(2);
ildStyle2->SetTextFont(42);

//set line widths
ildStyle2->SetFrameLineWidth(2);
ildStyle2->SetFuncWidth(2);
ildStyle2->SetHistLineWidth(2);

//set the number of divisions to show
ildStyle2->SetNdivisions(506, "xy");

//turn off xy grids
ildStyle2->SetPadGridX(0);
ildStyle2->SetPadGridY(0);

//set the tick mark style
ildStyle2->SetPadTickX(1);
ildStyle2->SetPadTickY(1);

//turn off stats
//ildStyle2->SetOptStat(0);
//ildStyle2->SetOptFit(0);

//marker settings
ildStyle2->SetMarkerColor(kBlue);
ildStyle2->SetMarkerStyle(22);
ildStyle2->SetMarkerSize(0.7);
ildStyle2->SetLineWidth(0); 

//ild TStyle
TStyle* ildStyle3 = new  TStyle("ildStyle3", "ILD Style 3 (from ~/rootlogon.C)");

//set the background color to white
ildStyle3->SetFillColor(10);
ildStyle3->SetFrameFillColor(10);
ildStyle3->SetCanvasColor(10);
ildStyle3->SetPadColor(10);
ildStyle3->SetTitleFillColor(0);
ildStyle3->SetStatColor(10);

//dont put a colored frame around the plots
ildStyle3->SetFrameBorderMode(0);
ildStyle3->SetCanvasBorderMode(0);
ildStyle3->SetPadBorderMode(0);
ildStyle3->SetLegendBorderSize(0);

//use the primary color palette
ildStyle3->SetPalette(1,0);

//set the default line color for a histogram to be black
ildStyle3->SetHistLineColor(kBlack);

//set the default line color for a fit function to be red
ildStyle3->SetFuncColor(kRed);

//make the axis labels black
ildStyle3->SetLabelColor(kBlack,"xyz");

//set the default title color to be white (no title) 
ildStyle3->SetTitleColor(kBlack);
 
//set the margins
ildStyle3->SetPadBottomMargin(0.18);
ildStyle3->SetPadTopMargin(0.08);
ildStyle3->SetPadRightMargin(0.08);
ildStyle3->SetPadLeftMargin(0.17);

//set axis label and title text sizes
ildStyle3->SetLabelFont(42,"xyz");
ildStyle3->SetLabelSize(0.04,"xyz");
ildStyle3->SetLabelOffset(0.015,"xyz");
ildStyle3->SetTitleFont(42,"xyz");
ildStyle3->SetTitleSize(0.04,"xyz");
ildStyle3->SetTitleOffset(1.1,"yz");
ildStyle3->SetTitleOffset(1.0,"x");
ildStyle3->SetStatFont(42);
ildStyle3->SetStatFontSize(0.02);
ildStyle3->SetTitleBorderSize(0);
ildStyle3->SetStatBorderSize(0);
ildStyle3->SetTextFont(42);

//set line widths
ildStyle3->SetFrameLineWidth(2);
ildStyle3->SetFuncWidth(2);
ildStyle3->SetHistLineWidth(2);

//set the number of divisions to show
ildStyle3->SetNdivisions(506, "xy");

//turn off xy grids
ildStyle3->SetPadGridX(0);
ildStyle3->SetPadGridY(0);

//set the tick mark style
ildStyle3->SetPadTickX(1);
ildStyle3->SetPadTickY(1);

//turn off stats
ildStyle3->SetOptStat(0);
ildStyle3->SetOptFit(0);

//marker settings
ildStyle3->SetMarkerColor(kBlue);
ildStyle3->SetMarkerStyle(20);
ildStyle3->SetMarkerSize(0.4);
ildStyle3->SetLineWidth(0); 


//done
ildStyle->cd();
//ildStyle->cd();
gROOT->ForceStyle();
gStyle->ls();

}

