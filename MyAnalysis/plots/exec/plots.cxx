#include <iostream>
#include <string>
#include <sstream>
#include <array>

#include "Process.h"

#include "TFile.h"
#include "TTree.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TText.h"
#include "TLatex.h"
#include "THStack.h"
#include "TLine.h"
#include "TEllipse.h"


struct processInfos
{
	std::string displayName = "";
	Color_t color = kBlack;
};

/*
static const std::map<std::string , Color_t> processMap =
{
	{"1, signal" , kBlack } ,
	{"2, Other Higgs decays" , kCyan+2 } , 
	{"3, Higgsstrahlung other Z decays" , kCyan} ,
	{"5, 2 fermions" , kRed-3 } ,
	{"4, 4 fermions" , kRed-7 } 
} ;
*/

static const std::map<int , processInfos> processMap =
{
	{0,{"signal" , kBlack }} ,
	{1,{"Other Higgs decays" , kCyan+2 }} , 
	{2,{"Higgsstrahlung other Z decays" , kCyan}} ,
	{3,{"2 fermions" , kRed-3 }} ,
	{4,{"4 fermions" , kRed-7 }}
} ;

static const std::vector<std::string> processNames = {"signal","Other Higgs decays","Higgsstrahlung other Z decays","2 fermions","4 fermions"};



void setStyle(TH1* histo)
{
	histo->SetLineWidth(2) ;
	histo->GetXaxis()->SetLabelSize(0.025f) ;
	histo->GetYaxis()->SetLabelSize(0.025f) ;

	histo->GetYaxis()->SetTitleOffset(1.5f) ;

	histo->GetXaxis()->SetTitleFont(62) ;
	histo->GetYaxis()->SetTitleFont(62) ;
}
void setStyle(THStack* histo)
{
	histo->GetXaxis()->SetLabelSize(0.025f) ;
	histo->GetYaxis()->SetLabelSize(0.025f) ;

	histo->GetYaxis()->SetTitleOffset(1.55f) ;

	histo->GetXaxis()->SetTitleFont(62) ;
	histo->GetYaxis()->SetTitleFont(62) ;
}

void setStyle(TCanvas* canvas)
{
	canvas->SetTicks() ;
	canvas->SetTopMargin(0.05f) ;
	canvas->SetRightMargin(0.05f) ;
}

void addWIP(TCanvas* canvas, std::string title)
{
	canvas->cd() ;
	TText* t = new TText(.95 , .96 , title.c_str()) ;
	t->SetNDC() ;
	t->SetTextColor(kGray+1) ;
	t->SetTextAlign(31);
	t->SetTextSize(0.03f) ;
	t->Draw() ;
}

void addPolText(TCanvas* canvas , std::string pol)
{
	canvas->cd() ;
	TLatex* t = new TLatex(.10 , .96 , pol.c_str()) ;
	t->SetNDC() ;
	t->SetTextAlign(11) ;
	t->SetTextSize(0.03f) ;
	t->Draw() ;
}

void addEnergyText(TCanvas* canvas , double energy)
{
	canvas->cd() ;
	std::stringstream toto ;
	toto << "#sqrt{s} = " << energy << " GeV" ;
	TLatex* t = new TLatex(.50 , .96 , toto.str().c_str()) ;
	t->SetNDC() ;
	t->SetTextAlign(21) ;
	t->SetTextSize(0.03f) ;
	t->Draw() ;
}

void addLumiText(TCanvas* canvas, std::string lumi)
{
	canvas->cd() ;
	std::stringstream toto ;
	toto << "L_{int} = " << lumi.c_str() << " fb^{-1}" ;
	TLatex* latex = new TLatex(.10 , .03 , toto.str().c_str()) ;
	latex->SetNDC() ;
	latex->SetTextSize(0.03f) ;
	latex->SetTextAlign(11) ;
	latex->Draw() ;
}

TCanvas* drawCanvas(std::string canvasName , TH2D* histo , std::string polText , std::string lumi , double energy)
{
	TCanvas* canvas = new TCanvas(canvasName.c_str() , canvasName.c_str() , 1000 , 1000) ;
	canvas->cd() ;
	setStyle(canvas) ;

	histo->Draw("col") ;

	addWIP(canvas , lumi) ;
	addPolText(canvas , polText) ;
	addEnergyText(canvas , energy) ;

	return canvas ;
}

TCanvas* drawCanvas(std::string canvasName , std::map<int , TH1D*> histoMap , std::string polText , std::string lumi , std::string title , double energy , bool ZH = true , std::vector<double> lines = {})
{
	TCanvas* canvas = new TCanvas(canvasName.c_str() , canvasName.c_str() , 1000 , 1000) ;
	canvas->cd() ;
	setStyle(canvas) ;

	double maximum = std::numeric_limits<double>::min() ;
	TH1D* maxHisto = nullptr ;

	for (const auto& histo : histoMap)
	{
		double max = histo.second->GetBinContent( histo.second->GetMaximumBin() )/histo.second->Integral() ;
		if ( max > maximum )
		{
			maximum = max ;
			maxHisto = histo.second ;
		}
	}

	maxHisto->DrawNormalized("HIST") ;
	//histoMap.at("ZH")->DrawNormalized("HIST same") ;
	for ( const auto& histo : histoMap )
	{
		if ( histo.first == 0 )
			continue ;
		histo.second->DrawNormalized("HIST same") ;
	}
	histoMap.at(0)->DrawNormalized("HIST same") ;

	//legend
	TLegend* leg = new TLegend(0.65,0.65,0.88,0.88) ;
	leg->SetBorderSize(0) ;

	TLegendEntry* le = nullptr ;

	//le = leg->AddEntry(histoMap.at("ZH") , "ZH" , "f") ;
	//le->SetTextColor( histoMap.at("ZH")->GetLineColor() ) ;

	for ( const auto& histo : histoMap )
	{
		le = leg->AddEntry(histo.second , processNames[histo.first].c_str() , "f") ;
		le->SetTextColor( histo.second->GetLineColor() ) ;
	}
	leg->Draw() ;
	addWIP(canvas , title) ;
	addPolText(canvas , polText) ;
	addEnergyText(canvas , energy) ;
	addLumiText(canvas , lumi);

	std::array<double,4> bords = {} ;
	bords[0] = canvas->GetLeftMargin() ;
	bords[1] = canvas->GetBottomMargin() ;
	bords[2] = 1.0 - canvas->GetRightMargin() ;
	bords[3] = 1.0 - canvas->GetTopMargin() ;

	for ( const auto& val : lines )
	{
		double x = ( val - maxHisto->GetXaxis()->GetXmin() ) / (maxHisto->GetXaxis()->GetXmax() - maxHisto->GetXaxis()->GetXmin() ) ;
		x = x*(bords[2]-bords[0]) + bords[0] ;

		TLine* line = new TLine(x,bords[1],x,bords[3]) ;
		line->SetNDC() ;
		line->SetLineColor(kMagenta+2) ;
		line->SetLineWidth(2) ;
		line->SetLineStyle(7) ;
		line->Draw() ;
	}

	canvas->RedrawAxis() ;

	return canvas ;
}

/*TCanvas* drawCanvasStack(std::string canvasName , TH1D* stack , std::string polText , std::string lumi , std::string title , double energy , bool ZH = true , std::vector<double> lines = {})
{
	TCanvas* canvas = new TCanvas(canvasName.c_str() , canvasName.c_str() , 1000 , 1000) ;
	canvas->cd() ;
	setStyle(canvas) ;

	double maximum = std::numeric_limits<double>::min() ;
	TH1D* maxHisto = nullptr ;

	for (const auto& histo : stack)
	{
		double max = histo.second->GetBinContent( histo.second->GetMaximumBin() )/histo.second->Integral() ;
		if ( max > maximum )
		{
			maximum = max ;
			maxHisto = histo.second ;
		}
	}

	maxHisto->DrawNormalized("HIST") ;
	//histoMap.at("ZH")->DrawNormalized("HIST same") ;
	for ( const auto& histo : histoMap )
	{
		if ( histo.first == std::string("ZH") )
			continue ;
		histo.second->DrawNormalized("HIST same") ;
	}
	//if ( ZH )
		//histoMap.at("ZH")->DrawNormalized("HIST same") ;

	//legend
	TLegend* leg = new TLegend(0.65,0.65,0.88,0.88) ;
	leg->SetBorderSize(0) ;

	TLegendEntry* le = nullptr ;

	//le = leg->AddEntry(histoMap.at("ZH") , "ZH" , "f") ;
	//le->SetTextColor( histoMap.at("ZH")->GetLineColor() ) ;

	for ( const auto& histo : histoMap )
	{
		// ( histo.first.c_str() == std::string("ZH") )
			//continue ;
		le = leg->AddEntry(histo.second , histo.first.c_str() , "f") ;
		le->SetTextColor( histo.second->GetLineColor() ) ;
	}
	leg->Draw() ;
	addWIP(canvas , title) ;
	addPolText(canvas , polText) ;
	addEnergyText(canvas , energy) ;
	addLumiText(canvas , lumi) ;

	std::array<double,4> bords = {} ;
	bords[0] = canvas->GetLeftMargin() ;
	bords[1] = canvas->GetBottomMargin() ;
	bords[2] = 1.0 - canvas->GetRightMargin() ;
	bords[3] = 1.0 - canvas->GetTopMargin() ;

	for ( const auto& val : lines )
	{
		double x = ( val - maxHisto->GetXaxis()->GetXmin() ) / (maxHisto->GetXaxis()->GetXmax() - maxHisto->GetXaxis()->GetXmin() ) ;
		x = x*(bords[2]-bords[0]) + bords[0] ;

		TLine* line = new TLine(x,bords[1],x,bords[3]) ;
		line->SetNDC() ;
		line->SetLineColor(kMagenta+2) ;
		line->SetLineWidth(2) ;
		line->SetLineStyle(7) ;
		line->Draw() ;
	}

	canvas->RedrawAxis() ;

	return canvas ;
}*/

void drawSquare( TCanvas* canvas , double x1 , double y1 , double x2 , double y2 )
{
	canvas->cd() ;
	TLine* line = new TLine(x1,y1,x2,y1) ;
	line->SetLineColor(kGreen-4) ;
	line->SetLineWidth(2) ;

	line->Draw() ;

	line = new TLine(x1,y2,x2,y2) ;
	line->SetLineColor(kGreen-4) ;
	line->SetLineWidth(2) ;

	line->Draw() ;

	line = new TLine(x1,y1,x1,y2) ;
	line->SetLineColor(kGreen-4) ;
	line->SetLineWidth(2) ;

	line->Draw() ;

	line = new TLine(x2,y1,x2,y2) ;
	line->SetLineColor(kGreen-4) ;
	line->SetLineWidth(2) ;

	line->Draw() ;
}

void drawEllipse(TCanvas* canvas , double x , double y , double r)
{
	canvas->cd() ;

	TEllipse* ellipse = new TEllipse(x,y,r) ;
	ellipse->SetLineColor(kGreen-4) ;
	ellipse->SetLineWidth(2) ;
	ellipse->SetFillStyle(0) ;

	ellipse->Draw() ;
}

int main(int argc , char** argv)
{

	if (argc < 4)
	{
		std::cout << "ERROR : pfpdpfds" << std::endl;
		return 1;
	}
	
	const double ePol = std::atof(argv[1]);
	const double pPol = std::atof(argv[2]);
	const int integratedLumi = std::atoi(argv[3]);
	
	if ( std::abs(ePol) > 1 || std::abs(pPol) > 1)
	{
		std::cout << "ERROR : incorrect pol" << std::endl;
		return 1;
	}
	
	
	//	gStyle->SetPalette(99) ;
	const std::set<int> signal = { 402007 , 402008 } ;
	const std::set<int> higgs = { 402001 , 402002 , 402003 , 402004 , 402005, 402006 , 402009 , 402010 , 402011 , 402012 } ;
	const std::set<int> _2fList = { 500006 , 500008 , 500010 , 500012 } ;
	const std::set<int> _4fList = { 500070 , 500072 , 500074 , 500076 , 500082 , 500084 , 500098 , 500100 , 500101 , 500102 , 500103 , 500104 , 500105 , 500106 , 500107 , 500108 , 500110 , 					 500112 , 500113 , 500114 , 500115 , 500116 , 500117 , 500118 , 500119 , 500120 , 500122 , 500124 } ;


	gStyle->SetOptStat(0) ;
	std::string polText = "P(e^{-},e^{+}) = (" + std::string(argv[1]) + " , " + std::string(argv[2]) + ")" ;

	constexpr double energy = 250 ;

	double eR = 0.5*(ePol+1) ;
	double eL = 0.5*(1-ePol) ;

	double pR = 0.5*(pPol+1) ;
	double pL = 0.5*(1-pPol) ;

	std::map<Polarisation , double> polWeightMap = {} ;
	polWeightMap.insert( { LR , eL*pR } ) ;
	polWeightMap.insert( { RL , eR*pL } ) ;
	polWeightMap.insert( { LL , eL*pL } ) ;
	polWeightMap.insert( { RR , eR*pR } ) ;

	std::cout << "Luminosity : " << integratedLumi << std::endl ;
	std::cout << "Polarisation : (" << 100*ePol << "% , " << 100*pPol << "%)" << std::endl ;

	std::cout << "weightLR : " << polWeightMap.at(LR) << std::endl ;
	std::cout << "weightRL : " << polWeightMap.at(RL) << std::endl ;
	std::cout << "weightLL : " << polWeightMap.at(LL) << std::endl ;
	std::cout << "weightRR : " << polWeightMap.at(RR) << std::endl ;

	std::vector<int> processIDVec = {} ;

	for ( auto it = ChannelInfoMap.begin() ; it != ChannelInfoMap.end() ; ++it )
		processIDVec.push_back(it->first) ;

	const std::string dir = "/home/nicolas/Documents/ep->nnH_(H->WW->qqqq)/MyAnalysis/script";

	std::map<int , TH1D*> recoE4jetsHistoMap = {} ;
	std::map<int , TH1D*> recoPt4jetsHistoMap = {} ;
	std::map<int , TH1D*> recoM4jetsHistoMap = {} ;
	std::map<int , TH1D*> recoWBigMassEnergyHistoMap = {} ;
    	std::map<int , TH1D*> recoWBigMassPtHistoMap = {} ;
	std::map<int , TH1D*> recoWBigMassMassHistoMap = {} ;
	std::map<int , TH1D*> recoWBigMassCosJetsHistoMap = {} ;
	std::map<int , TH1D*> recoWSmallMassEnergyHistoMap = {} ;
	std::map<int , TH1D*> recoWSmallMassPtHistoMap = {} ;
	std::map<int , TH1D*> recoWSmallMassMassHistoMap = {} ;
	std::map<int , TH1D*> recoWSmallMassCosJetsHistoMap = {} ;
	std::map<int , TH1D*> recoCosHistoMap = {} ;
	std::map<int , TH1D*> reconjetsHistoMap = {} ;
	std::map<int , TH1D*> recoY12HistoMap = {} ;
	std::map<int , TH1D*> recoY13HistoMap = {} ;
	std::map<int , TH1D*> recoY14HistoMap = {} ;
	std::map<int , TH1D*> recoY23HistoMap = {} ;
	std::map<int , TH1D*> recoY24HistoMap = {} ;
	std::map<int , TH1D*> recoY34HistoMap = {} ;
	
	THStack* recoM4jetsStack = new THStack( "hs3" , ";m (GeV);#events" ) ;


	for ( const auto& process : processMap )
	{
		std::stringstream toto ; toto << process.first << "_e4jets" ;
		TH1D* recoE4jets = new TH1D( toto.str().c_str() , ";E (GeV);#events" , 250 , 0 , 250 ) ;
		recoE4jets->SetDirectory(nullptr) ;
		recoE4jets->SetLineColor(process.second.color) ;

		std::stringstream toto2 ; toto2 << process.first << "_pt4jets" ;
		TH1D* recoPt4jets = new TH1D( toto2.str().c_str() , ";pt (GeV);#events" , 150 , 0 , 150 ) ;
		recoPt4jets->SetDirectory(nullptr) ;
		recoPt4jets->SetLineColor(process.second.color) ;
		
		std::stringstream toto3 ; toto3 << process.first << "_m4jets" ;
		TH1D* recoM4jets = new TH1D( toto3.str().c_str() , ";m (GeV);#events" , 250 , 0 , 250 ) ;
		recoM4jets->SetDirectory(nullptr) ;
		recoM4jets->SetLineColor(process.second.color) ;
		
		std::stringstream toto4 ; toto4 << process.first << "_wbigmassenergy" ;
		TH1D* recoWBigMassEnergy = new TH1D( toto4.str().c_str() , ";E (GeV);#events" , 250 , 0 , 250 ) ;
		recoWBigMassEnergy->SetDirectory(nullptr) ;
		recoWBigMassEnergy->SetLineColor(process.second.color) ;
		
		std::stringstream toto5 ; toto5 << process.first << "_wbigmasspt" ;
		TH1D* recoWBigMassPt = new TH1D( toto5.str().c_str() , ";pt (GeV);#events" , 150 , 0 , 150 ) ;
		recoWBigMassPt->SetDirectory(nullptr) ;
		recoWBigMassPt->SetLineColor(process.second.color) ;
		
		std::stringstream toto6 ; toto6 << process.first << "_wbigmassmass" ;
		TH1D* recoWBigMassMass = new TH1D( toto6.str().c_str() , ";m (GeV);#events" , 250 , 0 , 250 ) ;
		recoWBigMassMass->SetDirectory(nullptr) ;
		recoWBigMassMass->SetLineColor(process.second.color) ;
		
		std::stringstream toto7 ; toto7 << process.first << "_wbigmasscosjets" ;
		TH1D* recoWBigMassCosJets = new TH1D( toto7.str().c_str() , ";cos;#events" , 100 , -1 , 1 ) ;
		recoWBigMassCosJets->SetDirectory(nullptr) ;
		recoWBigMassCosJets->SetLineColor(process.second.color) ;
		
		std::stringstream toto8 ; toto8 << process.first << "_wsmallmassenergy" ;
		TH1D* recoWSmallMassEnergy = new TH1D( toto8.str().c_str() , ";E (GeV);#events" , 250 , 0 , 250 ) ;
		recoWSmallMassEnergy->SetDirectory(nullptr) ;
		recoWSmallMassEnergy->SetLineColor(process.second.color) ;
		
		std::stringstream toto9 ; toto9 << process.first << "_wsmallmasspt" ;
		TH1D* recoWSmallMassPt = new TH1D( toto9.str().c_str() , ";pt (GeV);#events" , 150 , 0 , 150 ) ;
		recoWSmallMassPt->SetDirectory(nullptr) ;
		recoWSmallMassPt->SetLineColor(process.second.color) ;
		
		std::stringstream toto10 ; toto10 << process.first << "_wsmallmassmass" ;
		TH1D* recoWSmallMassMass = new TH1D( toto10.str().c_str() , ";m (GeV);#events" , 250 , 0 , 250 ) ;
		recoWSmallMassMass->SetDirectory(nullptr) ;
		recoWSmallMassMass->SetLineColor(process.second.color) ;
		
		std::stringstream toto11 ; toto11 << process.first << "_wsmallmasscosjets" ;
		TH1D* recoWSmallMassCosJets = new TH1D( toto11.str().c_str() , ";cos;#events" , 100 , -1 , 1 ) ;
		recoWSmallMassCosJets->SetDirectory(nullptr) ;
		recoWSmallMassCosJets->SetLineColor(process.second.color) ;
		
		std::stringstream toto12 ; toto12 << process.first << "_cos" ;
		TH1D* recoCos = new TH1D( toto12.str().c_str() , ";cos;#events" , 100 , -1 , 1 ) ;
		recoCos->SetDirectory(nullptr) ;
		recoCos->SetLineColor(process.second.color) ;
		
		std::stringstream toto13 ; toto13 << process.first << "_njets" ;
		TH1D* reconjets = new TH1D( toto13.str().c_str() , ";#jets;#events" , 12 , 0 , 12 ) ;
		reconjets->SetDirectory(nullptr) ;
		reconjets->SetLineColor(process.second.color) ;
		
		std::stringstream toto14 ; toto14 << process.first << "_y12" ;
		TH1D* recoY12 = new TH1D( toto14.str().c_str() , ";-ln(Y12);#events" , 100 , 0 , 20 ) ;
		recoY12->SetDirectory(nullptr) ;
		recoY12->SetLineColor(process.second.color) ;
		
		std::stringstream toto15 ; toto15 << process.first << "_y13" ;
		TH1D* recoY13 = new TH1D( toto15.str().c_str() , ";-ln(Y13);#events" , 100 , 0 , 20 ) ;
		recoY13->SetDirectory(nullptr) ;
		recoY13->SetLineColor(process.second.color) ;
		
		std::stringstream toto16 ; toto16 << process.first << "_y14" ;
		TH1D* recoY14 = new TH1D( toto16.str().c_str() , ";-ln(Y14);#events" , 100 , 0 , 20 ) ;
		recoY14->SetDirectory(nullptr) ;
		recoY14->SetLineColor(process.second.color) ;
		
		std::stringstream toto17 ; toto17 << process.first << "_y23" ;
		TH1D* recoY23 = new TH1D( toto17.str().c_str() , ";-ln(Y23);#events" , 100 , 0 , 20 ) ;
		recoY23->SetDirectory(nullptr) ;
		recoY23->SetLineColor(process.second.color) ;
		
		std::stringstream toto18 ; toto18 << process.first << "_y24" ;
		TH1D* recoY24 = new TH1D( toto18.str().c_str() , ";-ln(Y24);#events" , 100 , 0 , 20 ) ;
		recoY24->SetDirectory(nullptr) ;
		recoY24->SetLineColor(process.second.color) ;
		
		std::stringstream toto19 ; toto19 << process.first << "_y34" ;
		TH1D* recoY34 = new TH1D( toto19.str().c_str() , ";-ln(Y34);#events" , 100 , 0 , 20 ) ;
		recoY34->SetDirectory(nullptr) ;
		recoY34->SetLineColor(process.second.color) ;

/*
		if ( process.first == std::string("ZH") )
		{
			histoRecMassFinal->SetLineColor(process.second) ;
		}
		else
		{
			histoRecMassFinal->SetFillColor(process.second) ;
			histoRecMassFinal->SetLineColor(kBlack) ;
		}
*/

		setStyle(recoE4jets) ;
		setStyle(recoPt4jets) ;
		setStyle(recoM4jets) ;
		setStyle(recoWBigMassEnergy) ;
		setStyle(recoWBigMassPt) ;
		setStyle(recoWBigMassMass) ;
		setStyle(recoWBigMassCosJets) ;
		setStyle(recoWSmallMassEnergy) ;
		setStyle(recoWSmallMassPt) ;
		setStyle(recoWSmallMassMass) ;
		setStyle(recoWSmallMassCosJets) ;
		setStyle(recoCos) ;
		setStyle(reconjets) ;
		setStyle(recoY12) ;
		setStyle(recoY13) ;
		setStyle(recoY14) ;
		setStyle(recoY23) ;
		setStyle(recoY24) ;
		setStyle(recoY34) ;

		recoE4jetsHistoMap.insert( {process.first , recoE4jets} ) ;
		recoPt4jetsHistoMap.insert( {process.first , recoPt4jets} ) ;
		recoM4jetsHistoMap.insert( {process.first , recoM4jets} ) ;
		recoWBigMassEnergyHistoMap.insert( {process.first , recoWBigMassEnergy} ) ;
		recoWBigMassPtHistoMap.insert( {process.first , recoWBigMassPt} ) ;
		recoWBigMassMassHistoMap.insert( {process.first , recoWBigMassMass} ) ;
		recoWBigMassCosJetsHistoMap.insert( {process.first , recoWBigMassCosJets} ) ;
		recoWSmallMassEnergyHistoMap.insert( {process.first , recoWSmallMassEnergy} ) ;
		recoWSmallMassPtHistoMap.insert( {process.first , recoWSmallMassPt} ) ;
		recoWSmallMassMassHistoMap.insert( {process.first , recoWSmallMassMass} ) ;
		recoWSmallMassCosJetsHistoMap.insert( {process.first , recoWSmallMassCosJets} ) ;
		recoCosHistoMap.insert( {process.first , recoCos} ) ;
		reconjetsHistoMap.insert( {process.first , reconjets} ) ;
		recoY12HistoMap.insert( {process.first , recoY12} ) ;
		recoY13HistoMap.insert( {process.first , recoY13} ) ;
		recoY14HistoMap.insert( {process.first , recoY14} ) ;
		recoY23HistoMap.insert( {process.first , recoY23} ) ;
		recoY24HistoMap.insert( {process.first , recoY24} ) ;
		recoY34HistoMap.insert( {process.first , recoY34} ) ;
	}

	
	float reco_E4jets = 0;
	float reco_Pt4jets = 0;
	float reco_M4jets = 0;
	float reco_WBigMass_Energy = 0;
	float reco_WBigMass_Pt = 0;
	float reco_WBigMass_Mass = 0;
	float reco_WBigMass_CosJets = 0;
	float reco_WSmallMass_Energy = 0;
	float reco_WSmallMass_Pt = 0;
	float reco_WSmallMass_Mass = 0;
	float reco_WSmallMass_CosJets = 0;
	float reco_Cos = 0;
	int reco_njets = 0;
	float reco_Y12 = 0;
	float reco_Y13 = 0;
	float reco_Y14 = 0;
	float reco_Y23 = 0;
	float reco_Y24 = 0;
	float reco_Y34 = 0;
	int mc_HiggsDecay = 0;
	int mc_DecayProcess = 0;
	
	for ( const int& processID : processIDVec )
	{
		std::cout << "Processing " << processID << "..." << std::endl ;

		std::stringstream fileName ;
		fileName << dir << "/" << processID << ".root" ;
		TFile* file = TFile::Open( fileName.str().c_str() ) ;
		if ( !file )
		{
			std::cerr << "ERROR : could not open data file" << std::endl ;
			return 1 ;
		}
		TTree* tree = dynamic_cast<TTree*>( file->Get("tree") ) ;
		if ( !tree )
		{
			std::cerr << "ERROR : tree not present" << std::endl ;
			return 1 ;
		}

		tree->SetBranchAddress("reco_E4jets", &reco_E4jets);
		tree->SetBranchAddress("reco_Pt4jets", &reco_Pt4jets);
		tree->SetBranchAddress("reco_M4jets", &reco_M4jets);
		tree->SetBranchAddress("reco_WBigMass_Energy", &reco_WBigMass_Energy);
		tree->SetBranchAddress("reco_WBigMass_Pt", &reco_WBigMass_Pt);
		tree->SetBranchAddress("reco_WBigMass_Mass", &reco_WBigMass_Mass);
		tree->SetBranchAddress("reco_WBigMass_CosJets", &reco_WBigMass_CosJets);
		tree->SetBranchAddress("reco_WSmallMass_Energy", &reco_WSmallMass_Energy);
		tree->SetBranchAddress("reco_WSmallMass_Pt", &reco_WSmallMass_Pt);
		tree->SetBranchAddress("reco_WSmallMass_Mass", &reco_WSmallMass_Mass);
		tree->SetBranchAddress("reco_WSmallMass_CosJets", &reco_WSmallMass_CosJets);
		tree->SetBranchAddress("reco_Cos", &reco_Cos);
		tree->SetBranchAddress("reco_njets", &reco_njets);
		tree->SetBranchAddress("reco_Y12", &reco_Y12);
		tree->SetBranchAddress("reco_Y13", &reco_Y13);
		tree->SetBranchAddress("reco_Y14", &reco_Y14);
		tree->SetBranchAddress("reco_Y23", &reco_Y23);
		tree->SetBranchAddress("reco_Y24", &reco_Y24);
		tree->SetBranchAddress("reco_Y34", &reco_Y34);
		if ( signal.count(processID) )
		{
			tree->SetBranchAddress("mc_HiggsDecay", &mc_HiggsDecay);
			tree->SetBranchAddress("mc_DecayProcess", &mc_DecayProcess);
		}

		double nEventsExpected = ChannelInfoMap.at(processID).xSect*integratedLumi ;
		nEventsExpected *= polWeightMap.at( ChannelInfoMap.at(processID).pol ) ;
		double weight = nEventsExpected/tree->GetEntries() ;

		int iEntry = 0 ;
		while( tree->GetEntry( iEntry++ ) )
		{
			bool isSignal = false;
			if (signal.count(processID))
			{
				if (mc_HiggsDecay == 2424)
					if (mc_DecayProcess == 808)
						isSignal = true; 
			}
			auto process = Process::getProcess(processID, isSignal) ;
			recoE4jetsHistoMap.at( process )->Fill(reco_E4jets , weight) ;
			recoPt4jetsHistoMap.at( process )->Fill(reco_Pt4jets , weight) ;
			recoM4jetsHistoMap.at( process )->Fill(reco_M4jets , weight) ;
			recoWBigMassEnergyHistoMap.at( process )->Fill(reco_WBigMass_Energy , weight) ;
			recoWBigMassPtHistoMap.at( process )->Fill(reco_WBigMass_Pt , weight) ;
			recoWBigMassMassHistoMap.at( process )->Fill(reco_WBigMass_Mass , weight) ;
			recoWBigMassCosJetsHistoMap.at( process )->Fill(reco_WBigMass_CosJets , weight) ;
			recoWSmallMassEnergyHistoMap.at( process )->Fill(reco_WSmallMass_Energy , weight) ;
			recoWSmallMassPtHistoMap.at( process )->Fill(reco_WSmallMass_Pt , weight) ;
			recoWSmallMassMassHistoMap.at( process )->Fill(reco_WSmallMass_Mass , weight) ;
			recoWSmallMassCosJetsHistoMap.at( process )->Fill(reco_WSmallMass_CosJets , weight) ;
			recoCosHistoMap.at( process )->Fill(reco_Cos , weight) ;
			reconjetsHistoMap.at( process )->Fill(reco_njets , weight) ;
			recoY12HistoMap.at( process )->Fill(reco_Y12 , weight) ; 
			recoY13HistoMap.at( process )->Fill(reco_Y13 , weight) ;
			recoY14HistoMap.at( process )->Fill(reco_Y14 , weight) ;
			recoY23HistoMap.at( process )->Fill(reco_Y23 , weight) ;
			recoY24HistoMap.at( process )->Fill(reco_Y24 , weight) ;
			recoY34HistoMap.at( process )->Fill(reco_Y34 , weight) ;
		}

		file->Close() ;
	}
	//std::vector<std::string> order = {"1, signal" , "2, Other Higgs decays" , "3, Higgsstrahlung other Z decays" , "4, 4 fermions" , "5, 2 fermions"} ;
	for ( const auto& [process , osef] : processMap)
	{ 
		recoM4jetsStack->Add(recoM4jetsHistoMap.at( process )) ;
	}
	TCanvas* c1000 = new TCanvas("c1000" , "c1000" , 1000 , 1000) ;
	c1000->cd() ;
	recoM4jetsStack->Draw("HIST") ;

	std::array<double,2> recoE4jetsCuts = {126.418,175.326};
	std::array<double,2> recoPt4jetsCuts = {11.4895,74.3785};
	std::array<double,2> recoM4jetsCuts = {116.87,162.354};
	std::array<double,2> recoWBigMassEnergyCuts = {73.2245,117.171};
	std::array<double,2> recoWBigMassPtCuts = {8.84967,72.5748};
	std::array<double,2> recoWBigMassMassCuts = {56.8366,95.8543};
	std::array<double,2> recoWBigMassCosJetsCuts = {-0.971238,0.226284};
	std::array<double,2> recoWSmallMassEnergyCuts = {30.0771,77.6923};
	std::array<double,2> recoWSmallMassPtCuts = {5.50957,50.2233};
	std::array<double,2> recoWSmallMassMassCuts = {16.3641,64.2769};
	std::array<double,2> recoWSmallMassCosJetsCuts = {-0.947432,0.807242};
	std::array<double,2> recoCosCuts = {-0.907282,0.897197};
	std::array<double,2> recoY12Cuts = {1.76816,4.3256};
	std::array<double,2> recoY13Cuts = {2.05077,4.77881};
	std::array<double,2> recoY14Cuts = {2.41655,5.24884};
	std::array<double,2> recoY23Cuts = {2.32113,5.31321};
	std::array<double,2> recoY24Cuts = {2.78922,5.64597};
	std::array<double,2> recoY34Cuts = {3.01806,6.07034};

	TCanvas* c1 = drawCanvas("c1" , recoE4jetsHistoMap , polText , argv[3] , "reco_E4jets" , energy ,  false , { recoE4jetsCuts[0] , recoE4jetsCuts[1] }) ;
	TCanvas* c2 = drawCanvas("c2" , recoPt4jetsHistoMap , polText , argv[3] , "reco_Pt4jets" , energy , false , { recoPt4jetsCuts[0] ,recoPt4jetsCuts[1] }) ;
	TCanvas* c3 = drawCanvas("c3" , recoM4jetsHistoMap , polText , argv[3] , "reco_M4jets" , energy , false , { recoM4jetsCuts[0] , recoM4jetsCuts[1] }) ;
	TCanvas* c4 = drawCanvas("c4" , recoWBigMassEnergyHistoMap , polText , argv[3] , "reco_WBigMass_Energy" , energy , false , { recoWBigMassEnergyCuts[0] , recoWBigMassEnergyCuts[1] }) ;
	TCanvas* c5 = drawCanvas("c5" , recoWBigMassPtHistoMap , polText , argv[3] , "reco_WBigMass_Pt" , energy , false , { recoWBigMassPtCuts[0] , recoWBigMassPtCuts[1] }) ;
	TCanvas* c6 = drawCanvas("c6" , recoWBigMassMassHistoMap , polText , argv[3] , "reco_WBigMass_Mass" , energy , false , { recoWBigMassMassCuts[0] , recoWBigMassMassCuts[1] }) ;
	TCanvas* c7 = drawCanvas("c7" , recoWBigMassCosJetsHistoMap , polText , argv[3] , "reco_WBigMass_CosJets" , energy , false , { recoWBigMassCosJetsCuts[0] , recoWBigMassCosJetsCuts[1] }) ;
	TCanvas* c8 = drawCanvas("c8" , recoWSmallMassEnergyHistoMap , polText , argv[3] , "reco_WSmallMass_Energy" , energy , false , { recoWSmallMassEnergyCuts[0] , recoWSmallMassEnergyCuts[1] }) ;
	TCanvas* c9 = drawCanvas("c9" , recoWSmallMassPtHistoMap , polText , argv[3] , "reco_WSmallMass_Pt" , energy , false , { recoWSmallMassPtCuts[0] , recoWSmallMassPtCuts[1] }) ;
	TCanvas* c10 = drawCanvas("c10" , recoWSmallMassMassHistoMap , polText , argv[3] , "reco_WSmallMass_Mass" , energy , false , { recoWSmallMassMassCuts[0] , recoWSmallMassMassCuts[1] }) ;
	TCanvas* c11 = drawCanvas("c11" , recoWSmallMassCosJetsHistoMap , polText , argv[3] , "reco_WSmallMass_CosJets" , energy , false , { recoWSmallMassCosJetsCuts[0] , recoWSmallMassCosJetsCuts[1] }) ;
	TCanvas* c12 = drawCanvas("c12" , recoCosHistoMap , polText , argv[3] , "reco_Cos" , energy , false , { recoCosCuts[0] , recoCosCuts[1] }) ;
	TCanvas* c13 = drawCanvas("c13" , reconjetsHistoMap , polText , argv[3] , "reco_njets" , energy , false , {4,7}) ;
	TCanvas* c14 = drawCanvas("c14" , recoY12HistoMap , polText , argv[3] , "reco_Y12" , energy , false , { recoY12Cuts[0] , recoY12Cuts[1] }) ;
	TCanvas* c15 = drawCanvas("c15" , recoY13HistoMap , polText , argv[3] , "reco_Y13" , energy , false , { recoY13Cuts[0] , recoY13Cuts[1] }) ;
	TCanvas* c16 = drawCanvas("c16" , recoY14HistoMap , polText , argv[3] , "reco_Y14" , energy , false , { recoY14Cuts[0] , recoY14Cuts[1] }) ;
	TCanvas* c17 = drawCanvas("c17" , recoY23HistoMap , polText , argv[3] , "reco_Y23" , energy , false , { recoY23Cuts[0] , recoY23Cuts[1] }) ;
	TCanvas* c18 = drawCanvas("c18" , recoY24HistoMap , polText , argv[3] , "reco_Y24" , energy , false , { recoY24Cuts[0] , recoY24Cuts[1] }) ;
	TCanvas* c19 = drawCanvas("c19" , recoY34HistoMap , polText , argv[3] , "reco_Y34" , energy , false , { recoY34Cuts[0] , recoY34Cuts[1] }) ;

/*
	THStack* hs = new THStack("hs",";m_{recoil} (GeV);nEvents/2.5GeV") ;

	std::vector< std::string > chanVec = {"ZH" , "q#bar{q}" , "q#bar{q}q#bar{q}" , "q#bar{q}ll" , "q#bar{q}l#nu" , "q#bar{q}#nu#nu"} ;
	auto chanVecRev = chanVec ;
	std::reverse(chanVecRev.begin() , chanVecRev.end() ) ;

	for ( const auto& str : chanVecRev )
		hs->Add( recMassFinalHistoMap.at(str) ) ;

	hs->Draw("hist") ;
	setStyle(hs) ;

	c7->Modified() ;

	//legend
	TLegend* leg2 = new TLegend(0.65,0.65,0.88,0.88) ;
	leg2->SetBorderSize(0) ;
	for ( const auto& str : chanVec )
	{
		auto histo = recMassFinalHistoMap.at(str) ;
		auto le2 = leg2->AddEntry(histo , str.c_str() , "f") ;
		le2->SetTextColor( processMap.at(str) ) ;
	}
	leg2->Draw() ;

	addWIP(c7) ;
	addPolText(c7,polText) ;
	addLumiText(c7) ;
	addEnergyText(c7 , 250) ;
*/

	TFile* outputFile = new TFile("test.root" , "RECREATE") ;
	outputFile->cd() ;

	c1->Write("c1") ;
	c2->Write("c2") ;
	c3->Write("c3") ;
	c4->Write("c4") ;
	c5->Write("c5") ;
	c6->Write("c6") ;
	c7->Write("c7") ;
	c8->Write("c8") ;
	c9->Write("c9") ;
	c10->Write("c10") ;
	c11->Write("c11") ;
	c12->Write("c12") ;
	c13->Write("c13") ;
	c14->Write("c14") ;
	c15->Write("c15") ;
	c16->Write("c16") ;
	c17->Write("c17") ;
	c18->Write("c18") ;
	c19->Write("c19") ;
	c1000->Write("c1000") ;

	outputFile->Close() ;

	return 0 ;
}
