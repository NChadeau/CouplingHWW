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
#include "TGraph.h"


std::map<double, std::array<double,2>> getQuantiles(TH1* histo)
{
	const double begin = 0;
	const double end = 0.15;
	const int nPoints = 150;
	const double step = (end-begin)/nPoints;
	std::array<double, nPoints*2> quantiles {};
	for(auto i = 0 ; i < nPoints ; i++)
	{
		quantiles[i] = begin + 0.5*i*step;
		quantiles[nPoints*2-1-i] = 1-0.5*i*step;
	}

	std::array<double, nPoints*2> results {};
	histo->GetQuantiles(quantiles.size(), results.data(), quantiles.data());
	
	std::map<double,std::array<double,2>> toReturn {};
	for(auto i = 0 ; i < nPoints ; i++)
	{
		toReturn[begin+i*step] = {results[i], results[nPoints*2-1-i]};
	}
	
	return toReturn;
}

struct processInfos
{
	std::string displayName = "";
	Color_t color = kBlack;
};

static const std::map<std::string , Color_t> processMap =
{
	{"1, signal" , kBlack } ,
	{"2, Other Higgs decays" , kCyan+2 } , 
	{"3, Higgsstrahlung other Z decays" , kCyan} ,
	{"5, 2 fermions" , kRed-3 } ,
	{"4, 4 fermions" , kRed-7 } 
} ;

/*
static const std::map<int , processInfos> processMap =
{
	{1,{"1, signal" , kBlack }} ,
	{2,{"2, Other Higgs decays" , kCyan+2 }} , 
	{3,{"3, Higgsstrahlung other Z decays" , kCyan}} ,
	{4,{"5, 2 fermions" , kRed-3 }} ,
	{5,{"4, 4 fermions" , kRed-7 }}
} ;
*/

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
	TLatex* latex = new TLatex(.85 , .96 , toto.str().c_str()) ;
	latex->SetNDC() ;
	latex->SetTextSize(0.03f) ;
	latex->SetTextAlign(11) ;
	latex->Draw() ;
}

TCanvas* drawCanvas(std::string canvasName , TGraph* graph , std::string polText , std::string lumi , double energy, std::string axis)
{
	TCanvas* canvas = new TCanvas(canvasName.c_str() , canvasName.c_str() , 1000 , 1000) ;
	canvas->cd() ;
	setStyle(canvas) ;

	graph->GetYaxis()->SetTitle( axis.c_str() ) ;
	graph->GetXaxis()->SetTitle( "quantile (%)" ) ;
	graph->Draw("AP") ;

	addLumiText(canvas , lumi) ;
	addPolText(canvas , polText) ;
	addEnergyText(canvas , energy) ;

	return canvas ;
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
	const std::set<int> signal = {402007 , 402008} ;
	const std::set<int> higgs = {402001 , 402002 , 402003 , 402004 , 402005 , 402006 , 402009 , 402010 , 402011 , 402012} ;
	const std::set<int> _2fList = {500006 , 500008 , 500010 , 500012} ;
	const std::set<int> _4fList = {500070 , 500072 , 500074 , 500076 , 500082 , 500084 , 500098 , 500100 , 500101 , 500102 , 500103 , 500104 , 500105 , 500106 , 500107 , 500108 , 500110 , 					 500112 , 500113 , 500114 , 500115 , 500116 , 500117 , 500118 , 500119 , 500120 , 500122 , 500124} ;


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
	for (const auto& processID : signal)
		processIDVec.push_back(processID) ;	
	for (const auto& processID : higgs)
		processIDVec.push_back(processID) ;	
	for (const auto& processID : _2fList)
		processIDVec.push_back(processID) ;	
	for (const auto& processID : _4fList)
		processIDVec.push_back(processID) ;		

	//for ( auto it = ChannelInfoMap.begin() ; it != ChannelInfoMap.end() ; ++it )
		//processIDVec.push_back(it->first) ;

	const std::string dir = "/home/nicolas/Documents/ep->nnH_(H->WW->qqqq)/MyAnalysis/script";

		TH1D* recoE4jets = new TH1D( "_e4jets" , ";E (GeV);#events" , 250 , 0 , 250 ) ;
		recoE4jets->SetDirectory(nullptr) ;

		TH1D* recoPt4jets = new TH1D( "_pt4jets" , ";pt (GeV);#events" , 150 , 0 , 150 ) ;
		recoPt4jets->SetDirectory(nullptr) ;

		TH1D* recoM4jets = new TH1D( "_m4jets" , ";m (GeV);#events" , 250 , 0 , 250 ) ;
		recoM4jets->SetDirectory(nullptr) ;

		TH1D* recoWBigMassEnergy = new TH1D( "_wbigmassenergy" , ";E (GeV);#events" , 250 , 0 , 250 ) ;
		recoWBigMassEnergy->SetDirectory(nullptr) ;

		TH1D* recoWBigMassPt = new TH1D( "_wbigmasspt" , ";pt (GeV);#events" , 150 , 0 , 150 ) ;
		recoWBigMassPt->SetDirectory(nullptr) ;

		TH1D* recoWBigMassMass = new TH1D( "_wbigmassmass" , ";m (GeV);#events" , 250 , 0 , 250 ) ;
		recoWBigMassMass->SetDirectory(nullptr) ;

		TH1D* recoWBigMassCosJets = new TH1D( "_wbigmasscosjets" , ";cos;#events" , 100 , -1 , 1 ) ;
		recoWBigMassCosJets->SetDirectory(nullptr) ;

		TH1D* recoWSmallMassEnergy = new TH1D( "_wsmallmassenergy" , ";E (GeV);#events" , 250 , 0 , 250 ) ;
		recoWSmallMassEnergy->SetDirectory(nullptr) ;

		TH1D* recoWSmallMassPt = new TH1D( "_wsmallmasspt" , ";pt (GeV);#events" , 150 , 0 , 150 ) ;
		recoWSmallMassPt->SetDirectory(nullptr) ;

		TH1D* recoWSmallMassMass = new TH1D( "_wsmallmassmass" , ";m (GeV);#events" , 250 , 0 , 250 ) ;
		recoWSmallMassMass->SetDirectory(nullptr) ;

		TH1D* recoWSmallMassCosJets = new TH1D( "_wsmallmasscosjets" , ";cos;#events" , 100 , -1 , 1 ) ;
		recoWSmallMassCosJets->SetDirectory(nullptr) ;

		TH1D* recoCos = new TH1D( "_cos" , ";cos;#events" , 100 , -1 , 1 ) ;
		recoCos->SetDirectory(nullptr) ;

		TH1D* reconjets = new TH1D( "_njets" , ";#jets;#events" , 12 , 0 , 12 ) ;
		reconjets->SetDirectory(nullptr) ;

		TH1D* recoY12 = new TH1D( "_y12" , ";-ln(Y12);#events" , 100 , 0 , 20 ) ;
		recoY12->SetDirectory(nullptr) ;

		TH1D* recoY13 = new TH1D( "_y13" , ";-ln(Y13);#events" , 100 , 0 , 20 ) ;
		recoY13->SetDirectory(nullptr) ;

		TH1D* recoY14 = new TH1D( "_y14" , ";-ln(Y14);#events" , 100 , 0 , 20 ) ;
		recoY14->SetDirectory(nullptr) ;

		TH1D* recoY23 = new TH1D( "_y23" , ";-ln(Y23);#events" , 100 , 0 , 20 ) ;
		recoY23->SetDirectory(nullptr) ;

		TH1D* recoY24 = new TH1D( "_y24" , ";-ln(Y24);#events" , 100 , 0 , 20 ) ;
		recoY24->SetDirectory(nullptr) ;

		TH1D* recoY34 = new TH1D( "_y34" , ";-ln(Y34);#events" , 100 , 0 , 20 ) ;
		recoY34->SetDirectory(nullptr) ;
		
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
	
	for ( const int& processID : signal )
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
		tree->SetBranchAddress("mc_HiggsDecay", &mc_HiggsDecay);
		tree->SetBranchAddress("mc_DecayProcess", &mc_DecayProcess);

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
			
			if (!isSignal)
				continue;
	
			recoE4jets->Fill(reco_E4jets , weight) ;
			recoPt4jets->Fill(reco_Pt4jets , weight) ;
			recoM4jets->Fill(reco_M4jets , weight) ;
			recoWBigMassEnergy->Fill(reco_WBigMass_Energy , weight) ;
			recoWBigMassPt->Fill(reco_WBigMass_Pt , weight) ;
			recoWBigMassMass->Fill(reco_WBigMass_Mass , weight) ;
			recoWBigMassCosJets->Fill(reco_WBigMass_CosJets , weight) ;
			recoWSmallMassEnergy->Fill(reco_WSmallMass_Energy , weight) ;
			recoWSmallMassPt->Fill(reco_WSmallMass_Pt , weight) ;
			recoWSmallMassMass->Fill(reco_WSmallMass_Mass , weight) ;
			recoWSmallMassCosJets->Fill(reco_WSmallMass_CosJets , weight) ;
			recoCos->Fill(reco_Cos , weight) ;
			reconjets->Fill(reco_njets , weight) ;
			recoY12->Fill(reco_Y12 , weight) ; 
			recoY13->Fill(reco_Y13 , weight) ;
			recoY14->Fill(reco_Y14 , weight) ;
			recoY23->Fill(reco_Y23 , weight) ;
			recoY24->Fill(reco_Y24 , weight) ;
			recoY34->Fill(reco_Y34 , weight) ;
		}

		file->Close() ;
	}
	
	auto recoE4jetsCuts = getQuantiles(recoE4jets) ;
	auto recoPt4jetsCuts = getQuantiles(recoPt4jets);
	auto recoM4jetsCuts = getQuantiles(recoM4jets);
	auto recoWBigMassEnergyCuts = getQuantiles(recoWBigMassEnergy);
	auto recoWBigMassPtCuts = getQuantiles(recoWBigMassPt);
	auto recoWBigMassMassCuts = getQuantiles(recoWBigMassMass);
	auto recoWBigMassCosJetsCuts = getQuantiles(recoWBigMassCosJets);
	auto recoWSmallMassEnergyCuts = getQuantiles(recoWSmallMassEnergy);
	auto recoWSmallMassPtCuts = getQuantiles(recoWSmallMassPt);
	auto recoWSmallMassMassCuts = getQuantiles(recoWSmallMassMass);
	auto recoWSmallMassCosJetsCuts = getQuantiles(recoWSmallMassCosJets);
	auto recoCosCuts = getQuantiles(recoCos) ;
	auto recoY12Cuts = getQuantiles(recoY12) ; 
	auto recoY13Cuts = getQuantiles(recoY13) ;
	auto recoY14Cuts = getQuantiles(recoY14) ;
	auto recoY23Cuts = getQuantiles(recoY23) ;
	auto recoY24Cuts = getQuantiles(recoY24) ;
	auto recoY34Cuts = getQuantiles(recoY34) ;
	
	// percentage, nEvents
	std::map<double,double> nSelectedSignal = {};
	std::map<double,double> nSelectedHiggs1 = {};
	std::map<double,double> nSelectedHiggs2 = {};
	std::map<double,double> nSelected_2f = {};
	std::map<double,double> nSelected_4f = {};
	std::map<double,double> SignificanceMap = {};
	std::map<double,double> SignalEfficiencyMap = {};
	std::map<double,double> Higgs1EfficiencyMap = {};
	std::map<double,double> Higgs2EfficiencyMap = {};
	std::map<double,double> _2fEfficiencyMap = {};
	std::map<double,double> _4fEfficiencyMap = {};
	std::map<double,double> PurityMap = {};
	float nSignal = 0;
	float nHiggs1 = 0;
	float nHiggs2 = 0;
	float n_2f = 0;
	float n_4f = 0;
	
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
				
		
		auto applyCut = [&](const double percentage) -> bool
		{
		if (!(reco_E4jets > recoE4jetsCuts[percentage][0] && reco_E4jets < recoE4jetsCuts[percentage][1])) 
				return false;
		if (!(reco_Pt4jets > recoPt4jetsCuts[percentage][0] && reco_Pt4jets < recoPt4jetsCuts[percentage][1])) 
			return false;
		if (!(reco_M4jets > recoM4jetsCuts[percentage][0] && reco_M4jets < recoM4jetsCuts[percentage][1])) 
			return false;
		if (!(reco_WBigMass_Energy > recoWBigMassEnergyCuts[percentage][0] && reco_WBigMass_Energy < recoWBigMassEnergyCuts[percentage][1]))
			return false;
		if (!(reco_WBigMass_Pt > recoWBigMassPtCuts[percentage][0] && reco_WBigMass_Pt < recoWBigMassPtCuts[percentage][1]))
			return false;
		if (!(reco_WBigMass_Mass > recoWBigMassMassCuts[percentage][0] && reco_WBigMass_Mass < recoWBigMassMassCuts[percentage][1]))
			return false;
		if (!(reco_WBigMass_CosJets > recoWBigMassCosJetsCuts[percentage][0] && reco_WBigMass_CosJets < recoWBigMassCosJetsCuts[percentage][1]))
			return false;
		if (!(reco_WSmallMass_Energy > recoWSmallMassEnergyCuts[percentage][0] && reco_WSmallMass_Energy < recoWSmallMassEnergyCuts[percentage][1]))
			return false;
		if (!(reco_WSmallMass_Pt > recoWSmallMassPtCuts[percentage][0] && reco_WSmallMass_Pt < recoWSmallMassPtCuts[percentage][1]))
			return false;
		if (!(reco_WSmallMass_Mass > recoWSmallMassMassCuts[percentage][0] && reco_WSmallMass_Mass < recoWSmallMassMassCuts[percentage][1]))
			return false;
		if (!(reco_WSmallMass_CosJets > recoWSmallMassCosJetsCuts[percentage][0] && reco_WSmallMass_CosJets < recoWSmallMassCosJetsCuts[percentage][1]))
			return false;
		if (!(reco_Cos > recoCosCuts[percentage][0] && reco_Cos < recoCosCuts[percentage][1]))
			return false;
		if (!(reco_njets > 3 && reco_njets < 7))
			return false;
		if (!(reco_Y12 > recoY12Cuts[percentage][0] && reco_Y12 < recoY12Cuts[percentage][1]))
			return false;
		if (!(reco_Y13 > recoY13Cuts[percentage][0] && reco_Y13 < recoY13Cuts[percentage][1]))
			return false;
		if (!(reco_Y14 > recoY14Cuts[percentage][0] && reco_Y14 < recoY14Cuts[percentage][1]))
			return false;
		if (!(reco_Y23 > recoY23Cuts[percentage][0] && reco_Y23 < recoY23Cuts[percentage][1]))
			return false;
		if (!(reco_Y24 > recoY24Cuts[percentage][0] && reco_Y24 < recoY24Cuts[percentage][1]))
			return false;
		if (!(reco_Y34 > recoY34Cuts[percentage][0] && reco_Y34 < recoY34Cuts[percentage][1]))
			return false;
					
		return true;
		};
		
		int iEntry = 0 ;
		while( tree->GetEntry( iEntry++ ) )
		{
			bool isSignal = false;
			bool isHiggs1 = false;
			bool isHiggs2 = false;
			bool is_2f = false;
			bool is_4f = false;
			if ( signal.count(processID) )
			{
				if (mc_HiggsDecay == 2424 && mc_DecayProcess == 808)
				{
					isSignal = true;
					nSignal += weight;
				}
				else 
				{
					isHiggs1 = true;
					nHiggs1 += weight;
				}
			}
			else if ( higgs.count(processID) )
			{
				isHiggs2 = true;
				nHiggs2 += weight;
			}
			else if ( _2fList.count(processID) )
			{
				is_2f = true;
				n_2f += weight;
			}
			else if ( _4fList.count(processID) )
			{
				is_4f = true;
				n_4f += weight;
			}
			auto process = Process::getProcess(processID, isSignal) ;
			
			for (const auto& [percentage, osef] : recoE4jetsCuts)
			{
				if ( applyCut(percentage) )
				{
					if (isSignal)
						nSelectedSignal[percentage] += weight;
					else if (isHiggs1)
						nSelectedHiggs1[percentage] += weight;	
					else if (isHiggs2)
						nSelectedHiggs2[percentage] += weight;
					else if (is_2f)
						nSelected_2f[percentage] += weight;
					else if (is_4f)
						nSelected_4f[percentage] += weight;
				}
			}
		}

		file->Close() ;
	}
	
	for (const auto& [percentage, osef] : recoE4jetsCuts)
	{
		double S = nSelectedSignal[percentage]/std::sqrt(nSelectedSignal[percentage] + nSelectedHiggs1[percentage] + nSelectedHiggs2[percentage] + nSelected_2f[percentage] + nSelected_4f[percentage]);
		double eSignal = 100*nSelectedSignal[percentage]/nSignal;
		double eHiggs1 = 100*nSelectedHiggs1[percentage]/nHiggs1;
		double eHiggs2 = 100*nSelectedHiggs2[percentage]/nHiggs2;
		double e_2f = 100*nSelected_2f[percentage]/n_2f;
		double e_4f = 100*nSelected_4f[percentage]/n_4f;
		double p = 100*nSelectedSignal[percentage]/(nSelectedSignal[percentage] + nSelectedHiggs1[percentage] + nSelectedHiggs2[percentage] + nSelected_2f[percentage] + nSelected_4f[percentage]);
		SignificanceMap[percentage] = S;
		SignalEfficiencyMap[percentage] = eSignal;
		Higgs1EfficiencyMap[percentage] = eHiggs1;
		Higgs2EfficiencyMap[percentage] = eHiggs2;
		_2fEfficiencyMap[percentage] = e_2f;
		_4fEfficiencyMap[percentage] = e_4f;
		PurityMap[percentage] = p;
	}
	
	TGraph* SignificanceGraph = new TGraph();
	SignificanceGraph->SetMarkerStyle(20);
	TGraph* EfficiencyGraph = new TGraph();
	EfficiencyGraph->SetMarkerStyle(20);
	TGraph* PurityGraph = new TGraph();
	PurityGraph->SetMarkerStyle(20);
	
	double maxSignificance = -1;
	double maxPercentage = -1;
	for (const auto& [percentage, osef] : recoE4jetsCuts)
	{
		auto s = SignificanceMap[percentage];
		if (s > maxSignificance)
		{
			maxSignificance = s;
			maxPercentage = percentage;
		}
		SignificanceGraph->SetPoint(SignificanceGraph->GetN(), percentage*100, SignificanceMap[percentage]);
		EfficiencyGraph->SetPoint(EfficiencyGraph->GetN(), percentage*100, SignalEfficiencyMap[percentage]);
		PurityGraph->SetPoint(PurityGraph->GetN(), percentage*100, PurityMap[percentage]);
	}
	
	TCanvas* c1 = drawCanvas("c1" , SignificanceGraph , polText , argv[3] , energy , "significance") ;
	TCanvas* c2 = drawCanvas("c2" , EfficiencyGraph , polText , argv[3] , energy , "signal selection efficiency (%)") ;
	TCanvas* c3 = drawCanvas("c3" , PurityGraph , polText , argv[3] , energy, "purity (%)") ;

	TFile* outputFile = new TFile("testSigni.root" , "RECREATE") ;
	outputFile->cd() ;

	c1->Write("c1") ;
	c2->Write("c2") ;
	c3->Write("c3") ;

	outputFile->Close() ;
	
	
	
    std::cout << "std::array<double,2> recoE4jetsCuts = {" << recoE4jetsCuts[maxPercentage][0] << "," << recoE4jetsCuts[maxPercentage][1] << "};" << std::endl;  
    std::cout << "std::array<double,2> recoPt4jetsCuts = {" << recoPt4jetsCuts[maxPercentage][0] << "," << recoPt4jetsCuts[maxPercentage][1] << "};" << std::endl;
    std::cout << "std::array<double,2> recoM4jetsCuts = {" << recoM4jetsCuts[maxPercentage][0] << "," << recoM4jetsCuts[maxPercentage][1] << "};" << std::endl;  
    std::cout << "std::array<double,2> recoWBigMassEnergyCuts = {" << recoWBigMassEnergyCuts[maxPercentage][0] << "," << recoWBigMassEnergyCuts[maxPercentage][1] << "};" << std::endl;
    std::cout << "std::array<double,2> recoWBigMassPtCuts = {" << recoWBigMassPtCuts[maxPercentage][0] << "," << recoWBigMassPtCuts[maxPercentage][1] << "};" << std::endl;   
    std::cout << "std::array<double,2> recoWBigMassMassCuts = {" << recoWBigMassMassCuts[maxPercentage][0] << "," << recoWBigMassMassCuts[maxPercentage][1] << "};" << std::endl;   
    std::cout << "std::array<double,2> recoWBigMassCosJetsCuts = {" << recoWBigMassCosJetsCuts[maxPercentage][0] << "," << recoWBigMassCosJetsCuts[maxPercentage][1] << "};" << std::endl; 
    std::cout << "std::array<double,2> recoWSmallMassEnergyCuts = {" << recoWSmallMassEnergyCuts[maxPercentage][0] << "," << recoWSmallMassEnergyCuts[maxPercentage][1] << "};" << std::endl;  
    std::cout << "std::array<double,2> recoWSmallMassPtCuts = {" << recoWSmallMassPtCuts[maxPercentage][0] << "," << recoWSmallMassPtCuts[maxPercentage][1] << "};" << std::endl;  
    std::cout << "std::array<double,2> recoWSmallMassMassCuts = {" << recoWSmallMassMassCuts[maxPercentage][0] << "," << recoWSmallMassMassCuts[maxPercentage][1] << "};" << std::endl;   
    std::cout << "std::array<double,2> recoWSmallMassCosJetsCuts = {" << recoWSmallMassCosJetsCuts[maxPercentage][0] << "," << recoWSmallMassCosJetsCuts[maxPercentage][1] << "};" << std::endl;
    std::cout << "std::array<double,2> recoCosCuts = {" << recoCosCuts[maxPercentage][0] << "," << recoCosCuts[maxPercentage][1] << "};" << std::endl;
    std::cout << "std::array<double,2> recoY12Cuts = {" << recoY12Cuts[maxPercentage][0] << "," << recoY12Cuts[maxPercentage][1] << "};" << std::endl;
    std::cout << "std::array<double,2> recoY13Cuts = {" << recoY13Cuts[maxPercentage][0] << "," << recoY13Cuts[maxPercentage][1] << "};" << std::endl;
    std::cout << "std::array<double,2> recoY14Cuts = {" << recoY14Cuts[maxPercentage][0] << "," << recoY14Cuts[maxPercentage][1] << "};" << std::endl;
    std::cout << "std::array<double,2> recoY23Cuts = {" << recoY23Cuts[maxPercentage][0] << "," << recoY23Cuts[maxPercentage][1] << "};" << std::endl;
    std::cout << "std::array<double,2> recoY24Cuts = {" << recoY24Cuts[maxPercentage][0] << "," << recoY24Cuts[maxPercentage][1] << "};" << std::endl;
    std::cout << "std::array<double,2> recoY34Cuts = {" << recoY34Cuts[maxPercentage][0] << "," << recoY34Cuts[maxPercentage][1] << "};" << std::endl;
    
    std::cout << "quantile: " << maxPercentage*100 << "%" << std::endl;
    std::cout << "significance: " << maxSignificance << std::endl;
    std::cout << "evenements de signal attendus: " << nSignal << ", evenements de signal selectionnés: " << nSelectedSignal[maxPercentage] << std::endl;
    std::cout << "signal selection efficiency: " << SignalEfficiencyMap[maxPercentage] << "%" << std::endl;
    std::cout << "evenements de Higgs1 attendus: " << nHiggs1 << ", evenements de Higgs1 selectionnés: " << nSelectedHiggs1[maxPercentage] << std::endl;
    std::cout << "Higgs1 selection efficiency: " << Higgs1EfficiencyMap[maxPercentage] << "%" << std::endl;
    std::cout << "evenements de Higgs2 attendus: " << nHiggs2 << ", evenements de Higgs2 selectionnés: " << nSelectedHiggs2[maxPercentage] << std::endl;
    std::cout << "Higgs2 selection efficiency: " << Higgs2EfficiencyMap[maxPercentage] << "%" << std::endl;
    std::cout << "evenements de 2f attendus: " << n_2f << ", evenements de 2f selectionnés: " << nSelected_2f[maxPercentage] << std::endl;
    std::cout << "2f selection efficiency: " << _2fEfficiencyMap[maxPercentage] << "%" << std::endl;
    std::cout << "evenements de 4f attendus: " << n_4f << ", evenements de 4f selectionnés: " << nSelected_4f[maxPercentage] << std::endl;
    std::cout << "4f selection efficiency: " << _4fEfficiencyMap[maxPercentage] << "%" << std::endl;
    std::cout << "purity: " << PurityMap[maxPercentage] << "%" << std::endl;

	return 0 ;
}
