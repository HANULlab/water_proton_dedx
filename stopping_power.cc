// written by Hyunmin Yang, HANUL, Korea University
// This root macro draws stopping power of proton in water, reading data from SRIM.

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>

#include "TStyle.h"
#include "TAxis.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TF1.h"
#include "TBox.h"

using namespace std;
using namespace TMath;

double unitConversion(const string &unit)
{
    if(unit == "eV")
        return 1e-6;
    else if(unit == "keV")
        return 1e-3;
    else if(unit == "MeV") // Unit of energy : MeV
        return 1.;
    else if(unit == "GeV")
        return 1e3;
    else if(unit == "A")
        return 1e-8;
    else if(unit == "nm")
        return 1e-7;
    else if(unit == "um")
        return 1e-4;
    else if(unit == "mm")
        return 1e-1;
    else if(unit == "cm") // Unit of length : cm
        return 1.;
    else if(unit == "m")
        return 1e2;
    else
    {
        cerr << "Invalid unit name : " << unit << endl;
        return 0.;
    }
}


// not going to use this function because there are a misterious constant offset with Stopping data from SRIM.
double betheBloch(double *x, double *p)
{
    static constexpr double constant = 0.307075;
    static constexpr double me = 0.5109989461;
    const double Ek = x[0];
    const double m = p[0];
    const double z = p[1];

    const double I = p[2];
    const double Z = p[3];
    const double A = p[4];

    const double gamma = Ek/m + 1.;
    const double beta = Sqrt(gamma*gamma - 1.)/gamma;
    const double Wmax = 2*me*beta*beta*gamma*gamma/(1 + 2*gamma*me/m + me*me/(m*m));
    return constant*z*z*Z*
        (Log(2*me*beta*beta*gamma*gamma*Wmax/(I*I))/2 - beta*beta)
        /(A*beta*beta);
}

int stopping_power()
{

    gStyle->SetLineWidth(2);
    ifstream txtFile("elos_proton_water.txt", std::ios_base::in);
    stringstream ss;
    string line, unit;
    double value;
	// dEdX = edEdX(electrical losses) + ndEdX(nuclear losses)
    vector<double> Ek, edEdX, ndEdX, dEdX, range, lStrag, tStrag;

	// To limit ranges of the plots.
	int n = 0;    // total # of print
	
	// boundary of low and intermediate energy region
	int n_lm = 0;
	double Ek_lm = 0.025;
	
	// boundary of intermediate and high energy region
	int n_mh = 0;
	double Ek_mh = 2.;

	// loop for reading text file data from SRIM
    while(!txtFile.eof())
    {
        getline(txtFile, line);
        if(line.front() == '#')
            continue;
        else
        {
            auto pos = line.find_first_not_of(' ');
            if(pos == string::npos)
                continue;
            line = line.substr(pos);
        }
        ss.clear();
        ss.str(line);
		//--------------- Kinetic Energy ---------------------------
        ss >> value >> unit; Ek.push_back(value*unitConversion(unit));
		if(Ek.back() < Ek_lm)
            n_lm = n;
        else if(Ek.back() < Ek_mh)
            n_mh = n;
        // cout << value << " " << unit << " ";
       
		//--------------- Electrical stopping power ---------------------------
        ss >> value; edEdX.push_back(1000*value);
        // cout << value << " ";
		
		//--------------- Nuclear stopping power---------------------------
        ss >> value; ndEdX.push_back(1000*value);
        // cout << value << " ";
       
		//--------------- Total stopping power---------------------------
        dEdX.push_back(edEdX.back() + ndEdX.back());

		//--------------- Range ---------------------------
        ss >> value >> unit; range.push_back(value*unitConversion(unit));
        // cout << value << " " << unit << " ";
        
		//--------------- Longitudinal straggling ---------------------------
        ss >> value >> unit; lStrag.push_back(value*unitConversion(unit));
        // cout << value << " " << unit << " ";
        
		//--------------- Transversal straggling ---------------------------
        ss >> value >> unit; tStrag.push_back(value*unitConversion(unit));
        // cout << value << " " << unit << endl;
        ++n;
    }
    txtFile.close();

    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 600);

    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle("Mass stopping power of proton in liquid water;#it{E}_{k} [MeV];Mass Stopping Power[MeV cm^{2}/g]");
    
	mg->GetXaxis()->SetLimits(1e-3, 1e4);
    mg->GetXaxis()->SetTitleOffset(1.);mg->GetXaxis()->SetTitleSize(0.05);
    
	mg->GetYaxis()->SetRangeUser(1, 2.e3);
    mg->GetYaxis()->SetTitleOffset(0.6);mg->GetYaxis()->SetTitleSize(0.05);

    TGraph *g1 = new TGraph(n_mh, Ek.data(), dEdX.data()); // low and intermediate energy region 
    TGraph *g2 = new TGraph(n_lm, Ek.data(), edEdX.data()); // electrical losses in low energy region
    TGraph *g3 = new TGraph(n_mh, Ek.data(), ndEdX.data()); // nuclear losses in low and intermediate energy region 
	TGraph *g4 = new TGraph(n - n_mh, Ek.data() + n_mh, dEdX.data() + n_mh); // Bethe-Bloch region
    
	g1->SetLineWidth(3);g1->SetLineColor(kRed - 2);
    g2->SetLineWidth(2);g2->SetLineColor(6);g2->SetLineStyle(9);
    g3->SetLineWidth(2);g3->SetLineColor(29);g3->SetLineStyle(9);
    g4->SetLineWidth(3);g4->SetLineColor(kRed);

    mg->Add(g1);
    mg->Add(g2);
    mg->Add(g3);
    mg->Add(g4);

    c1->SetLogx();
    c1->SetLogy();
    c1->SetLeftMargin(0.08);
    c1->SetRightMargin(0.05);
    c1->SetBottomMargin(0.12);

    mg->Draw("ACSAME");
    
	// TF1 *bethe = new TF1("bethe", betheBloch, 2.1, 1e4, 5);
    // bethe->SetParameters(938.282, 1., 70.e-6, 7.42, 18.01528);
	// bethe->SetLineWidth(3);
    // bethe->Draw("SAME");
	
	// region indicator
    TBox *box = new TBox();box->SetFillStyle(4010);box->SetFillColorAlpha(kBlue - 8, 0.2);
    box->DrawBox(2.e-2, 1., 3.e-2, 2.e3);
    box->DrawBox(1.6, 1., 2.40, 2.e3);
    return 0;
}

