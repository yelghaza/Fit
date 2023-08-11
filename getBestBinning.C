#include <iostream>
#include <cmath>
//#ifndef  __ATLASSTYLE_H
//#include "AtlasStyle.C"
//#endif // __ATLASSTYLE_H

std::vector<std::string> inputFiles = {"./"};
std::string channel = "1Lep";

std::vector<std::string> BkgList; 


TF1 *linearfit= new TF1("linearfit","[0]+x*[1]",0,10); //fit width vs mass
void drawTex(const std::string& name, double x, double y){
	TLatex * tex4 = new TLatex(x,y,name.c_str() );
	tex4->SetNDC();
	tex4->SetTextFont(72);
	tex4->SetTextSize(0.03);
	tex4->SetLineWidth(2);
	tex4->Draw();

}
std::vector<std::string> getListOfSamples(const std::string& signalModel){
	std::vector<std::string> SignalList;
	for(int j=0;j<inputFiles.size();j++){
		TSystemDirectory histSearchDir((inputFiles[j]).c_str(),(inputFiles[j]).c_str());
		TList* files = histSearchDir.GetListOfFiles();
		if(files){
			TSystemFile *file;
			TIter next(files);
			while((file=(TSystemFile*)next())){
				string fname= file->GetName();
				if(fname.find(".root") == std::string::npos) continue;
				if(fname.find("Wjets") != std::string::npos) BkgList.push_back(inputFiles[j]+"/"+fname);
				if(fname.find("ttbar") != std::string::npos) BkgList.push_back(inputFiles[j]+"/"+fname);
				if(fname.find("stop") != std::string::npos) BkgList.push_back(inputFiles[j]+"/"+fname);
				if(fname.find("Diboson") != std::string::npos) BkgList.push_back(inputFiles[j]+"/"+fname);
				if(fname.find("Zjets") != std::string::npos) BkgList.push_back(inputFiles[j]+"/"+fname);
				if(fname.find("SMVH") != std::string::npos) BkgList.push_back(inputFiles[j]+"/"+fname);
				if(fname.find("ttH") != std::string::npos) BkgList.push_back(inputFiles[j]+"/"+fname);
				if(fname.find("ttV") != std::string::npos) BkgList.push_back(inputFiles[j]+"/"+fname);
				if(fname.find(signalModel) != std::string::npos) {
					
					if(signalModel.find("VBF") != std::string::npos){
						if(fname.find("VBF") != std::string::npos){
							SignalList.push_back(inputFiles[j]+"/"+fname);
						}
					}   
					else{
						if(fname.find("VBF") == std::string::npos) SignalList.push_back(inputFiles[j]+"/"+fname);
					}

				}
			}
			delete file;
		}
		delete files;
	}//for Loop
		
	return SignalList;
}
//get sum of the bkg
TH1F* getBkgHistSum(const std::vector<std::string>&  ListOfFiles, const std::string &regionName){
	std::string variable ="";
	if (channel =="0Lep") variable= regionName.find("Merg") != std::string::npos ? "_vvJ_m" : "_vvjj_m";
	if (channel =="1Lep") variable= regionName.find("Merg") != std::string::npos ? "_lvJ_m" : "_lvjj_m";
	if (channel =="2Lep") variable= regionName.find("Merg") != std::string::npos ? "_llJ_m" : "_lljj_m";

	std::string fullRegionName  = (regionName.find("Res") != std::string::npos) ? regionName + variable : regionName + variable;
	TH1F* hist_sum=NULL;;
	for (auto file : ListOfFiles){
				
		std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "READ"));

		if (!f) continue;
		TIter next(f->GetListOfKeys());
		TKey* key;
		//loop over keys
		while((key=(TKey*)next())){
			if(!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH1") ) continue;
			std::string name= key->GetName();

			if(name.find(fullRegionName)==std::string::npos) continue;
			//std::cout << name << std::endl;
			TH1F* hist = (TH1F*)key->ReadObj();
			if(hist_sum==NULL)
					hist_sum=(TH1F*)hist->Clone();
			else
					hist_sum->Add(hist);
			hist_sum->SetDirectory(0);
			f->Close();
		}//loop over keys
		delete key;
		//delete f;
	}// loop over files
	if (hist_sum){
		hist_sum->SetName(( "Bkg_" + regionName).c_str());
		//std::cout << "---------------------- Sum of Bkgs " << hist_sum ->Integral() << std::endl	
	}
	std::cout << "sum " << hist_sum->GetName() << " " << hist_sum ->Integral() << std::endl;
	return hist_sum;
}
//Get signal histograms 
std::map<int,TH1F* > getSigHistSum(const std::vector<std::string>& ListOfFiles, std::string sampleName, const std::string &regionName, const std::vector<int>& masspoints){
	std::string variable ="";
	if (channel =="0Lep") variable= regionName.find("Merg") != std::string::npos ? "_vvJ_m" : "_vvjj_m";
	if (channel =="1Lep") variable= regionName.find("Merg") != std::string::npos ? "_lvJ_m" : "_lvjj_m";
	if (channel =="2Lep") variable= regionName.find("Merg") != std::string::npos ? "_llJ_m" : "_lljj_m";


	std::map<int,TH1F* >signalHists;
	for (auto mpoint: masspoints){
		std::string masspoint = to_string(mpoint);
		std::string  histSigName = sampleName + "lvqq" + masspoint;
		TH1F* hist_sum=NULL;
		std::string fullRegionName = regionName.find("Merg") != std::string::npos ? ( histSigName + "_" + regionName+  variable ) : (histSigName + "_" + regionName + variable );
		std::cout << "fulleregionName " << fullRegionName << std::endl;
		for (auto file : ListOfFiles){
		

			std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "READ"));

			if (!f) continue;
			TIter next(f->GetListOfKeys());
			TKey* key;
			//loop over keys
			while((key=(TKey*)next())){
					if(!gROOT->GetClass(key->GetClassName())->InheritsFrom("TH1") ) continue;
					std::string name= key->GetName();
			
					if(name.find(fullRegionName)==std::string::npos) continue;
					std::cout << name << std::endl;
					TH1F* hist = (TH1F*)key->ReadObj();
					if(hist_sum==NULL)
							hist_sum=(TH1F*)hist->Clone();
					else
							hist_sum->Add(hist);
					hist_sum->SetDirectory(0);
					f->Close();
			}//loop over keys
			delete key;
			//delete f;
		}// loop over files 

		if (hist_sum){
				hist_sum->SetName((histSigName + "_" + sampleName + "_" + regionName).c_str());
				//std::cout << "---------------------- Sum of Bkgs " << hist_sum ->Integral() << std::endl;
				
		}
		//hist_sum->Scale(1./hist_sum->Integral());
		signalHists[mpoint] = hist_sum ;
	}
	return signalHists;
}
/// \param _x    The variable.
/// \param _Xp   The peak position.
/// \param _sigp The peak width as FWHM divided by 2*sqrt(2*log(2))=2.35
/// \param _xi   Peak asymmetry. Use values around 0.
/// \param _rho1 Left tail. Use slightly negative starting values.
/// \param _rho2 Right tail. Use slightly positive starting values.

double RooBukinPdf(double *X, double *par) //source: https://github.com/root-project/root/blob/master/roofit/roofit/src/RooBukinPdf.cxx
{   
		double x = X[0];
		double Xp = par[0]; // x_p
		double sigp = par[1]; // \sigma_p
		double xi = par[2]; // \xi
		double rho1 = par[3]; // \rho_1
		double rho2 = par[4]; // \rho_2
		double ap = par[5]; // \ \A_p //Normalisation factor

		const double consts = 2*sqrt(2*log(2.0));
		double r1=0,r2=0,r3=0,r4=0,r5=0,hp=0;
		double x1 = 0,x2 = 0;
		double fit_result = 0;

		hp=sigp*consts;
		r3=log(2.);
		r4=sqrt(xi*xi+1);
		r1=xi/r4;

		if(std::abs(xi) > exp(-6.)){
		r5=xi/log(r4+xi);
		}
		else
		r5=1;

		x1 = Xp + (hp / 2) * (r1-1);
		x2 = Xp + (hp / 2) * (r1+1);

		//--- Left Side
		if(x < x1){
		r2=rho1*(x-x1)*(x-x1)/(Xp-x1)/(Xp-x1)-r3 + 4 * r3 * (x-x1)/hp * r5 * r4/(r4-xi)/(r4-xi);
		}


		//--- Center
		else if(x < x2) {
		if(std::abs(xi) > exp(-6.)) {
				r2=log(1 + 4 * xi * r4 * (x-Xp)/hp)/log(1+2*xi*(xi-r4));
				r2=-r3*r2*r2;
		}
		else{
				r2=-4*r3*(x-Xp)*(x-Xp)/hp/hp;
		}
		}


		//--- Right Side
		else {
		r2=rho2*(x-x2)*(x-x2)/(Xp-x2)/(Xp-x2)-r3 - 4 * r3 * (x-x2)/hp * r5 * r4/(r4+xi)/(r4+xi);
		}

		if(std::abs(r2) > 100){
		fit_result = 0;
		}
		else{
		//---- Normalize the result
		fit_result = exp(r2);
		}

		return fit_result *ap;
}

// \param _x	The variable of the PDF.
// \param _x0	Location parameter of the Gaussian component.
// \param _sigmaL	Width parameter of the left side of the Gaussian component.
// \param _sigmaR	Width parameter of the right side of the Gaussian component.
// \param _alphaL	Location of transition to a power law on the left, in standard deviations away from the mean.
// \param _nL	Exponent of power-law tail on the left.
// \param _alphaR	Location of transition to a power law on the right, in standard deviations away from the mean.
// \param _\nR	Exponent of power-law tail on the right. 
Double_t DoubleSidedCB2(double *X, double *par) //https://root.cern/doc/master/classRooCrystalBall.html
{

		double x = X[0];
		double mu = par[0];
		double width = par[1];
		double a1 = par[2];
		double p1 = par[3];
		double a2 = par[4]; //alpha 
		double p2 = par[5]; // n

		double u   = (x-mu)/width;
		double A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
		double A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
		double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
		double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

		double result(1);
		if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
		else if (u<a2)  result *= TMath::Exp(-u*u/2);
		else            result *= A2*TMath::Power(B2+u,-p2);
		return result;
}


Double_t ExpGaussExp(double *X, double *par) {//from https://arxiv.org/pdf/1603.08591.pdf

		double x = X[0]; 
		double mean = par[0];
		double sigma = par[1];
		double kL = par[2];
		double kH = par[3];

		double u = (x - mean) / sigma;
		double result = 1.0;

		if (u <= -kL) {
				result *= exp((kL * kL / 2) + kL * u);
		} else if (u > -kL && u <= kH) {
				result *= exp(-u * u / 2);
		} else if (u > kH) {
				result *= exp(kH * u - (kH * kH / 2));
		}

		return result;
}
TF1* getFitFunction(int masspoint, double h_sigma){

		//////////////////////// BUKIN ///////////////////////////////////////////////////
		//TF1* f = new TF1("BukinPdf",RooBukinPdf, 0,1000,6);
		//f->SetParameters(masspoint,100,-0.2,-0.5,0,40); // 
		

		//////////////////// Double crystal ball ///////////////////////////////////////////
		TF1* f = new TF1("DoubleSidedCB",DoubleSidedCB2, 400.,700.,5);
		f->SetParameters(500, 10, 0.7, 1, 0.7,1); //  mean, sigma, aleft, pleft, aright, pright

		//////////////////////// Gaussian ///////////////////////////
		//TF1* f= new TF1("fit","[0]*TMath::Gaus(x,[1],[2])",0,3000);
		//f->SetParameters(1,masspoint,100);
		////////////////////////////////////////////////////////////
		//f->SetParLimits(1, 1., 10.);
		//f->SetParLimits(2, 1., 10.);
		//f->SetParLimits(3, 1., 10.);
		//f->SetParLimits(4, 1., 10.);
		//f->SetParLimits(5, 1., 10.);

		//////////////////////ExpGaussExp ////////////////////////////////////////////////
		//TF1* f = new TF1("ExpGaussExp", ExpGaussExp, 0., 3000., 7);
		//f->SetParameters(masspoint, 2, 1, 1); // mean, sigma KL, KH
		return f;
}
void getBestBin(std::string signal, std::string region){
	//SetAtlasStyle();
	std::vector<std::string> signalFiles = getListOfSamples(signal);

	
	std::vector<int> masses = {500, 600, 700, 800, 1000};
	std::map<int, TH1F*> sigHists = getSigHistSum(signalFiles, signal,region , masses);

	//for (auto sigHist: sigHists) std::cout << sigHist->GetName() << "  " << sigHist->Integral() << std::endl;


	TCanvas *c1 = new TCanvas( "c1", "A Simple Graph Example", 200, 10, 800, 600 );

	TGraph* WidthVsMass = new TGraph(sigHists.size());
	c1->Print("summaryFit.pdf[");
	for ( unsigned short i = 0;const auto &p : sigHists ){//for (unsigned short i = 0;auto sigHist: sigHists)	 
		int mass = p.first;
		TH1F *h_sig = p.second; //(TH1F*)sigHist->Clone();

		std::cout << mass << "  " << h_sig->Integral() << std::endl;
		
		std::cout << "Fitting signal " << signal << " in region " << region << std::endl; 
		std::cout << mass << " stdv " << h_sig->GetStdDev() << std::endl;
		TF1* fit = getFitFunction(mass, h_sig->GetStdDev());
		h_sig->Fit(fit, "E", "R", 0., 3000.);
		double mean  = fit->GetParameter(1);
		double sigma = fit->GetParameter(2);

		WidthVsMass->SetPoint(i, mass, sigma);

		h_sig->Draw("hist");
		fit->SetLineColor(kRed);
		fit->Draw("same");
		// show fir parameters 
		Float_t value = fit->GetParameter(i);
		

		gStyle->SetOptStat(0000);
		gStyle->SetOptFit(1111);

		c1->Print("summaryFit.pdf");

		delete h_sig;
		++i;
			
	}
	
	WidthVsMass->Draw("AP*");
	WidthVsMass->Fit(linearfit, "E");
	WidthVsMass->GetYaxis()->SetTitle("Width");
	WidthVsMass->GetXaxis()->SetTitle("m(X)[GeV]");
	c1->Print("summaryFit.pdf");
	
	
	c1->Print("summaryFit.pdf]");
	delete c1;
	delete WidthVsMass;
	sigHists.clear();
}


void getFitFunctionDis(){
		// plot the fit function
		TCanvas *c1 = new TCanvas( "c1", "A Simple Graph Example", 200, 10, 800, 600 );
		//TF1* f1 = new TF1("BukinPdf",RooBukinPdf, 0,3000.,7);
		//f1->SetParameters(500,100,-0.2,0.1,0,5);//SetParameters(masspoint, 10, 3, 12, 12,1); // 

        TF1* f1 = new TF1("DoubleSidedCB",DoubleSidedCB2, 0.,3000.,7);
		f1->SetParameters(500, 74, 10, 1, 10,1); //  mean, sigma, aleft, pleft, aright, pright

		f1->Draw();
		c1->Print("TESTFit.png");
}
void getBestBinning(){

		//getFitFunctionDis();
		getBestBin("HVTWZ", "VV1Lep_Res_GGF_WZ_01btag_SR");
		

	 	
}
