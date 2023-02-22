#define plotone_cxx
#include "plotone.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TVirtualFFT.h"

void plotone::LoopFFT(int ichoice, int chchoice)
{

  if (fChain == 0) return;
  Long64_t jentry = ichoice;

  TCanvas *c1=new TCanvas("c1",Form("EventFFT%d.pdf",(int)jentry),1500,900);
  c1->Divide(4,4); // divides the canvas into two rows and three columns
  
  std::cout << " nsamples per waveform setting " << nsamples << " " << fRun << std::endl;

  TH1F *wf[16];
  for (int i=0;i<16;++i) {
    TString title = Form("ch%02d",i);
    wf[i] = new TH1F(title,title,nsamples,0,nsamples);
    wf[i]->GetXaxis()->SetTitle("time (1 tick = 2 ns)");
  }

  Long64_t nbytes = 0, nb = 0;

  Long64_t ientry = LoadTree(jentry);
  if (ientry < 0) return;
  nb = fChain->GetEntry(jentry);   nbytes += nb;
 
  // make pictures of waveforms
  for (auto ich = 0;ich<fWvfmsVec->size();++ich) {
    vector<unsigned short> thiswf = fWvfmsVec->at(ich);
    if (jentry==0)	{ if (nsamples!=thiswf.size()) std::cout << "Number of samples set at " << nsamples <<
              " and waveforms are " << thiswf.size() << std::endl;}
    for (auto itick=0;itick<thiswf.size();++itick) {
      wf[ich]->SetBinContent(itick+1,thiswf[itick]);
    }
  }

  /*int wvsize=fWvfmsVec->at(chchoice).size();
  int n=wvsize+1;
  Double_t *in = new Double_t[2*((n+1)/2+1)];
  for(int i=0; i<wvsize; i++ ){
    in[i] = fWvfmsVec->at(chchoice).at(i);
    std::cout << in[i] << std::endl;
  }



  wvsize=2*wvsize+1;
  TVirtualFFT *fft_own = TVirtualFFT::FFT(1, &wvsize, "R2C ES K");
  fft_own->SetPoints(in);
  fft_own->Transform();
  //Copy all the output points:
  fft_own->GetPoints(in);

  TH1 *hr = 0;
  hr = TH1::TransformHisto(fft_own, hr, "RE");*/

  c1->cd(1);
  wf[0]->Draw();

  // TH1F *wfFFT[16]; for (auto ich = 0;ich<wf->size();++ich) {}
  
  c1->cd(2);
  TH1 *hm =0;
  TVirtualFFT::SetTransform(0);
  hm = wf[0]->FFT(hm, "MAG");
  hm->SetTitle("Magnitude of the 1st transform");
  hm->Draw();
 
  c1->cd();c1->Update(); c1->WaitPrimitive();
  //c1->Print(Form("Event%d.pdf",(int)jentry));

}

void plotone::Loop(int ichoice)
{
//   In a ROOT session, you can do:
//      root> .L plotone.C
//      root> plotone t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
    Long64_t jentry = ichoice;

   TCanvas *c1=new TCanvas("c1",Form("Event%d.pdf",(int)jentry),1500,900);
   c1->Divide(4,4); // divides the canvas into two rows and three columns


   //int nsamples = 5500;//15000;
   /* if (fRun>2350 ) nsamples=15000; */
   std::cout << " nsamples per waveform setting " << nsamples << " " << fRun << std::endl;

   TH1F *wf[16];
   for (int i=0;i<16;++i) {
     TString title = Form("ch%02d",i);
     wf[i] = new TH1F(title,title,nsamples,0,nsamples);
     wf[i]->GetXaxis()->SetTitle("time (1 tick = 2 ns)");
   }

   /* Long64_t nentries = fChain->GetEntriesFast(); */
   /* std::cout << nentries << std::endl; */

    Long64_t nbytes = 0, nb = 0;
   /* for (Long64_t jentry=0; jentry<nentries;jentry++) { */

   Long64_t ientry = LoadTree(jentry);
   if (ientry < 0) return;
   nb = fChain->GetEntry(jentry);   nbytes += nb;
   // if (Cut(ientry) < 0) continue;

   // make pictures of waveforms
   //      std::cout << "Number of waveforms" << fWvfmsVec->size() << std::endl;
   for (auto ich = 0;ich<fWvfmsVec->size();++ich) {
     vector<unsigned short> thiswf = fWvfmsVec->at(ich);
     if (jentry==0)	{ if (nsamples!=thiswf.size()) std::cout << "Number of samples set at " << nsamples <<
							 " and waveforms are " << thiswf.size() << std::endl;}
     //	std::cout << "Number of samples " << thiswf.size() << std::endl;
     for (auto itick=0;itick<thiswf.size();++itick) {
       wf[ich]->SetBinContent(itick+1,thiswf[itick]);
     }
     // determine if channel has signal

     //  determine turn on time
   }

   //   wf[0]->SetTitle("counter 1");
   //wf[0]->SetTitle("channel 0");
//   wf[1]->SetTitle("Trigger");//channel 1");
   //wf[2]->SetTitle("counter 1");
   //wf[2]->SetTitle("channel 2");
//   wf[3]->SetTitle("RWM");//channel 3");
   //wf[4]->SetTitle("counter 2");
   //wf[4]->SetTitle("channel 4");
//   wf[5]->SetTitle("RWM");//channel 5");
   //wf[6]->SetTitle("channel 6");
//   wf[7]->SetTitle("RWM");//channel 7");
   //wf[8]->SetTitle("counter 3");
   //wf[8]->SetTitle("channel 8");
//   wf[9]->SetTitle("RWM");//channel 9");
   //wf[10]->SetTitle("counter 4");
   //wf[10]->SetTitle("channel 10");
//   wf[11]->SetTitle("RWM");//channel 11");
   //wf[12]->SetTitle("counter 5");
   //wf[12]->SetTitle("channel 12");
//   wf[13]->SetTitle("RWM");//channel 13");
   //wf[14]->SetTitle("channel 14");
//   wf[15]->SetTitle("Trigger");//channel 15");
   //wf[7]->SetTitle("event trigger");
   //wf[14]->SetTitle("light trigger");
   //wf[15]->SetTitle("RWM");

   c1->cd(1);  wf[0]->Draw();
   c1->cd(1);  wf[0]->Draw();
   c1->cd(2);  wf[1]->Draw();
   c1->cd(3);  wf[2]->Draw();
   c1->cd(4);  wf[3]->Draw();
   c1->cd(5);  wf[4]->Draw();
   c1->cd(6);  wf[5]->Draw();
   c1->cd(7);  wf[6]->Draw();
   c1->cd(8);  wf[7]->Draw();
   c1->cd(9);  wf[8]->Draw();
   c1->cd(10);  wf[9]->Draw();
   c1->cd(11);  wf[10]->Draw();
   c1->cd(12);  wf[11]->Draw();
   c1->cd(13);  wf[12]->Draw();
   c1->cd(14);  wf[13]->Draw();
   c1->cd(15);  wf[14]->Draw();
   c1->cd(16);  wf[15]->Draw();
   /*c1->cd(1);  wf[2]->Draw();
   c1->cd(2);  wf[4]->Draw();
   c1->cd(3);  wf[8]->Draw();
   c1->cd(4);  wf[10]->Draw();
   c1->cd(5);  wf[12]->Draw();
   c1->cd(6);  wf[7]->Draw();
   c1->cd(7);  wf[14]->Draw();
   c1->cd(8);  wf[15]->Draw();*/

   c1->cd();c1->Update(); c1->WaitPrimitive();

    //c1->Print(Form("../../6Desktop/PlotsSBND/Event%02d.png",(int)jentry));
    c1->Print(Form("Event%d.pdf",(int)jentry));

}

void plotone::Average()
{
//   In a ROOT session, you can do:
//      root> .L plotone.C
//      root> plotone t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Average();       // Loop on all entries and get average
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TCanvas *c1=new TCanvas("c1","c1",1500,900);
   c1->Divide(1,2); // divides the canvas into two rows and three columns


   //int nsamples = 15000;
   /* if (fRun>2350 ) nsamples=15000; */
   std::cout << " nsamples per waveform setting " << nsamples << " " << fRun << std::endl;

   TH1F *wf[16];
   for (int i=0;i<16;++i) {
     TString title = Form("ch%02d",i);
     wf[i] = new TH1F(title,title,280,0,280);//nsamples,0,nsamples);
     wf[i]->GetXaxis()->SetTitle("time around peak (1 tick=2 ns)");
     wf[i]->SetBit(TH1::kIsAverage);
   }

   TH1F *wf_temp[16];
   for (int i=0;i<16;++i) {
     TString title_temp = Form("ch%02d",i);
     wf_temp[i] = new TH1F(title_temp,title_temp,280,0,280);//nsamples,0,nsamples);
     wf_temp[i]->GetXaxis()->SetTitle("time around peak (1 tick=2 ns)");
     wf_temp[i]->SetBit(TH1::kIsAverage);
   }

    Long64_t nentries = fChain->GetEntriesFast();
    std::cout << nentries << std::endl;

    Long64_t nbytes = 0, nb = 0;
 for (Long64_t jentry=0; jentry<nentries;jentry++) {
   //Long64_t jentry = ichoice;
   Long64_t ientry = LoadTree(jentry);
   if (ientry < 0) return;
   nb = fChain->GetEntry(jentry);   nbytes += nb;
   // if (Cut(ientry) < 0) continue;


   // make pictures of waveforms
   //      std::cout << "Number of waveforms" << fWvfmsVec->size() << std::endl;
   for (auto ich = 0;ich<fWvfmsVec->size();++ich) {
     vector<unsigned short> thiswf = fWvfmsVec->at(ich);
     if (jentry==0)	{ if (nsamples!=thiswf.size()) std::cout << "Number of samples set at " << nsamples <<
							 " and waveforms are " << thiswf.size() << std::endl;}
     //	std::cout << "Number of samples " << thiswf.size() << std::endl;
     auto starttick=0;
     for (auto itick=0;itick<thiswf.size();++itick) {
       if (thiswf[itick]>4000) {
         starttick=itick-100;
         break;
       }
     }
     if (starttick<0)starttick=0;
     auto endtick = starttick+280;
     if (endtick>thiswf.size()-1)endtick=thiswf.size()-1;
     auto bintick = 1;
     for (auto itick=starttick;itick<=endtick;++itick) {
       wf_temp[ich]->SetBinContent(bintick,thiswf[itick]);
       bintick++;
     }
     // determine if channel has signal
     wf[ich]->Add(wf[ich],wf_temp[ich]);
     wf_temp[ich]->Reset();
     //  determine turn on time
   }
 }

   //   wf[0]->SetTitle("counter 1");
   wf[0]->SetTitle("channel 0");
   wf[1]->SetTitle("channel 1");
   wf[2]->SetTitle("counter 1");
   wf[3]->SetTitle("channel 3");
   wf[4]->SetTitle("counter 2");
   wf[5]->SetTitle("channel 5");
   wf[6]->SetTitle("channel 6");
   //wf[7]->SetTitle("channel 7");
   wf[8]->SetTitle("counter 3");
   wf[9]->SetTitle("channel 9");
   wf[10]->SetTitle("counter 4");
   wf[11]->SetTitle("channel 11");
   wf[12]->SetTitle("counter 5");
   wf[13]->SetTitle("channel 13");
   wf[7]->SetTitle("event trigger");
   wf[14]->SetTitle("light trigger");
   wf[15]->SetTitle("RWM");

   //   c1->cd(1);  wf[0]->Draw();
   c1->cd(1);  wf[0]->Draw();
   c1->cd(2);  wf[1]->Draw();
   //c1->cd(3);  wf[2]->Draw();
   //c1->cd(4);  wf[3]->Draw();
   //c1->cd(5);  wf[4]->Draw();
   //c1->cd(6);  wf[5]->Draw();
   //c1->cd(7);  wf[6]->Draw();
   //c1->cd(8);  wf[7]->Draw();
   //c1->cd(9);  wf[8]->Draw();
   //c1->cd(10);  wf[9]->Draw();
   //c1->cd(11);  wf[10]->Draw();
   //c1->cd(12);  wf[11]->Draw();
   //c1->cd(13);  wf[12]->Draw();
   //c1->cd(14);  wf[13]->Draw();
   //c1->cd(15);  wf[14]->Draw();
   //c1->cd(16);  wf[15]->Draw();
   /*c1->cd(1);  wf[2]->Draw();
   c1->cd(2);  wf[4]->Draw();
   c1->cd(3);  wf[8]->Draw();
   c1->cd(4);  wf[10]->Draw();
   c1->cd(5);  wf[12]->Draw();
   c1->cd(6);  wf[7]->Draw();
   c1->cd(7);  wf[14]->Draw();
   c1->cd(8);  wf[15]->Draw();*/

      //      c1->Print(Form("../../6Desktop/PlotsSBND/Event%02d.png",(int)jentry));
      //c1->Print(Form("Event%d.png",(int)jentry));

}
