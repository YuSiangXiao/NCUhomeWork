#include "/home/7202018/personalLib/Math/DBMMath.h"
#include "/home/7202018/personalLib/Math/DBMMath.h"
#include "/home/7202018/personalLib/RPU/DBMTH1F.h"
//constant
const double pi =     3.1415926575; //pi
const double G =      6.67430e-11;  //Gravity Constant
//intial condition(SI)
const double Me =     5.972e+26;    //Mass of Earth
const double Ms =     1.989E+30;    //Mass of Sun
const double t0 =     0;            //initial time
const double pex0 =   1.47e+11;     //position of earth on x
const double pey0 =   0 ;           //position of earth on y
const double psx0 =   0;            //position of sun on x
const double psy0 =   0;            //position of sun on y
//const double vex0 =   0;
const double vex0 =   0;            //velocity of earth on x
const double vey0 =   3.03e+4;      //velocity of earth on y
const double vsx0 =   0;            //velocity of sun on x
const double vsy0 =   0;            //velocity of sun on y
//control const
const int ntime =  73      ;        //point program execution
const double dt = 86400*5;          // dt
//variable, operator
double t = t0;                      //time variable
double dist = 0;                    //dist. variable
double dp=0,da=0,Ecc = 0;           //da,dp,Eccentricity variable
double Ve[2]= {vex0,vey0};          //vector of earth velocity
double Vs[2] = {vsx0,vsy0};         //vector of sun velocity
double Vei[2]= {vex0,vey0};         //Instantaneous Speed of earth velocity
double Vsi[2] = {vsx0,vsy0};        //Instantaneous Speed of sun velocity
double Pe[2] = {pex0,pey0};         //momentum of earth
double Ps[2] = {psx0,psy0};         //momentum of sun
double Ae[2], As[2], vecES[2];      //Acceleration of sun and earth and vector of two star
double ptt[ntime],pxte[ntime],pyte[ntime],pxts[ntime],pyts[ntime];
double Ee[ntime],Es[ntime],EkTot[ntime],EuTot[ntime],ETot[ntime];
double Vte[ntime],Vts[ntime];
//data array of time, position, energy, speed of two star on x-y, total energy, potential, kinetic
//fucnctions
double Dot(double *v1, double *v2){ //Dot
  //cout<<  v1[0]*v2[0]+v1[1]*v2[1]<<endl;
  return =v1[0]*v2[0]+v1[1]*v2[1];
}
//Gravity Calculate
void Gravity(double m1,double m2,double *p1,double *p2,double *A1,double *A2,double &Eu,double &a,double &p){
  vecES[0] = p2[0]-p1[0]; vecES[1] = p2[1]-p1[1]; //direction vector
  dist  = sqrt(Dot(vecES,vecES));                 //dist cal
  A1[0] = (G*m2)*(+vecES[0])/pow(dist,3);         //acc on x for earth
  A1[1] = (G*m2)*(+vecES[1])/pow(dist,3);         //acc on y for earth
  A2[0] = (G*m1)*(-vecES[0])/pow(dist,3);         //acc on x for sun
  A2[1] = (G*m1)*(-vecES[1])/pow(dist,3);         //acc on y for sun
  Eu    = -(G*m1*m2)/dist;                        //potential
  if (dist>a){
    a = dist;                                     //find max dist
  }else if(dist<p){
    p = dist;                                     //find min dist
  }
}
void Orbit(){ //main program
  dist = sqrt((pex0-psx0)*(pex0-psx0)+(pey0-psy0)*(pey0-psy0));//initial dist
  da = dist; dp = dist;                                        //initial da,dp
  for (int j=0;j<ntime;j++){                                   //for loop for track
    Gravity(Me,Ms,Pe,Ps,Ae,As,EuTot[j],da,dp);                 //calculate the gravity
    //cout<<EuTot[j]<<endl;
    for(int i=0;i<2;i++){                                      //calculate the track
      Ve[i] = Ve[i] + Ae[i]*dt;
      Vs[i] = Vs[i] + As[i]*dt;
      Vei[i] = Ve[i] + Ae[i];
      Vsi[i] = Vs[i] + As[i];
      Pe[i] = Pe[i] + Ve[i]*dt;
      Ps[i] = Ps[i] + Vs[i]*dt;
    }
                                                               //record data
    pxte[j]=Pe[0];pyte[j]=Pe[1];pxts[j]=Ps[0];pyts[j]=Ps[1];
    ptt[j]=t;
    Vte[j]= sqrt(Dot(Vei,Vei));
    Vts[j]= sqrt(Dot(Vsi,Vsi));
    Ee[j] = 0.5*Me*Vte[j]*Vte[j];
    Es[j] = 0.5*Ms*Vts[j]*Vts[j];
    EkTot[j] = Ee[j]+Es[j];
    ETot[j] = EkTot[j]+EuTot[j];
    t+=dt;  
    //cout<<t<<"\t"<<Ae[0]<<"\t"<<Ae[1]<<endl;
  }
  Ecc = (da-dp)/(da+dp);                                       //calculate Ecc
  
  TCanvas *c1 = new TCanvas("c1","Orbit_Track",810,810);
  {//drawing seting
  TGraph *OrbitE = new TGraph(ntime,pxte,pyte);
  TGraph *OrbitS = new TGraph(ntime,pxts,pyts);
  c1->cd();
  OrbitE->SetTitle("Sun-Earth Binary star system orbit prediction(Eular)");
  OrbitE->GetXaxis()->SetTitle("X(m)");
  OrbitE->GetYaxis()->SetTitle("Y(m)");
  OrbitE->GetYaxis()->SetTitleOffset(1.3);
  OrbitS->SetTitle("");

  OrbitE->SetMarkerStyle(4);
  OrbitE->SetMarkerColor(1);
  OrbitS->SetMarkerStyle(3);
  OrbitS->SetMarkerColor(1);
  OrbitE->Draw("ap");
  //OrbitE->GetXaxis()->SetRangeUser(-1.0E+9, +1.0E+9);
  //OrbitE->GetYaxis()->SetRangeUser(-1.0E+9, +1.0E+9);
  OrbitS->Draw("samep");
  TLegend *LOrbitS = new TLegend(0.7,0.8,0.89,0.89);
  LOrbitS->AddEntry(OrbitE,"Earth trace","p");
  LOrbitS->AddEntry(OrbitE,Form("e =%.4f",Ecc),"");
  LOrbitS->AddEntry(OrbitS,"Sun trace","p");
  
  setLegendDefault(LOrbitS);
  LOrbitS->Draw();
  TLegend *Ldp = new TLegend(0.435,0.76,0.645,0.79);
  Ldp->AddEntry(OrbitE,Form("d_{p}=%.3e",dp),"");
  setLegendDefault(Ldp);
  Ldp->Draw();
  TLegend *Lda = new TLegend(0.435,0.20,0.645,0.23);
  Lda->AddEntry(OrbitE,Form("d_{a}=%.3e",da),"");
  setLegendDefault(Lda);
  Lda->Draw();
  }
  TCanvas *c2 = new TCanvas("c2","Orbit_PvsT",810,810);
  {//drawing seting
  c2->Divide(1,2);
  TGraph *XTE = new TGraph(ntime,ptt,pxte);
  TGraph *XTS = new TGraph(ntime,ptt,pxts);
  TGraph *YTE = new TGraph(ntime,ptt,pyte);
  TGraph *YTS = new TGraph(ntime,ptt,pyts);
  c2->cd(1);
  XTE->SetMarkerStyle(4);
  XTE->SetMarkerColor(1);
  XTS->SetMarkerStyle(3);
  XTS->SetMarkerColor(1);
  XTE->Draw("AP");
  XTS->Draw("p");
  XTE->SetTitle("Sun-Earth Binary star system X-T(Eular)");
  XTE->GetXaxis()->SetTitle("T(s)");
  XTE->GetYaxis()->SetTitle("X(m)");
  TLegend *LXTE = new TLegend(0.55,0.8,0.74,0.89);
  LXTE->AddEntry(OrbitE,"Earth trace","p");
  LXTE->AddEntry(OrbitS,"Sun trace","p");
  setLegendDefault(LXTE);
  LXTE->Draw();
  c2->cd(2);
  YTE->SetMarkerStyle(4);
  YTE->SetMarkerColor(1);
  YTS->SetMarkerStyle(3);
  YTS->SetMarkerColor(1);
  YTE->Draw("AP");
  YTS->Draw("p");
  YTE->SetTitle("Sun-Earth Binary star system Y-T(Eular)");
  YTE->GetXaxis()->SetTitle("T(s)");
  YTE->GetYaxis()->SetTitle("Y(m)");
  TLegend *LYTE = new TLegend(0.55,0.8,0.74,0.89);
  LYTE->AddEntry(OrbitE,"Earth trace","p");
  LYTE->AddEntry(OrbitS,"Sun trace","p");
  setLegendDefault(LYTE);
  LYTE->Draw();
  }
  TCanvas *c3 = new TCanvas("c3","Orbit_EkvsT",810,810);
  {//drawing seting
  c3->Divide(2,2);
  TGraph *VTE = new TGraph(ntime,ptt,Vte);
  TGraph *VTS = new TGraph(ntime,ptt,Vts);
  VTE->SetMarkerStyle(4);
  VTE->SetMarkerColor(1);
  VTS->SetMarkerStyle(3);
  VTS->SetMarkerColor(1);
  c3->cd(1);
  VTE->Draw("pa");
  VTE->SetTitle("Sun-Earth Binary star system |V_{e}|-T(Eular)");
  VTE->GetXaxis()->SetTitle("T(s)");
  VTE->GetYaxis()->SetTitle("|v|(m/s)");
  VTE->GetYaxis()->SetTitleOffset(1.3);
  c3->cd(2);
  VTS->Draw("pa");
  VTS->SetTitle("Sun-Earth Binary star system |v_{s}|-T(Eular)");
  VTS->GetXaxis()->SetTitle("T(s)");
  VTS->GetYaxis()->SetTitle("|v|(m/s)");
  VTS->GetYaxis()->SetTitleOffset(1.3);
  
  TLegend *LVTS = new TLegend(0.7,0.8,0.89,0.89);
  LVTS->AddEntry(OrbitE,"Earth speed","p");
  LVTS->AddEntry(OrbitS,"Sun speed","p");
  setLegendDefault(LVTS);
  LVTS->Draw();
  
  TGraph *ETE = new TGraph(ntime,ptt,Ee);
  TGraph *ETS = new TGraph(ntime,ptt,Es);
  ETE->SetMarkerStyle(4);
  ETE->SetMarkerColor(1);
  ETS->SetMarkerStyle(3);
  ETS->SetMarkerColor(1);
  c3->cd(3);
  ETE->Draw("pa");
  ETE->SetTitle("Sun-Earth Binary star system Ek_{e}-T(Eular)");
  ETE->GetXaxis()->SetTitle("T(s)");
  ETE->GetYaxis()->SetTitle("Ek(J)");
  ETE->GetYaxis()->SetTitleOffset(1.3);
  
  c3->cd(4);
  ETS->Draw("pa");
  ETS->SetTitle("Sun-Earth Binary star system Ek_{s}-T(Eular)");
  ETS->GetXaxis()->SetTitle("T(s)");
  ETS->GetYaxis()->SetTitle("Ek(J)");
  ETS->GetYaxis()->SetTitleOffset(1.3);
  TLegend *LETS = new TLegend(0.7,0.8,0.89,0.89);
  LETS->AddEntry(OrbitE,"Earth energy","p");
  LETS->AddEntry(OrbitS,"Sun energy","p");
  setLegendDefault(LETS);
  LETS->Draw();
  }
  TCanvas *c4 = new TCanvas("c4","Orbit_TotalEvsT",810,1000);
  {//drawing seting
  c4->Divide(1,3);
  TGraph *TTotEk = new TGraph(ntime,ptt,EkTot);
  TGraph *TTotEu = new TGraph(ntime,ptt,EuTot);
  TGraph *TTotE  = new TGraph(ntime,ptt,ETot);
  TTotEk->SetMarkerStyle(4);
  TTotEk->SetMarkerColor(1);
  TTotEu->SetMarkerStyle(3);
  TTotEu->SetMarkerColor(1);
  TTotE ->SetMarkerStyle(5);
  TTotE ->SetMarkerColor(1);
  TTotEk->SetTitle("Sun-Earth Binary star system total Ek-T(Eular)");
  TTotEu->SetTitle("Sun-Earth Binary star system total Eu-T(Eular)");
  TTotE ->SetTitle("Sun-Earth Binary star system total E-T(Eular)");
  TTotEk->GetXaxis()->SetTitle("T(s)");
  TTotEk->GetYaxis()->SetTitle("Ek(J)");
  TTotEu->GetXaxis()->SetTitle("T(s)");
  TTotEu->GetYaxis()->SetTitle("Eu(J)");
  TTotE->GetXaxis()->SetTitle("T(s)");
  TTotE->GetYaxis()->SetTitle("E(J)");
  c4->cd(1);
  TTotEk->Draw("Ap");
  TLegend *LTTotEk = new TLegend(0.7,0.91,0.89,1);
  LTTotEk->AddEntry(TTotEk,"Kinetic energy","p");
  setLegendDefault(LTTotEk);
  LTTotEk->Draw();
  c4->cd(2);
  TTotEu->Draw("Ap");
  TLegend *LTTotEu = new TLegend(0.7,0.91,0.89,1);
  LTTotEu->AddEntry(TTotEu,"Potential energy","p");
  setLegendDefault(LTTotEu);
  LTTotEu->Draw();
  c4->cd(3);
  TTotE->Draw("Ap");
  TLegend *LTTotE = new TLegend(0.7,0.91,0.89,1);
  LTTotE->AddEntry(TTotE,"Total energy","p");
  setLegendDefault(LTTotE);
  LTTotE->Draw();
  }
  c1->Print("Orbit_Track.pdf");     //show the track and da,dp,Ecc of two star
  c2->Print("Orbit_PvsT.pdf");      //show the momentum Vs Time of two star
  c3->Print("Orbit_EkvsT.pdf");     //show the Kinetic energy and speed of two star
  c4->Print("Orbit_TotalEvsT.pdf"); //show the Mechanical energy of all system
}



//if ((HLTEleMuX >> 41 & 1) == 0 && (HLTEleMuX >> 42 & 1) == 0) continue; 
