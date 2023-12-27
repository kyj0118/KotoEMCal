#include "TFile.h"
#include "TH1D.h"
#include "TRotation.h"
#include "TVector3.h"

void mkHist(double randpos = 0) {
  double pxpz;
  double pypz;
  double pxpz_true;
  double pypz_true;
  double x, y;
  double rad2deg = 180.0 / 3.14159265;
  for (string pos : {"position0", "position5mm"}) {
    for (string particle : {"electron", "positron"}) {
      for (string energy : {"200MeV_", "500MeV_", "1000MeV_"}) {
        string fname = "result/" + particle + energy + pos + ".txt";
        cout << fname << endl;
        TString foutname = "hist_" + particle + energy + pos + ".root";
        auto tf = new TFile(foutname, "RECREATE");
        ifstream ifile(fname);

        TH1D* hangle[9] = {NULL};
        TH2D* hangle2[9] = {NULL};
        for (int i = 0; i < 9; i++) {
          TString str_hname = Form("hangle_%d", i);
          TString str_hname2 = Form("hangle2_%d", i);

          hangle[i] = new TH1D(str_hname, ";#Delta #theta (deg)", 300, -15, 15);
          hangle2[i] = new TH2D(str_hname2, ";#phi (deg);#Delta #theta (deg)", 100, -180, 180, 150, -15, 15);
        }
        int ievt = 0;
        while (ifile >> pxpz && ifile >> pypz && ifile >> pxpz_true && ifile >> pypz_true) {
          ievt++;

          double theta = atan(sqrt(pxpz * pxpz + pypz * pypz));
          double phi = atan2(pypz, pxpz);
          double theta_true = atan(sqrt(pxpz_true * pxpz_true + pypz_true * pypz_true));
          double phi_true = atan2(pypz_true, pxpz_true);

          TVector3 RecVec(sin(theta) * cos(phi),
                          sin(theta) * sin(phi),
                          cos(theta));

          double val = 0;

          TVector3 GenVec(sin(theta_true) * cos(phi_true),
                          sin(theta_true) * sin(phi_true),
                          cos(theta_true));
          TVector3 unitZ(0, 0, 1);
          TRotation RotGenVecToZ;
          TVector3 v = GenVec.Cross(unitZ);
          RotGenVecToZ.Rotate(GenVec.Theta(), v);

          TVector3 RotToZVec = RotGenVecToZ * RecVec;
          TVector3 rotRecVec = RotGenVecToZ * RecVec;

          rotRecVec.RotateZ(-phi_true);

          if (!(pxpz_true == 0 && pypz_true == 0)) {
            val = asin(rotRecVec.Y());
          } else {
            val = asin(RecVec.Y());
          }
          val *= rad2deg;

          int i = (int)(theta_true * rad2deg + 0.1) / 5.0;
          // cout << i << ", " << theta << ", " << val << endl;
          hangle[i]->Fill(val);
          hangle2[i] -> Fill(phi_true*rad2deg,val);
        }
        ifile.close();
        for (int i = 0; i < 9; i++) {
          hangle[i]->Write();
          hangle2[i]->Write();
        }
        tf->Close();
      }
    }
  }
}
