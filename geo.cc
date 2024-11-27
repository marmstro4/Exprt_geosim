#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TEllipse.h>
#include <TLine.h>
#include <tuple>
#include <TGraph.h>
#include <vector>
#include <cmath>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <stdexcept>
#include "TVectorD.h"
#include <fstream>
#include <random>
#include <thread>
#include <chrono>
#include<TF1.h>

#include "Garfield/ComponentAnsys123.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewFEMesh.hh"
#include "Garfield/ViewMedium.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/AvalancheMicroscopic.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/Random.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/ViewDrift.hh"
#include "Garfield/ViewCell.hh"

using namespace Garfield;
using namespace std;

const double rCell = 0.5;

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>> GetStrawCenters() {

    std::vector<double> xCenter, yCenter, zCenter;

    double offset = 2.5;

    int nTubes = 10;

    for (int i = 0; i <= 2 * nTubes; ++i) {
                xCenter.push_back(offset);
                yCenter.push_back(-nTubes * 2 * rCell + i * 2 * rCell);
                zCenter.push_back(0);

    }

    /*

    // Calculate the number of tubes needed to cover +/-10 cm
    int nTubes = static_cast<int>(5 / strawRadius) + 1;
    //int nTubesH = static_cast<int>(2 / strawRadius) + 1;
    int nTubesH = rowN;

    for (int j = 0; j<=2*nTubesH; j++) {

        for (int i = 0; i <= 2 * nTubes; ++i) {

            if (j==0) {
                yCenter.push_back(verticalOffset);
                zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius);
            }
        }

        if (j==1) {
            for (int i = 0; i <= 2 * nTubes; ++i) {
                yCenter.push_back(verticalOffset + std::sqrt(3) * strawRadius);
                zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius +     strawRadius);
            }
        }

        if ((j!=0) && (j!=1) && (j % 2 == 0)) {
            for (int i = 0; i <= 2 * nTubes; ++i) {
                int k = (j-2);
                yCenter.push_back(verticalOffset + 2 * std::sqrt(3) * strawRadius+ 2*k*strawRadius);
                zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius);
            }
        }

        if ((j!=0) && (j!=1) && (j % 2 != 0)) {
            for (int i = 0; i <= 2 * nTubes; ++i) {
                int k = (j-2);
                yCenter.push_back(verticalOffset + 2 * std::sqrt(3) * strawRadius + 2*k*strawRadius);
                zCenter.push_back(-nTubes * 2 * strawRadius + i * 2 * strawRadius +
            strawRadius);
            }
        }

    }

    */

    return std::make_tuple(xCenter, yCenter, zCenter);

}

int main(int argc, char* argv[]) {

  TApplication app("app", &argc, argv);


  auto [xCenter,yCenter,zCenter] = GetStrawCenters();

  for (int i=0; i<xCenter.size(); i++) {
      cout<<xCenter[i]<<","<<yCenter[i]<<","<<zCenter[i]<<endl;
  }


  app.Run();
  return 0;
}
