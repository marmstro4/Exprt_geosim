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
#include <TBox.h>
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
double offset = 2.5;
double length = 12;
int layers = 10;
int rows = static_cast<int>(length / (2*rCell));

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>> GetStrawCentersLongitudinal() {

    std::vector<double> xCenter, yCenter, zCenter;

    double offset = 2.5;

    double Vspace = 2.5;
    double Hspace = 10;

    int nTubesV = static_cast<int>(0.5*Vspace / rCell);
    int nTubesH = static_cast<int>(0.5*Hspace / rCell);

    //Original
    for (int j = 0; j<nTubesH; ++j) {
        for (int i = -j-1; i <= 2*nTubesV+j; ++i) {
                xCenter.push_back(offset+j*2*rCell);
                yCenter.push_back(-nTubesV*2*rCell+i*2*rCell+rCell);
                zCenter.push_back(0);
        }
    }

    //Rotate 90
    std::vector<double> x_rotated,y_rotated,z_rotated;

    for (int j = 0; j<nTubesH; ++j) {
        for (int i = -j; i <= 2*nTubesV+j-1; ++i) {
                y_rotated.push_back(-(offset+j*2*rCell));
                x_rotated.push_back(-nTubesV*2*rCell+i*2*rCell+rCell);
                z_rotated.push_back(0);
        }
    }

    //Reflect verticals
    std::vector<double> x_flip,y_flip,z_flip;

    for (size_t i = 0; i < xCenter.size(); ++i) {
        x_flip.push_back(-xCenter[i]);
        y_flip.push_back(yCenter[i]);
        z_flip.push_back(zCenter[i]);  // z remains unchanged
    }

      // Append the rotated coordinates to the original vectors
    xCenter.insert(xCenter.end(), x_flip.begin(), x_flip.end());
    yCenter.insert(yCenter.end(), y_flip.begin(), y_flip.end());
    zCenter.insert(zCenter.end(), z_flip.begin(), z_flip.end());

    //Add rotated
    for (size_t i = 0; i < x_rotated.size(); ++i) {
        xCenter.push_back(x_rotated[i]);
        yCenter.push_back(y_rotated[i]);
        zCenter.push_back(z_rotated[i]);  // z remains unchanged
    }


    for (size_t i = 0; i < x_rotated.size(); ++i) {
        xCenter.push_back(x_rotated[i]);
        yCenter.push_back(-y_rotated[i]);
        zCenter.push_back(z_rotated[i]);  // z remains unchanged
    }

    return std::make_tuple(xCenter, yCenter, zCenter);

}

std::tuple<std::vector<double>, std::vector<double>,std::vector<double>> GetStrawCentersTransverse() {

    std::vector<double> xCenter, yCenter, zCenter;

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(-2*offset);
            yCenter.push_back(offset+j*2*rCell+rCell);
            zCenter.push_back(i*2*rCell+rCell);
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(offset+j*2*rCell+rCell);
            yCenter.push_back(2*offset);
            zCenter.push_back(i*2*rCell+rCell);
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(2*offset);
            yCenter.push_back(-(offset+j*2*rCell+rCell));
            zCenter.push_back(i*2*rCell+rCell);
        }
    }

    for (int i = -rows; i<rows; ++i) {
        for (int j = 0; j<layers; ++j) {
            xCenter.push_back(-(offset+j*2*rCell+rCell));
            yCenter.push_back(-2*offset);
            zCenter.push_back(i*2*rCell+rCell);
        }
    }

    return std::make_tuple(xCenter, yCenter, zCenter);

}

void PlotMidYZ(const std::vector<double>& yCenter, const std::vector<double>& zCenter) {
    TCanvas* canvas = new TCanvas("canvas", "Straws and Track", 800, 800);
    // Create a canvas
    canvas->SetFixedAspectRatio(); // Ensure equal scaling on both axes
    canvas->DrawFrame(-20, -20, 20, 20); // Set the drawing frame (xMin, yMin, xMax, yMax)

    TGraph* graph = new TGraph(yCenter.size(), zCenter.data(), yCenter.data());
    graph->SetLineStyle(0);
    graph->SetMarkerStyle(20); // Set marker style (20 is a solid circle)
    graph->SetMarkerSize(1.0); // Default size, can be adjusted
    graph->SetMarkerColor(kRed);  // Red points
    graph->Draw("AP");

    canvas->Update();
}

void PlotXYCells(const std::vector<double>& xCenter, const std::vector<double>& yCenter) {
    // Check if the inputs are valid

    TCanvas* canvas = new TCanvas("canvas", "Straws and Track", 800, 800);
    // Create a canvas
    canvas->SetFixedAspectRatio(); // Ensure equal scaling on both axes
    canvas->DrawFrame(-20, -20, 20, 20); // Set the drawing frame (xMin, yMin, xMax, yMax)


    // Draw the box (straws)
    for (int i = 0; i < rows*2*layers; ++i) {
        TBox* straw = new TBox(xCenter[i]-3*offset, yCenter[i]-rCell, xCenter[i]+3*offset, yCenter[i]+rCell);
        straw->SetFillStyle(0); // No fill
        straw->SetLineColor(kBlack);
        straw->Draw("same");
    }

    for (int i = rows*2*layers; i < 2*rows*2*layers; ++i) {
        TBox* straw1 = new TBox(xCenter[i]-rCell, yCenter[i]-3*offset, xCenter[i]+rCell, yCenter[i]+3*offset);
        straw1->SetFillStyle(0); // No fill
        straw1->SetLineColor(kRed);
        straw1->Draw("same");
    }

    for (int i = 2*rows*2*layers; i < 2*rows*3*layers; ++i) {
        TBox* straw2 = new TBox(xCenter[i]-3*offset, yCenter[i]-rCell, xCenter[i]+3*offset, yCenter[i]+rCell);
        straw2->SetFillStyle(0); // No fill
        straw2->SetLineColor(kBlue);
        straw2->Draw("same");
    }

    for (int i = 2*rows*3*layers; i < 2*rows*4*layers; ++i) {
        TBox* straw3 = new TBox(xCenter[i]-rCell, yCenter[i]-3*offset, xCenter[i]+rCell, yCenter[i]+3*offset);
        straw3->SetFillStyle(0); // No fill
        straw3->SetLineColor(kGreen);
        straw3->Draw("same");

    }

    // Update the canvas to render the drawing
    canvas->Update();
}

void PlotTrack(std::vector<double> trk_x,std::vector<double> trk_y) {

  TGraph* graph = new TGraph(trk_x.size(), trk_x.data(), trk_y.data());

  // Set the marker style for the points (e.g., a circle)
  graph->SetMarkerStyle(21);  // 21: circle
  graph->SetMarkerColor(kRed);  // Red points
  graph->SetMarkerSize(3);  // Marker size

  // Draw the graph on the canvas
  graph->Draw("same");

}

std::vector<double> generateDir() {
    // Seed the random number generator
    srand(static_cast<unsigned>(time(nullptr)));

    static std::random_device rd; // Random device for seeding
    static std::mt19937 gen(rd()); // Mersenne Twister generator
    static std::uniform_real_distribution<> dis(0,2);

    // Generate uniform random angles
    double phi = M_PI*dis(gen);  // Azimuthal angle (0 to 2π)
    double theta = acos(1.0 - dis(gen));  // Polar angle (0 to π)

    // Convert spherical coordinates to Cartesian coordinates
    double x = sin(theta) * cos(phi);
    double y = sin(theta) * sin(phi);
    double z = cos(theta);

    return {x,y,z};

}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>,
           std::vector<double>, std::vector<double>, std::vector<double>>
hits(std::vector<double> dirvect, std::vector<double> xcell, std::vector<double> ycell, std::vector<double> zcell) {

    std::vector<double> xcells, ycells, zcells;
    std::vector<double> xhits, yhits, zhits;

    double L = 12.5;
    double ux = dirvect[0], uy = dirvect[1], uz = dirvect[2];

    for (size_t i = 0; i < zcell.size(); i++) {

        double x0 = xcell[i], y0 = ycell[i], z0 = zcell[i];

        double A = uy * uy + uz * uz;
        double B = 2 * (y0 * uy + z0 * uz);
        double C = y0 * y0 + z0 * z0 - rCell * rCell;

        // Calculate the discriminant
        double discriminant = B * B - 4 * A * C;

        if (discriminant < 0) {
            //std::cout << "No intersection" << std::endl;
            continue;
        }

        // Calculate the two solutions for t
        double t1 = (-B - std::sqrt(discriminant)) / (2 * A);
        double t2 = (-B + std::sqrt(discriminant)) / (2 * A);

        // Compute intersection points
        double x1 = x0 + t1 * ux;
        double y1 = y0 + t1 * uy;
        double z1 = z0 + t1 * uz;

        double x2 = x0 + t2 * ux;
        double y2 = y0 + t2 * uy;
        double z2 = z0 + t2 * uz;

        double diff1 = std::sqrt(x1 * x1 + y1 * y1 + z1 * z1);
        double diff2 = std::sqrt(x2 * x2 + y2 * y2 + z2 * z2);

        if (x1 >= x0 - L / 2 && x1 <= x0 + L / 2) {
            if (diff2 > diff1) {
                xhits.push_back(x1);
                yhits.push_back(y1);
                zhits.push_back(z1);

                xcells.push_back(x0);
                ycells.push_back(y0);
                zcells.push_back(z0);
            } else if (x2 >= x0 - L / 2 && x2 <= x0 + L / 2) {
                xhits.push_back(x2);
                yhits.push_back(y2);
                zhits.push_back(z2);

                xcells.push_back(x0);
                ycells.push_back(y0);
                zcells.push_back(z0);
            }
        }

        if (x2 >= x0 - L / 2 && x2 <= x0 + L / 2) {
            if (diff2 < diff1) {
                xhits.push_back(x2);
                yhits.push_back(y2);
                zhits.push_back(z2);

                xcells.push_back(x0);
                ycells.push_back(y0);
                zcells.push_back(z0);
            } else if (x1 >= x0 - L / 2 && x1 <= x0 + L / 2) {
                xhits.push_back(x1);
                yhits.push_back(y1);
                zhits.push_back(z1);

                xcells.push_back(x0);
                ycells.push_back(y0);
                zcells.push_back(z0);
            }
        }
    }

    return std::make_tuple(xhits, yhits, zhits, xcells, ycells, zcells);
}

int main(int argc, char* argv[]) {

  //I need to check the algorithm for the cylinder intersection... can probably solve 2D
  //Then check if within length

  TApplication app("app", &argc, argv);

  auto [xCenter,yCenter,zCenter] = GetStrawCentersTransverse();

  PlotXYCells(xCenter,yCenter);
  //PlotMidYZ(yCenter,zCenter);

 std::vector<double> trk_x,trk_y,trk_z, dirvect={1,2,0};//generateDir();

  for (double t = 0; t<20; t=t+0.1) {
    trk_x.push_back(dirvect[0]*t);
    trk_y.push_back(dirvect[1]*t);
    trk_z.push_back(dirvect[2]*t);
  }

  PlotTrack(trk_x,trk_y);

  auto [xhits, yhits, zhits, xcells, ycells, zcells] = hits(dirvect, xCenter, yCenter, zCenter);

  TGraph *graph = new TGraph();

  for (int i = 0; i < xhits.size(); ++i) {
        graph->SetPoint(i, xcells[i], ycells[i]); // Add the i-th point
    }

    graph->SetMarkerStyle(20); // Set marker style
    graph->SetMarkerSize(1);   // Set marker size
    graph->SetTitle("Example TGraph;X-axis;Y-axis");
    graph->Draw("same, P");         // "A" for axes, "P" for points


  app.Run();
  return 0;
}
