#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <TROOT.h>
#include <TH1D.h>
#include <TFile.h>

void plot() {
    // Input file name
    std::string infile = "out.dat"; // Replace with your file name

    // Define the histogram
    int bins = 30;           // Number of bins
    double xMin = -7.5;        // Minimum x-axis value
    double xMax = 7.5;       // Maximum x-axis value (adjust to your data)
    TH2D *hist = new TH2D("hist", "hist",30,-7.5,7.5,10,-2.5,2.5);

    // Open the CSV file
    std::ifstream file("out.dat");
    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file!" << std::endl;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string number;
        std::vector<double> values;

        // Parse the line to extract numbers
        while (std::getline(ss, number, ',')) {
            values.push_back(std::stod(number));
        }

        // Ensure there are at least two values in the row
        if (values.size() >= 2) {
            double x = values[0]; // First number (X-coordinate)
            double y = values[1]; // Second number (Y-coordinate)
            double z = values[2];
            cout<<(x+7.5)/0.5<<","<<(y+2.5)/0.5<<","<<z<<endl;
            hist->SetBinContent((x+8.0)/0.5, (y+3.0)/0.5, z);     // Fill the histogram
        }
    }

    hist->Draw("surf1");

}
