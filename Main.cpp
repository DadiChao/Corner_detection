//
//  main.cpp
//  corner_detection
//
//  Created by Devin on 6/26/15.
//  Copyright (c) 2015 image analysis and understanding. All rights reserved.
//

#include "CImg.h"
#include "Header.h"
#include <vector>
#include <cmath>
#include <iostream>

#define PI 3.14159265

using namespace cimg_library;
using namespace std;

/* SHIFT LOGO = SL */
/* SHIFT RICE = SR */
/* ROTATE LOGO = RL */
/* ROTATE RICE = RR */

// My project folder: /Users/Devin/Documents/IAU/Project/corner_detection/corner_detection

string infile_sl = ".../logo_for_shift.bmp";
string infile_sr = ".../rice_for_shift.bmp";
string infile_rl = ".../logo.bmp";
string infile_rr = ".../rice.bmp";

string outfile_sl = ".../result_logo_shift.bmp";
string outfile_sr = ".../result_rice_shift.bmp";
string outfile_rl = ".../result_logo_rotated.bmp";
string outfile_rr = ".../result_rice_rotated.bmp";

void SL() {
    CImg<float> img(infile_sl.c_str());
    for (int i=0; i<5; i++) { for (int j=0; j<5; j++) {int k1 = -10+i*5, k2 = -10+j*5; displacement (img, k1, k2, outfile_sl);} }
}

void SR() {
    CImg<float> img(infile_sr.c_str());img = thresholding(img, 150);
    for (int i=0; i<5; i++) { for (int j=0; j<5; j++) {int k1 = -10+i*5, k2 = -10+j*5; displacement (img, k1, k2, outfile_sr);} }
}

void RL() {
    CImg<float> img(infile_rl.c_str());
    for (int i=1; i<=1; i++) {int k = i*5; rotation(img, k, outfile_rl);}
}

void RR() {
    CImg<float> img(infile_rr.c_str());img = thresholding(img, 150);
    for (int i=1; i<=1; i++) {int k = i*5; rotation(img, k, outfile_rr);}
}


int main () {
    RR();
    return 0;
}