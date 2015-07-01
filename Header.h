//
//  Header.h
//  corner_detection
//
//  Created by Devin on 6/27/15.
//  Copyright (c) 2015 image analysis and understanding. All rights reserved.
//

#ifndef corner_detection_Header_h
#define corner_detection_Header_h

#include <vector>
#include <iostream>
#include "CImg.h"

#define PI 3.14159265

using namespace cimg_library;
using namespace std;

/* START THRESHOLDING */
CImg<float> thresholding (CImg<float> img, int thresholding_value=0) {
    //img.resize(-100,-100,-100,3);
    for (int x=0; x<img.width(); x++) {
        for (int y=0; y<img.height(); y++) {
            if (img(x,y,0,0)>thresholding_value) {
                img(x,y,0,0)=255;
                img(x,y,0,1)=255;
                img(x,y,0,2)=255;
            }
            else {
                img(x,y,0,0)=0;
                img(x,y,0,1)=0;
                img(x,y,0,2)=0;
            }
        }
    }
    return img;
}
/* END THRESHOLDING */

/* START OVERLAY */
CImg<float> overlay (CImg<float> img, CImg<float>& highlight, int color=0) {
    img.resize(-100,-100,-100,3);
    for (int x=0; x<img.width(); x++) {
        for (int y=0; y<img.height(); y++) {
            img(x,y,0,color) = max(img(x,y,0,color),highlight(x,y));
        }
    }
    return img;
}
/* END OVERLAY */


/* START CONVERSION */
/** ensure 3-color channels -> RGB */
CImg<float> RGB(CImg<float> img) { return img.resize(-100,-100,-100,3); }
/** ensure 1-color channels -> GRAY */
CImg<float> GRAY(CImg<float> img) {
    CImg<float> res = img;
    for (int x=0; x<img.width(); x++) {
        for (int y=0; y<img.height(); y++) {
            res(x,y,0,0)=img(x,y,0,0);
            res(x,y,0,1)=img(x,y,0,0);
            res(x,y,0,2)=img(x,y,0,0);
        }
    }
    return res;
}
/* END CONVERSION */


/* START DETECTION */
typedef std::vector<std::pair<int,int> > TVectorOfPairs;
void detection (CImg<float>& img, TVectorOfPairs& mark, float thresh, int halfwidth) {
    mark.clear();
    for (int y=halfwidth; y<img.height()-halfwidth; y++) {
        for (int x=halfwidth; x<img.width()-halfwidth; x++) {
            float value = img(x,y);
            if (value<thresh) { continue; }
            
            bool ismax = true;
            for (int ny=y-halfwidth; ny<=y+halfwidth; ny++) {
                for (int nx=x-halfwidth; nx<=x+halfwidth; nx++) {
                    ismax = ismax && (img(nx,ny)<=value);
                }}
            if (!ismax) continue;
            
            mark.push_back (make_pair(x,y));
        }}
}
/* END DETECTION */


/* START GAUSS_FILTER */
void gauss_filter (CImg<float>& filter, float sigma=1.0f, int deriv=0) {
    float width = 3*sigma;
    float sigma2 = sigma*sigma;
    filter.assign(int(2*width)+1);
    int i=0;
    for (float x=-width; x<=width; x+=1.0f) {
        float g = exp(-0.5*x*x/sigma2) / sqrt(2*M_PI) / sigma;
        if (deriv==1) g *= -x/sigma2;
        if (deriv==2) g *= (x*x/sigma2 - 1.0f)/sigma2;
        filter[i] = g ;
        i++;
    }
}
/* END GAUSS_FILTER */


/* START HARRIS */
//kappa (about 0.04 -- 0.15)
void harris (const CImg<float>& greyImg, CImg<float>& har,
             const float kappa=0.05, const float sigmaD=5.0, const float sigmaA=5.0) {
    // compute first derivatives
    CImg<float> Dx(greyImg), Dy(greyImg);
    CImg<float> filter, filterD;
    gauss_filter(filter, sigmaD, 0);
    gauss_filter(filterD, sigmaD, 1);
    Dx.convolve(filter).convolve(filterD.get_transpose());
    Dy.convolve(filterD).convolve(filter.get_transpose());
    
    // compute autocorrelation
    CImg<float> har11(greyImg), har12(greyImg), har22(greyImg);
    cimg_forXY (Dx, x, y) {
        const float& t1 = Dx(x,y);
        const float& t2 = Dy(x,y);
        har11(x,y) = t1*t1;
        har12(x,y) = t1*t2;
        har22(x,y) = t2*t2;
    }
    // may be done with self-defined convolutions
    gauss_filter(filter, sigmaA);
    har11.convolve(filter).convolve(filter.get_transpose());
    har12.convolve(filter).convolve(filter.get_transpose());
    har22.convolve(filter).convolve(filter.get_transpose());
    
    // compute harris response
    har.assign(greyImg); // allocate space
    //har = greyImg;
    cimg_forXY (har11, x, y) {
        float A = har11(x,y);
        float B = har22(x,y);
        float C = har12(x,y);
        float det = A*B-C*C;
        float trc = A+B;
        har(x,y) = det - kappa*trc*trc;
    }
}
/* END HARRIS */

/* START SHIFT OPERATION */
float displacement (CImg<float> img, int para_x, int para_y, string outfile_link) {
    float sigmaA = 5.0f;
    float sigmaD = 5.0f;
    float kappa  = 0.05f;
    // post processing stuf
    float thresh = 800.0f;
    
    // load image and ensure greyscale img!
    CImg<float> input_primary = img;
    input_primary = GRAY(input_primary);
    CImg<float> input_shifted = input_primary;
    input_shifted.shift(para_x,para_y,0,0,0);
    (input_primary,input_shifted).display ("Primary and Shifted Image");
    
    // apply corner detector and display result
    CImg<float> output_primary;
    CImg<float> output_shifted;
    harris (input_primary, output_primary, kappa, sigmaD, sigmaA);
    harris (input_shifted, output_shifted, kappa, sigmaD, sigmaA);
    (output_primary,output_shifted).display ("Harris Raw of Primary and Shifted Image");
    
    // for primary image
    TVectorOfPairs mark_primary;
    detection(output_primary, mark_primary, thresh, 3);
    CImg<float> corners_primary (output_primary);
    corners_primary.fill(0.f);
    float kpx=0, kpy=0;
    for (unsigned int i=0; i<mark_primary.size(); i++) {
        float color = 255.0f;
        corners_primary.draw_circle (mark_primary[i].first, mark_primary[i].second, 2, &color);
        kpx=(kpx*i+mark_primary[i].first)/(i+1);
        kpy=(kpy*i+mark_primary[i].second)/(i+1);
    }
    
    // for rotated image
    TVectorOfPairs mark_shifted;
    detection(output_shifted, mark_shifted, thresh, 3);
    CImg<float> corners_shifted (output_shifted);
    corners_shifted.fill(0.f);
    float ksx=0, ksy=0;
    for (unsigned int i=0; i<mark_shifted.size(); i++) {
        float color = 255.0f;
        corners_shifted.draw_circle (mark_shifted[i].first, mark_shifted[i].second, 2, &color);
        ksx=(ksx*i+mark_shifted[i].first)/(i+1);
        ksy=(ksy*i+mark_shifted[i].second)/(i+1);
    }
    
    float est_shift_x = ksx-kpx;
    float est_shift_y = ksy-kpy;
    cout << '(' << est_shift_x << ',' << est_shift_y << ')' << "; " << '(' << abs(est_shift_x-para_x) << ',' << abs(est_shift_y-para_y) << ')' << endl;
    cout << "Error=" << abs(est_shift_x-para_x)+abs(est_shift_y-para_y) << endl;
    (RGB(corners_primary),RGB(corners_shifted)).display("Corners");
    (overlay(input_primary,corners_primary,0), overlay(input_shifted,corners_shifted,0)).display ("Corner Detecture Results");
    
    CImgList<float> imgl (input_primary, input_shifted, output_shifted, RGB(corners_shifted), overlay(input_shifted,corners_shifted), 0);
    imgl.display ("Result");
    imgl.get_append('x').save (outfile_link.c_str());
    
    return 0;
}
/* END SHIFT OPERATION */

/* START ROTATION OPERATION */
float rotation (CImg<float> img, int para_rotate, string outfile_link) {
    float sigmaA = 5.0f;
    float sigmaD = 5.0f;
    float kappa  = 0.05f;
    // pre-scaling/rotation operations
    int angle = para_rotate;
    // post processing stuf
    float thresh = 800.0f;
    
    // load image and ensure greyscale img!
    CImg<float> input_primary = img;
    input_primary = GRAY(input_primary);
    CImg<float> input_rotated = input_primary;
    if (angle != 0) {input_rotated.rotate(angle);}
    (input_primary,input_rotated).display ("Primary and Rotated Image");
    
    // apply corner detector and display result
    CImg<float> output_primary;
    CImg<float> output_rotated;
    harris (input_primary, output_primary, kappa, sigmaD, sigmaA);
    harris (input_rotated, output_rotated, kappa, sigmaD, sigmaA);
    (output_primary,output_rotated).display ("Harris Raw of Primary and Rotated Image");
    
    // for primary image
    TVectorOfPairs mark_primary;
    detection(output_primary, mark_primary, thresh, 3);
    CImg<float> corners_primary (output_primary);
    corners_primary.fill(0.f);
    float kpx=0, kpy=0, cpx=output_primary.width()/2, cpy=output_primary.height()/2;
    for (unsigned int i=0; i<mark_primary.size(); i++) {
        float color = 255.0f;
        corners_primary.draw_circle (mark_primary[i].first, mark_primary[i].second, 2, &color);
        kpx=(kpx*i+mark_primary[i].first)/(i+1);
        kpy=(kpy*i+mark_primary[i].second)/(i+1);
    }
    
    // for rotated image
    TVectorOfPairs mark_rotated;
    detection(output_rotated, mark_rotated, thresh, 3);
    CImg<float> corners_rotated (output_rotated);
    corners_rotated.fill(0.f);
    float krx=0, kry=0, crx=output_rotated.width()/2, cry=output_rotated.height()/2;
    for (unsigned int i=0; i<mark_rotated.size(); i++) {
        float color = 255.0f;
        corners_rotated.draw_circle (mark_rotated[i].first, mark_rotated[i].second, 2, &color);
        krx=(krx*i+mark_rotated[i].first)/(i+1);
        kry=(kry*i+mark_rotated[i].second)/(i+1);
    }
    
    float vpx=kpx-cpx, vpy=kpy-cpy, vrx=krx-crx, vry=kry-cry;
    float cross_product = vpx*vrx + vpy*vry;
    float module_product = sqrt((vpx*vpx+vpy*vpy)*(vrx*vrx+vry*vry));
    float cos_result = cross_product/module_product;
    float result_rotation = acos (cos_result)/PI*180;
    cout << result_rotation << endl;
    float error = 1 - abs(result_rotation-(float)angle) / (float)angle;
    cout<<error<<endl;
    
    (RGB(corners_primary),RGB(corners_rotated)).display("Corners");
    (overlay(input_primary,corners_primary,0), overlay(input_rotated,corners_rotated,0)).display ("Corner Detecture Results");
    
    CImgList<float> imgl (input_primary, input_rotated, output_rotated, RGB(corners_rotated), overlay(input_rotated,corners_rotated), 0);
    imgl.display ("Result");
    imgl.get_append('x').save (outfile_link.c_str());
    
    return error;
}
/* END ROTATION OPERATION */


#endif
