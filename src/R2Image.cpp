// Source file for image class



// Include files 

#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"
#include "iostream"
#include "vector"



////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////


R2Image::
R2Image(void)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
}



R2Image::
R2Image(const char *filename)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
  // Read image
  Read(filename);
}



R2Image::
R2Image(int width, int height)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
  : pixels(NULL),
    npixels(image.npixels),
    width(image.width), 
    height(image.height)
    
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
  // Delete previous pixels
  if (pixels) { delete [] pixels; pixels = NULL; }

  // Reset width and height
  npixels = image.npixels;
  width = image.width;
  height = image.height;

  // Allocate new pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];

  // Return image
  return *this;
}


void R2Image::
svdTest(void)
{
  // fit a 2D conic to five points
  R2Point p1(1.2,3.5);
  R2Point p2(2.1,2.2);
  R2Point p3(0.2,1.6);
  R2Point p4(0.0,0.5);
  R2Point p5(-0.2,4.2);

  // build the 5x6 matrix of equations
  double** linEquations = dmatrix(1,5,1,6);

  linEquations[1][1] = p1[0]*p1[0];
  linEquations[1][2] = p1[0]*p1[1];
  linEquations[1][3] = p1[1]*p1[1];
  linEquations[1][4] = p1[0];
  linEquations[1][5] = p1[1];
  linEquations[1][6] = 1.0;

  linEquations[2][1] = p2[0]*p2[0];
  linEquations[2][2] = p2[0]*p2[1];
  linEquations[2][3] = p2[1]*p2[1];
  linEquations[2][4] = p2[0];
  linEquations[2][5] = p2[1];
  linEquations[2][6] = 1.0;

  linEquations[3][1] = p3[0]*p3[0];
  linEquations[3][2] = p3[0]*p3[1];
  linEquations[3][3] = p3[1]*p3[1];
  linEquations[3][4] = p3[0];
  linEquations[3][5] = p3[1];
  linEquations[3][6] = 1.0;
  
  linEquations[4][1] = p4[0]*p4[0];
  linEquations[4][2] = p4[0]*p4[1];
  linEquations[4][3] = p4[1]*p4[1];
  linEquations[4][4] = p4[0];
  linEquations[4][5] = p4[1];
  linEquations[4][6] = 1.0;

  linEquations[5][1] = p5[0]*p5[0];
  linEquations[5][2] = p5[0]*p5[1];
  linEquations[5][3] = p5[1]*p5[1];
  linEquations[5][4] = p5[0];
  linEquations[5][5] = p5[1];
  linEquations[5][6] = 1.0;

  printf("\n Fitting a conic to five points:\n");
  printf("Point #1: %f,%f\n",p1[0],p1[1]);
  printf("Point #2: %f,%f\n",p2[0],p2[1]);
  printf("Point #3: %f,%f\n",p3[0],p3[1]);
  printf("Point #4: %f,%f\n",p4[0],p4[1]);
  printf("Point #5: %f,%f\n",p5[0],p5[1]);

  // compute the SVD
  double singularValues[7]; // 1..6
  double** nullspaceMatrix = dmatrix(1,6,1,6);
  

  // get the result
  printf("\n Singular values: %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6]);

  // find the smallest singular value:
  int smallestIndex = 1;
  for(int i=2;i<7;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

  // solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
  printf("Conic coefficients: %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex]);

  // make sure the solution is correct:
  printf("Equation #1 result: %f\n",  p1[0]*p1[0]*nullspaceMatrix[1][smallestIndex] + 
                    p1[0]*p1[1]*nullspaceMatrix[2][smallestIndex] + 
                    p1[1]*p1[1]*nullspaceMatrix[3][smallestIndex] + 
                    p1[0]*nullspaceMatrix[4][smallestIndex] + 
                    p1[1]*nullspaceMatrix[5][smallestIndex] + 
                    nullspaceMatrix[6][smallestIndex]);

  printf("Equation #2 result: %f\n",  p2[0]*p2[0]*nullspaceMatrix[1][smallestIndex] + 
                    p2[0]*p2[1]*nullspaceMatrix[2][smallestIndex] + 
                    p2[1]*p2[1]*nullspaceMatrix[3][smallestIndex] + 
                    p2[0]*nullspaceMatrix[4][smallestIndex] + 
                    p2[1]*nullspaceMatrix[5][smallestIndex] + 
                    nullspaceMatrix[6][smallestIndex]);

  printf("Equation #3 result: %f\n",  p3[0]*p3[0]*nullspaceMatrix[1][smallestIndex] + 
                    p3[0]*p3[1]*nullspaceMatrix[2][smallestIndex] + 
                    p3[1]*p3[1]*nullspaceMatrix[3][smallestIndex] + 
                    p3[0]*nullspaceMatrix[4][smallestIndex] + 
                    p3[1]*nullspaceMatrix[5][smallestIndex] + 
                    nullspaceMatrix[6][smallestIndex]);

  printf("Equation #4 result: %f\n",  p4[0]*p4[0]*nullspaceMatrix[1][smallestIndex] + 
                    p4[0]*p4[1]*nullspaceMatrix[2][smallestIndex] + 
                    p4[1]*p4[1]*nullspaceMatrix[3][smallestIndex] + 
                    p4[0]*nullspaceMatrix[4][smallestIndex] + 
                    p4[1]*nullspaceMatrix[5][smallestIndex] + 
                    nullspaceMatrix[6][smallestIndex]);

  printf("Equation #5 result: %f\n",  p5[0]*p5[0]*nullspaceMatrix[1][smallestIndex] + 
                    p5[0]*p5[1]*nullspaceMatrix[2][smallestIndex] + 
                    p5[1]*p5[1]*nullspaceMatrix[3][smallestIndex] + 
                    p5[0]*nullspaceMatrix[4][smallestIndex] + 
                    p5[1]*nullspaceMatrix[5][smallestIndex] + 
                    nullspaceMatrix[6][smallestIndex]);

  R2Point test_point(0.34,-2.8);

  printf("A point off the conic: %f\n", test_point[0]*test_point[0]*nullspaceMatrix[1][smallestIndex] + 
                      test_point[0]*test_point[1]*nullspaceMatrix[2][smallestIndex] + 
                      test_point[1]*test_point[1]*nullspaceMatrix[3][smallestIndex] + 
                      test_point[0]*nullspaceMatrix[4][smallestIndex] + 
                      test_point[1]*nullspaceMatrix[5][smallestIndex] + 
                      nullspaceMatrix[6][smallestIndex]);

  return; 
}



////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Brighten(double factor)
{
  // Brighten the image by multiplying each pixel component by the factor.
  // This is implemented for you as an example of how to access and set pixels
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j) *= factor;
      Pixel(i,j).Clamp();
    }
  }
}

void R2Image::
SobelX(void)
{
  // Apply the Sobel oprator to the image in X direction
  R2Image tempImage(*this); // Copy the image

  for (int y = 1; y < height - 1; y++){
    for (int x = 1; x < width - 1; x++){

      Pixel(x, y) = tempImage.Pixel(x-1, y-1) + tempImage.Pixel(x-1, y) + tempImage.Pixel(x-1, y)+ tempImage.Pixel(x-1, y+1);

      Pixel(x, y) = Pixel(x, y) - tempImage.Pixel(x+1, y-1) - tempImage.Pixel(x+1, y) - tempImage.Pixel(x+1, y) - tempImage.Pixel(x+1, y+1);

      //Pixel(x, y).Clamp();
    }

  }

}

void R2Image::
SobelY(void)
{
  // Apply the Sobel oprator to the image in Y direction
  R2Image tempImage(*this);

  for (int y = 1; y < height - 1; y++){
    for (int x = 1; x < width - 1; x++){

      Pixel(x, y) = tempImage.Pixel(x-1, y-1) + tempImage.Pixel(x, y-1) + tempImage.Pixel(x, y-1)+ tempImage.Pixel(x+1, y-1);

      Pixel(x, y) = Pixel(x, y) - tempImage.Pixel(x-1, y+1) - tempImage.Pixel(x, y+1) - tempImage.Pixel(x, y+1) - tempImage.Pixel(x+1, y+1);

      //Pixel(x, y).Clamp();

    }
  }
}

//sigmaR is for color difference, sigmaD is for distance difference

void R2Image::
Bilateral(double sigmaR, double sigmaD)
{

  R2Image temp(*this);

  int kernelDSize = sigmaD * 3;
  double kernelR[901];
  double kernelD[kernelDSize];

  for (int i = 0; i< 901; i++){
    kernelR[i] =  (1/(sqrt(2*acos(-1)*pow(sigmaR, 2))))*exp(-pow(i,2)/(2*pow(sigmaR,2)));
  }

  for (int i = -kernelDSize; i < kernelDSize; i++){
    kernelD[i+kernelDSize] = (1/(sqrt(2*acos(-1)*pow(sigmaD, 2))))*exp(-pow(i,2)/(2*pow(sigmaD,2)));
  }

  for (int x = kernelDSize; x < (width - kernelDSize)/3; x++){
    printf("x: %d\n",x);
    for (int y = kernelDSize; y < height - kernelDSize; y++){
      double total = 0.0;
      Pixel(x,y).Reset(0.0,0.0,0.0,0.0);
      double pixelIntensity = temp.Pixel(x,y)[0] + temp.Pixel(x,y)[1] + temp.Pixel(x,y)[2];
      for (int a = -kernelDSize; a < kernelDSize; a++){
        for (int b = -kernelDSize; b < kernelDSize; b++){
          double currIntensity = temp.Pixel(x+a, y+b)[0] + temp.Pixel(x+a,y+b)[1] + temp.Pixel(x+a, y+b)[2];
          double colorDifference = pow(pixelIntensity - currIntensity,2);
          double weight = kernelD[a+kernelDSize]*kernelD[b+kernelDSize]*kernelR[int(colorDifference*100)];
          total += weight;
          Pixel(x,y) += temp.Pixel(x+a,y+b)*weight;

        }
      }

      Pixel(x,y) /= total;
      

    }
  }
}

struct PixelDifference
{
  R2Pixel pixel;
  double difference;

  PixelDifference(R2Pixel newP, double diff){
    pixel = newP;
    difference = diff;
  }

  bool operator<(const PixelDifference& pd) const{
    return difference < pd.difference;
  }
};

void R2Image::
Median(int kernelSize, double sigma)
{
  R2Image temp(*this);

  // int cutSize = int(kernelSize * sigma);
  int halfKernel = int(kernelSize/2);

  // for (int x = halfKernel; x < (int)(width - halfKernel)/3; x++){
  //   for (int y = halfKernel; y < height - halfKernel; y++){
  //     std::vector<PixelDifference> currKernel;
  //     for (int a = -halfKernel; a < halfKernel+1; a++){
  //       for (int b = -halfKernel; b < halfKernel+1; b++){
  //         R2Pixel curr = temp.Pixel(x+a, y+b);
  //         double difference = pow(curr[0],2)+pow(curr[1],2)+pow(curr[2],2);
  //         currKernel.push_back(PixelDifference(curr, difference));
  //       }
  //     }
  //     std::sort(currKernel.begin(), currKernel.end());

  //     std::vector<R2Pixel> acceptPixels;
  //     for (int i = cutSize; i< currKernel.size()-cutSize; i++){
  //       acceptPixels.push_back(currKernel.at(i).pixel);
  //     }

  //     Pixel(x,y).Reset(0.0,0.0,0.0,0.0);
  //     int num = 0;
  //     for (int j = 0; j<acceptPixels.size(); j++){
  //       Pixel(x,y) += acceptPixels.at(j);
  //       num +=1;
  //     }
  //     Pixel(x,y)/=num;
  //   }
  // }

  for (int x = halfKernel; x <(int)(width - halfKernel)/3; x++){
    for (int y = halfKernel; y < height - halfKernel; y++){
      std::vector<PixelDifference> currKernel;
      for (int a = -halfKernel; a < halfKernel+1; a++){
        for (int b = -halfKernel; b<halfKernel+1; b++){
          R2Pixel curr = temp.Pixel(x+a,y+b);
          double difference = pow(curr[0],2)+pow(curr[1],2)+pow(curr[2],2);
          currKernel.push_back(PixelDifference(curr, difference));
        }
      }
      std::sort(currKernel.begin(), currKernel.end());
      R2Pixel newPixel = currKernel.at((int)(currKernel.size()/2)).pixel;
      Pixel(x,y) = newPixel;
    }
  }


}

void R2Image::
LoG(void)
{
  // Apply the LoG oprator to the image
  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "LoG() not implemented\n");
}



void R2Image::
ChangeSaturation(double factor)
{
  // Changes the saturation of an image
  // Find a formula that changes the saturation without affecting the image brightness

  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "ChangeSaturation(%g) not implemented\n", factor);
}

void R2Image::
FeatureDetection(){
  return;
}

void R2Image::
blendOtherImageHomography(R2Image * otherImage)
{
  return;
}


// Linear filtering ////////////////////////////////////////////////
void R2Image::
Blur(double sigma)
{
  R2Image tempImage(*this);

  int kernelSize = sigma * 3;
  int kernelWidth = sigma * 6 + 1;
  double gaussianKernel[kernelWidth];
  double total = 0.0;
  double curr = 0.0;
  double offset = (1/(sqrt(2*acos(-1)*pow(sigma, 2))))*exp(-pow(kernelSize,2)/(2*pow(sigma,2)));

  for (int i = 1; i< kernelSize+1; i++){
    curr =  (1/(sqrt(2*acos(-1)*pow(sigma, 2))))*exp(-pow(i,2)/(2*pow(sigma,2))) - offset;
    total+= curr;
  }

  total *= 2;
  total += (1/(sqrt(2*acos(-1)*pow(sigma, 2))))*exp(-pow(0,2)/(2*pow(sigma,2))) - offset;

  for (int i = - kernelSize; i < kernelSize + 1; i++){
    gaussianKernel[i+kernelSize] = ((1/(sqrt(2*acos(-1)*pow(sigma, 2))))*exp(-pow(i,2)/(2*pow(sigma,2))) - offset)/total;
  }


  for (int x = kernelSize; x < width - kernelSize; x++){
    for (int y = kernelSize; y < height - kernelSize; y++){
      Pixel(x,y).Reset(0.0, 0.0, 0.0, 0.0);

      for (int a = -kernelSize; a < kernelSize + 1; a++){
        for (int b = -kernelSize; b<kernelSize + 1; b++){
          Pixel(x,y) = Pixel(x,y) + tempImage.Pixel(x+a,y+b)*gaussianKernel[a+kernelSize]*gaussianKernel[b+kernelSize]; 
        }  
      }

      //Pixel(x, y).Clamp();
    }
  }
}


void R2Image::
Harris(double sigma)
{
  // Harris corner detector. Make use of the previously developed filters, such as the Gaussian blur filter
  // Output should be 50% grey at flat regions, white at corners and black/dark near edges
  R2Image sobelx(*this);
  R2Image sobely(*this);

  R2Image Img1(*this);
  R2Image Img2(*this);
  R2Image Img3(*this);

  sobelx.SobelX();
  sobely.SobelY();

  for (int x = 0; x < width; x++){
    for (int y = 0; y < height; y++){
      Img3.Pixel(x,y) = sobelx.Pixel(x,y)*sobely.Pixel(x,y);
      Img1.Pixel(x,y) = sobelx.Pixel(x,y)*sobelx.Pixel(x,y);
      Img2.Pixel(x,y) = sobely.Pixel(x,y)*sobely.Pixel(x,y);
    }
  }
  Img1.Blur(sigma);
  Img2.Blur(sigma);
  Img3.Blur(sigma);

  R2Pixel *pixel = new R2Pixel(0.5, 0.5, 0.5, 1);

  for (int x = 0; x < width; x++){
    for (int y = 0; y < height; y++){
      Pixel(x,y) = Img1.Pixel(x,y)*Img2.Pixel(x,y)-Img3.Pixel(x,y)*Img3.Pixel(x,y)-0.04*((Img1.Pixel(x,y)+Img2.Pixel(x,y))*(Img1.Pixel(x,y)+Img2.Pixel(x,y)));
      Pixel(x,y) += *pixel;
      Pixel(x,y).Clamp();
    }
  }

  
}


struct Feature
{
  // create a structure with custom data, and how to sort them fast.
  int centerX;
  int centerY;
  R2Pixel HarrisValue;

  Feature(int x, int y, R2Pixel val){
    centerX = x;
    centerY = y;
    HarrisValue = val;
  }

  bool operator<(const Feature& feature) const{
    double valueIntensity = HarrisValue[0]+HarrisValue[1]+HarrisValue[2];
    double featureIntensity = feature.HarrisValue[0]+feature.HarrisValue[1]+feature.HarrisValue[2];
    return valueIntensity < featureIntensity;
  }

  double difference(const Feature& feature) const{
    double valueIntensity = HarrisValue[0]+HarrisValue[1]+HarrisValue[2];
    double featureIntensity = feature.HarrisValue[0]+feature.HarrisValue[1]+feature.HarrisValue[2];
    return (valueIntensity - featureIntensity)*(valueIntensity - featureIntensity);
  }

  double charScale(const Feature& feature) const{
    double valueIntensity = HarrisValue[0]+HarrisValue[1]+HarrisValue[2];
    double featureIntensity = feature.HarrisValue[0]+feature.HarrisValue[1]+feature.HarrisValue[2];
    return (double)pow(valueIntensity,2)/(double)pow(featureIntensity,2);
  }

};

void R2Image::
FirstFrameProcessing()
{
  // Find 150 features with high corner score, and make sure there is at least 10 pixel distance between them 
  std::vector<Feature> featureVec;

  R2Image harrisImage(*this);
  harrisImage.Harris(2);

  R2Pixel *redPixel = new R2Pixel(1.0, 0.0, 0.0, 1.0);
  bool flag;

  for (int x = 0; x < width; x++){
    for (int y = 0; y < height; y++){
      featureVec.push_back(Feature(x,y,harrisImage.Pixel(x,y)));
    }
  }

  std::sort(featureVec.begin(), featureVec.end());

  int i = 0;
  while (i<150 && featureVec.size() > 0){
    Feature curr = featureVec.back();
    featureVec.pop_back();

    flag = true;
    
    for (int b = -5; b < 6; b++){
      if (Pixel(curr.centerX+b, curr.centerY+b) == *redPixel){
        flag = false;
      }
    }
    
    if (flag){
      for (int a = -5; a < 6; a++){
        Pixel(curr.centerX + a, curr.centerY - 5) = *redPixel;
        Pixel(curr.centerX + a, curr.centerY + 5) = *redPixel;
        Pixel(curr.centerX - 5, curr.centerY + a) = *redPixel;
        Pixel(curr.centerX + 5, curr.centerY + a) = *redPixel;
        
      }
      i++;
    }
  }


}



void R2Image::
Sharpen()
{
  // Sharpen an image using a linear filter. Use a kernel of your choosing.
  R2Image blur(*this);
  R2Image temp(*this);

  // for (int x = 1; x < width -1; x++){
  //   for (int y = 1; y < height -1; y++){
  //     Pixel(x,y) *= 5;
  //     Pixel(x,y) = Pixel(x,y)- tempImage.Pixel(x, y-1);
  //     Pixel(x,y) = Pixel(x,y) - tempImage.Pixel(x-1,y) - tempImage.Pixel(x+1, y);
  //     Pixel(x,y) = Pixel(x,y) - tempImage.Pixel(x, y+1);

  //     Pixel(x,y).Clamp();
  //   }
  // }

  blur.Blur(5);
  for (int x = 0; x < width; x++){
    for (int y = 0; y < height; y++){
      Pixel(x,y) += temp.Pixel(x,y) - blur.Pixel(x,y);
    }
  }

}

void R2Image::
LensDistortion()
{
  R2Image temp(*this);

  for(int i = 0; i < width; i++){
    for (int j = 0; j < height; j++){
      int x = i - width/2;
      int y = j - height/2;

      double radius = sqrt(pow(x,2)+pow(y,2));
      double factor = sqrt(1-(pow(radius/width, 1.5)));

      int xold = (int)x/factor/1.25;
      int yold = (int)y/factor/1.25;
      int iold = width/2+xold;
      int jold = height/2+yold;
      if (iold<0){
        iold = 0;
      }
      if (iold >= width){
        iold = width -1;
      }
      if (jold < 0){
        jold =0;
      }
      if (jold>= height){
        jold = height -1;
      }
      Pixel(i, j) = temp.Pixel(iold,jold);
    }
  }

}

struct FeaturePair
{
  // create a structure with custom data, and how to sort them fast.
  int feat1X;
  int feat1Y;
  int feat2X;
  int feat2Y;
  int vectorX;
  int vectorY;
  double inlier;

  FeaturePair(int x1, int y1, int x2, int y2, int x, int y, double in){
    feat1X = x1;
    feat1Y = y1;
    feat2X = x2;
    feat2Y = y2;
    vectorX = x;
    vectorY = y;
    inlier = in;
  }

  double difference(const FeaturePair& pair) const{
    int subtractX = vectorX - pair.vectorX;
    int subtractY = vectorY - pair.vectorY;
    return subtractY * subtractY + subtractX * subtractX;
  }

  bool operator<(const FeaturePair& currPair) const{
    return inlier < currPair.inlier;
  }

};


void R2Image::
ScaleInvariantHarris(){
  for (int i = 5; i < 20; i++){
    R2Image harrisImage(*this);
    harrisImage.Harris(i);
    std::vector<Feature> featureVec;

    for (int x = 0; x < width; x++){
      for (int y = 0; y < height; y++){
        featureVec.push_back(Feature(x,y,harrisImage.Pixel(x,y)));
      }
    }
    std::sort(featureVec.begin(), featureVec.end());
    std::cout << "featureVec size is " << featureVec.size() << std::endl;
    Feature firstFeat = featureVec.back();
    featureVec.pop_back();
    // std::cout << "featureVec size is " << featureVec.size() << std::endl;

    int featureX = firstFeat.centerX;
    int featureY = firstFeat.centerY;
    std::cout << "featureX is " << featureX << std::endl;
    std::cout << "featureY is " << featureY << std::endl;
    Feature secondFeat = featureVec.back();
    std::cout << "featureX is " << secondFeat.centerX << std::endl;
    std::cout << "featureY is " << secondFeat.centerY << std::endl;
    double charScale = firstFeat.charScale(secondFeat);
    std::cout << "charScale is " << charScale << std::endl;
    square(featureX-charScale, featureY+charScale, featureX+charScale, featureY+charScale, featureX+charScale, featureY-charScale, featureX-charScale, featureY-charScale);
  }

}


void R2Image::
square(int a0, int a1, int b0, int b1, int c0, int c1, int d0, int d1)
{
  line(a0, a1, b0, b1, 1.0, 0.0, 0.0);
  line(b0, b1, c0, c1, 1.0, 0.0, 0.0);
  line(c0, c1, d0, d1, 1.0, 0.0, 0.0);
  line(d0, d1, a0, a1, 1.0, 0.0, 0.0);
}

void R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
  // find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
  // compute the matching translation (pixel precision is OK), and blend the translated "otherImage" 
  // into this image with a 50% opacity.

  R2Image image1 (*this);
  R2Image image2 (*otherImage);

  // finding features
  image1.Harris(2);

  std::vector<Feature> featureVec;
  std::vector<Feature> feature1;
  std::vector<Feature> feature2;
  std::vector<FeaturePair> featurePairCurr;

  R2Pixel *redPixel = new R2Pixel(1.0, 0.0, 0.0, 1.0);
  R2Pixel *greenPixel = new R2Pixel(0.0, 1.0, 0.0, 1.0);
  bool flag;

  // setting the size of search window
  int searchWidthHalf = 0.1*width;
  int searchHeightHalf = 0.1*height;

  for (int x = searchWidthHalf; x < width - searchWidthHalf; x++){
    for (int y = searchHeightHalf; y < height -searchHeightHalf; y++){
      featureVec.push_back(Feature(x,y,image1.Pixel(x,y)));
    }
  }

  std::sort(featureVec.begin(), featureVec.end());
  int i = 0;


  //for each feature detected
  while (i<150 && featureVec.size() > 0){
    
    // test printing
    std::cout << "featureVec size is %d" << featureVec.size() << std::endl;

    Feature curr = featureVec.back();
    featureVec.pop_back();
    flag = true;
    for (int b = -5; b < 6; b++){
      if (image1.Pixel(curr.centerX+b, curr.centerY+b) == *redPixel){
        flag = false;
      }
      if (image1.Pixel(curr.centerX+10-b, curr.centerY+b) == *redPixel){
        flag = false;
      }
    }
    
    if (flag){
      // drawing the first feature 
      for (int a = -5; a < 6; a++){
        image1.Pixel(curr.centerX + a, curr.centerY - 5) = *redPixel;
        image1.Pixel(curr.centerX + a, curr.centerY + 5) = *redPixel;
        image1.Pixel(curr.centerX - 5, curr.centerY + a) = *redPixel;
        image1.Pixel(curr.centerX + 5, curr.centerY + a) = *redPixel;
      }
        
      
      double bestDiff = 1000.0;  
      int bestX = -10;
      int bestY = -10;    
      // loop through all coordinates within the search window
      for (int searchX = -searchWidthHalf; searchX< searchWidthHalf; searchX++){
        for (int searchY = -searchHeightHalf; searchY<searchHeightHalf; searchY++){

          //the coordinates of second feature
          int currX = curr.centerX + searchX;
          int currY = curr.centerY + searchY;
          //compare a small window in A at (curr.centerX, curr.centerY) with a small window in B
          // at (currX, currY)

          double currDiff = 0;

          for (int ssdX = -5; ssdX < 6; ssdX++){
            for (int ssdY = -5; ssdY < 6; ssdY++){
              Feature second = Feature(currX+ssdX, currY+ssdY, image2.Pixel(currX+ssdX, currY+ssdY));
              Feature first = Feature(curr.centerX+ssdX,curr.centerY+ssdY, Pixel(curr.centerX+ssdX, curr.centerY+ssdY));
              currDiff += first.difference(second);
            }
          }

          if (currDiff < bestDiff){
            bestDiff = currDiff;
            bestX = currX;
            bestY = currY;
          }
        }
      }
      if (bestX >0 && bestY > 0){
        feature1.push_back(Feature(curr.centerX, curr.centerY, Pixel(curr.centerX, curr.centerY))); 
        feature2.push_back(Feature(bestX, bestY, image2.Pixel(bestX, bestY)));
        int vectorx = abs(curr.centerX - bestX);
        int vectory = abs(curr.centerY - bestY);
        featurePairCurr.push_back(FeaturePair(curr.centerX, curr.centerY, bestX, bestY, vectorx, vectory, 0));
      }
      i++;
    }
  }

  // draw the new graph
  for(int x = 0; x < width; x++){
    for (int y = 0; y < height; y++){
      Pixel(x,y) = 0.5*Pixel(x,y) + 0.5*image2.Pixel(x,y);
    }
  }

  // comparing feature pairs
  // picking a pair each time
  bool found = false;
  FeaturePair bestPair = featurePairCurr.back();
  for (int a = 0; a < featurePairCurr.size(); a++){
    FeaturePair first = featurePairCurr.at(a);
    int inlierNum = -1;

    for (int b = 0; b < featurePairCurr.size(); b++){
      FeaturePair curr = featurePairCurr.at(b);
      if (first.difference(curr)<10){
        inlierNum++;
      }
    }
    if (inlierNum > 130){
      found = true;
      bestPair = first;
    }
  }
  if(!found){
    std::sort(featurePairCurr.begin(), featurePairCurr.end());
    bestPair = featurePairCurr.at(0);
  }

  for (int i = 0; i < featurePairCurr.size(); i++){
    FeaturePair currPair= featurePairCurr.at(i);
    int x1 = currPair.feat1X;
    int y1 = currPair.feat1Y;
    int x2 = currPair.feat2X;
    int y2 = currPair.feat2Y;
    if (bestPair.difference(currPair)<10){
      // decide the color to draw
      line(x1, x2, y1, y2, 0.0, 1.0, 0.0);
      for(int a = -5; a < 6; a++){
        Pixel(x1 + a, y1 - 5) = *greenPixel;
        Pixel(x1 + a, y1 + 5) = *greenPixel;
        Pixel(x1 - 5, y1 + a) = *greenPixel;
        Pixel(x1 + 5, y1 + a) = *greenPixel;
        Pixel(x2 + a, y2 - 5) = *greenPixel;
        Pixel(x2 + a, y2 + 5) = *greenPixel;
        Pixel(x2 - 5, y2 + a) = *greenPixel;
        Pixel(x2 + 5, y2 + a) = *greenPixel;
      }
    }else{
      for(int a = -5; a < 6; a++){
        Pixel(x1 + a, y1 - 5) = *redPixel;
        Pixel(x1 + a, y1 + 5) = *redPixel;
        Pixel(x1 - 5, y1 + a) = *redPixel;
        Pixel(x1 + 5, y1 + a) = *redPixel;
        Pixel(x2 + a, y2 - 5) = *redPixel;
        Pixel(x2 + a, y2 + 5) = *redPixel;
        Pixel(x2 - 5, y2 + a) = *redPixel;
        Pixel(x2 + 5, y2 + a) = *redPixel;
      }
      line(x1, x2, y1, y2, 1.0, 0.0, 0.0);
    }
    
  }


}


void R2Image::
line(int x0, int x1, int y0, int y1, float r, float g, float b)
{
  // helper function to draw the lines
  if (x0 > x1){
    int x = y1;
    y1 = y0;
    y0 = x;

    x = x1;
    x1 = x0;
    x0 = x;
  }

  int deltax = x1 - x0;
  int deltay = y1 - y0;
  float error = 0;
  float deltaerr = 0.0;

  if (deltax != 0) deltaerr = fabs(float(float(deltay)/deltax));
  // note that this division needs to be done in a way that preserves the fractional part.
  int y = y0;

  for (int x = x0; x <= x1; x++){
    Pixel(x,y).Reset(r,g,b,1.0);
    error = error + deltaerr;
    if (error >= 0.5){
      if (deltay >0) y = y + 1;
      else y = y -1;
      error = error - 1.0;
    }
  } 
  if (x0 > 3 && x0 < width-3 && y0>3 && y0<height-3){
    for (int x = x0-3; x<= x0+3; x++){
      for (int y = y0 -3 ; y <= y0+3; y++){
        Pixel(x,y).Reset(r,g,b,1.0);
      }
    }
  }
}

// void R2Image::
// computeHomography(void)
// {
//   R2Point points[] = {p1, p2, p3, p4, p5, p6, p7, p8};
//   double** linEquations = dmatrix(1, 8, 1, 9);
//   }
//   for (int i = 0 ; i < 4; i++){
//     linEquations[i*2+1][1] = -points[i*2][0];
//     linEquations[i*2+1][2] = -points[i*2][1];
//     linEquations[i*2+1][3] = -1.0;
//     linEquations[i*2+1][4] = 0.0;
//     linEquations[i*2+1][5] = 0.0;
//     linEquations[i*2+1][6] = 0.0;
//     linEquations[i*2+1][7] = points[i*2][0]*points[i*2+1][0];
//     linEquations[i*2+1][8] = points[i*2][1]*points[i*2+1][0];
//     linEquations[i*2+1][9] = points[i*2+1][0];

//     linEquations[i*2+2][1] = 0.0;
//     linEquations[i*2+2][2] = 0.0;
//     linEquations[i*2+2][3] = 0.0;
//     linEquations[i*2+2][4] = -points[i*2][0];
//     linEquations[i*2+2][5] = -points[i*2][1];
//     linEquations[i*2+2][6] = -1.0;
//     linEquations[i*2+2][7] = points[i*2][0]*points[i*2+1][1];
//     linEquations[i*2+2][8] = points[i*2][1]*points[i*2+1][1];
//     linEquations[i*2+2][9] = points[i*2+1][1];
//   }

//   double singularValues[10];
//   double** nullspaceMatrix = dmatrix(1, 9, 1, 9);
//   svdcmp(linEquations, 8, 9, singularValues, nullspaceMatrix);

//   printf("\n Singular values: %f, %f, %f, %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],
//     singularValues[4],singularValues[5],singularValues[6],singularValues[7],singularValues[8],singularValues[9]);

//   // find the smallest singular value:
//   int smallestIndex = 1;
//   for(int i=2;i<10;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

//   // solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
//   printf("Conic coefficients: %f, %f, %f, %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],
//     nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex], 
//     nullspaceMatrix[7][smallestIndex],nullspaceMatrix[8][smallestIndex],nullspaceMatrix[9][smallestIndex]);

//   return;
// }

void R2Image::
FrameProcessing(R2Image * otherImage)
{
  // find at least 100 features on this image, and another 100 on the "otherImage". Based on these,
  // compute the matching homography, and blend the transformed "otherImage" into this image with a 50% opacity.

  R2Image image1 (*this);
  R2Image temp(*this);
  R2Image image2 (*otherImage);

  // finding features
  image1.Harris(2);

  std::vector<Feature> featureVec;
  std::vector<Feature> feature1;
  std::vector<Feature> feature2;

  R2Pixel *redPixel = new R2Pixel(1.0, 0.0, 0.0, 1.0);
  R2Pixel *greenPixel = new R2Pixel(0.0, 1.0, 0.0, 1.0);
  bool flag;
  int searchWidthHalf = 0.1*width;
  int searchHeightHalf = 0.1*height;
  for (int x = searchWidthHalf; x < width - searchWidthHalf; x++){
    for (int y = searchHeightHalf; y < height -searchHeightHalf; y++){
      featureVec.push_back(Feature(x,y,image1.Pixel(x,y)));
    }
  }
  std::sort(featureVec.begin(), featureVec.end());
  int i = 0;
  //for each feature detected
  while (i<150 && featureVec.size() > 0){
    
    // test printing
    //std::cout << "featureVec size is %d" << featureVec.size() << std::endl;

    Feature curr = featureVec.back();
    featureVec.pop_back();
    flag = true;
    for (int b = -5; b < 6; b++){
      if (image1.Pixel(curr.centerX+b, curr.centerY+b) == *redPixel){
        flag = false;
      }
      if (image1.Pixel(curr.centerX+10-b, curr.centerY+b) == *redPixel){
        flag = false;
      }
    }
    
    if (flag){
      // drawing the first feature 
      for (int a = -5; a < 6; a++){
        image1.Pixel(curr.centerX + a, curr.centerY - 5) = *redPixel;
        image1.Pixel(curr.centerX + a, curr.centerY + 5) = *redPixel;
        image1.Pixel(curr.centerX - 5, curr.centerY + a) = *redPixel;
        image1.Pixel(curr.centerX + 5, curr.centerY + a) = *redPixel;
      }
      double bestDiff = 1000.0;  
      int bestX = -10;
      int bestY = -10;    
      // loop through all coordinates within the search window
      for (int searchX = -searchWidthHalf; searchX< searchWidthHalf; searchX++){
        for (int searchY = -searchHeightHalf; searchY<searchHeightHalf; searchY++){
          int currX = curr.centerX + searchX;
          int currY = curr.centerY + searchY;
          double currDiff = 0;
          for (int ssdX = -5; ssdX < 6; ssdX++){
            for (int ssdY = -5; ssdY < 6; ssdY++){
              Feature second = Feature(currX+ssdX, currY+ssdY, image2.Pixel(currX+ssdX, currY+ssdY));
              Feature first = Feature(curr.centerX+ssdX,curr.centerY+ssdY, Pixel(curr.centerX+ssdX, curr.centerY+ssdY));
              currDiff += first.difference(second);
            }
          }
          if (currDiff < bestDiff){
            bestDiff = currDiff;
            bestX = currX;
            bestY = currY;
          }
        }
      }
      if (bestX >0 && bestY > 0){
        feature1.push_back(Feature(curr.centerX, curr.centerY, Pixel(curr.centerX, curr.centerY))); 
        feature2.push_back(Feature(bestX, bestY, image2.Pixel(bestX, bestY)));
      }
      i++;
    }
  }

  double** HMatrix = dmatrix(1,3,1,3);
  double** newHMatrix = dmatrix(1,3,1,3);

  int inlierNum = 0;

  for (int trial = 0; trial < 10000; trial++){
    std::vector<Feature> random1;
    std::vector<Feature> random2;
    for (int i = 0; i < 4; i++){
      int randomIndex = int(rand()%150);
      random1.push_back(feature1.at(randomIndex));
      random2.push_back(feature2.at(randomIndex));
    }

    R2Point p1(random1.at(0).centerX, random1.at(0).centerY);
    R2Point p2(random2.at(0).centerX, random2.at(0).centerY);
    R2Point p3(random1.at(1).centerX, random1.at(1).centerY);
    R2Point p4(random2.at(1).centerX, random2.at(1).centerY);
    R2Point p5(random1.at(2).centerX, random1.at(2).centerY);
    R2Point p6(random2.at(2).centerX, random2.at(2).centerY);
    R2Point p7(random1.at(3).centerX, random1.at(3).centerY);
    R2Point p8(random2.at(3).centerX, random2.at(3).centerY);

    R2Point points[] = {p1, p2, p3, p4, p5, p6, p7, p8};
    double** linEquations = dmatrix(1, 8, 1, 9);

    for (int i = 0 ; i < 4; i++){
      linEquations[i*2+1][1] = -points[i*2][0];
      linEquations[i*2+1][2] = -points[i*2][1];
      linEquations[i*2+1][3] = -1.0;
      linEquations[i*2+1][4] = 0.0;
      linEquations[i*2+1][5] = 0.0;
      linEquations[i*2+1][6] = 0.0;
      linEquations[i*2+1][7] = points[i*2][0]*points[i*2+1][0];
      linEquations[i*2+1][8] = points[i*2][1]*points[i*2+1][0];
      linEquations[i*2+1][9] = points[i*2+1][0];

      linEquations[i*2+2][1] = 0.0;
      linEquations[i*2+2][2] = 0.0;
      linEquations[i*2+2][3] = 0.0;
      linEquations[i*2+2][4] = -points[i*2][0];
      linEquations[i*2+2][5] = -points[i*2][1];
      linEquations[i*2+2][6] = -1.0;
      linEquations[i*2+2][7] = points[i*2][0]*points[i*2+1][1];
      linEquations[i*2+2][8] = points[i*2][1]*points[i*2+1][1];
      linEquations[i*2+2][9] = points[i*2+1][1];
    }

    double singularValues[10];
    double** nullspaceMatrix = dmatrix(1, 9, 1, 9);
    svdcmp(linEquations, 8, 9, singularValues, nullspaceMatrix);

    int smallestIndex = 1;
    for(int i=2;i<10;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

    double** homographyMatrix = dmatrix(1, 3, 1, 3);
    homographyMatrix[1][1] = nullspaceMatrix[1][smallestIndex];
    homographyMatrix[1][2] = nullspaceMatrix[2][smallestIndex];
    homographyMatrix[1][3] = nullspaceMatrix[3][smallestIndex];
    homographyMatrix[2][1] = nullspaceMatrix[4][smallestIndex];
    homographyMatrix[2][2] = nullspaceMatrix[5][smallestIndex];
    homographyMatrix[2][3] = nullspaceMatrix[6][smallestIndex];
    homographyMatrix[3][1] = nullspaceMatrix[7][smallestIndex];
    homographyMatrix[3][2] = nullspaceMatrix[8][smallestIndex];
    homographyMatrix[3][3] = nullspaceMatrix[9][smallestIndex];

    // find the number of inliers
    int inlier = 0;
    for (int num = 0; num < feature1.size(); num++){
      int x1 = feature1.at(num).centerX;
      int y1 = feature1.at(num).centerY;
      int x2 = feature2.at(num).centerX;
      int y2 = feature2.at(num).centerY;
      double hpointZ = homographyMatrix[3][1]*x1+homographyMatrix[3][2]*y1+homographyMatrix[3][3];
      double hpointX = (homographyMatrix[1][1]*x1+homographyMatrix[1][2]*y1+homographyMatrix[1][3])/hpointZ;
      double hpointY = (homographyMatrix[2][1]*x1+homographyMatrix[2][2]*y1+homographyMatrix[2][3])/hpointZ;

      if(sqrt(pow(hpointX-x2,2)+pow(hpointY-y2,2)) < 5){
        inlier ++;
      }
    }

    if (trial == 0 || inlier > inlierNum){
      HMatrix = homographyMatrix;
      inlierNum = inlier;
      std::cout << "inlier is " << inlierNum << std::endl;
    }
  }  
  std::vector<R2Point> matchPoints;

  for (int num = 0; num < feature1.size(); num++){
    int x1 = feature1.at(num).centerX;
    int y1 = feature1.at(num).centerY;
    int x2 = feature2.at(num).centerX;
    int y2 = feature2.at(num).centerY;
    double hpointZ = HMatrix[3][1]*x1 + HMatrix[3][2]*y1 + HMatrix[3][3];
    double hpointX = (HMatrix[1][1]*x1+HMatrix[1][2]*y1+HMatrix[1][3])/hpointZ;
    double hpointY = (HMatrix[2][1]*x1+HMatrix[2][2]*y1+HMatrix[2][3])/hpointZ;
    double difference = sqrt(pow(hpointX-x2,2)+pow(hpointY-y2,2));
    if (difference < 5){
      matchPoints.push_back(R2Point(x1,y1));
      matchPoints.push_back(R2Point(x2,y2));
    }

    int size = matchPoints.size();
    std::cout << "size is " << size << std::endl;
    double** linearEquations = dmatrix(1, size, 1, 9);
    for (int i = 0 ; i < (int)(size/2); i++){
      linearEquations[i*2+1][1] = -matchPoints.at(i*2)[0];
      linearEquations[i*2+1][2] = -matchPoints.at(i*2)[1];
      linearEquations[i*2+1][3] = -1.0;
      linearEquations[i*2+1][4] = 0.0;
      linearEquations[i*2+1][5] = 0.0;
      linearEquations[i*2+1][6] = 0.0;
      linearEquations[i*2+1][7] = matchPoints.at(i*2)[0]*matchPoints.at(i*2+1)[0];
      linearEquations[i*2+1][8] = matchPoints.at(i*2)[1]*matchPoints.at(i*2+1)[0];
      linearEquations[i*2+1][9] = matchPoints.at(i*2+1)[0];

      linearEquations[i*2+2][1] = 0.0;
      linearEquations[i*2+2][2] = 0.0;
      linearEquations[i*2+2][3] = 0.0;
      linearEquations[i*2+2][4] = -matchPoints.at(i*2)[0];
      linearEquations[i*2+2][5] = -matchPoints.at(i*2)[1];
      linearEquations[i*2+2][6] = -1.0;
      linearEquations[i*2+2][7] = matchPoints.at(i*2)[0]*matchPoints.at(i*2+1)[1];
      linearEquations[i*2+2][8] = matchPoints.at(i*2)[1]*matchPoints.at(i*2+1)[1];
      linearEquations[i*2+2][9] = matchPoints.at(i*2+1)[1];
      } 

    double singularVs[10];
    double** nullMatrix = dmatrix(1, 9, 1, 9);
    svdcmp(linearEquations, size, 9, singularVs, nullMatrix);
    

    int smallestIndex = 1;
    for(int i=2;i<10;i++) if(singularVs[i]<singularVs[smallestIndex]) smallestIndex=i;
    newHMatrix[1][1] = nullMatrix[1][smallestIndex];
    newHMatrix[1][2] = nullMatrix[2][smallestIndex];
    newHMatrix[1][3] = nullMatrix[3][smallestIndex];
    newHMatrix[2][1] = nullMatrix[4][smallestIndex];
    newHMatrix[2][2] = nullMatrix[5][smallestIndex];
    newHMatrix[2][3] = nullMatrix[6][smallestIndex];
    newHMatrix[3][1] = nullMatrix[7][smallestIndex];
    newHMatrix[3][2] = nullMatrix[8][smallestIndex];
    newHMatrix[3][3] = nullMatrix[9][smallestIndex];
  }

  // //draw the graph
  // for(int x = 0; x < width; x++){
  //   for(int y = 0; y < height; y++){
  //     double img2Z = newHMatrix[3][1]*x+newHMatrix[3][2]*y+newHMatrix[3][3];
  //     double img2X = ((newHMatrix[1][1]*x+newHMatrix[1][2]*y+newHMatrix[1][3])/img2Z);
  //     double img2Y = ((newHMatrix[2][1]*x+newHMatrix[2][2]*y+newHMatrix[2][3])/img2Z);
  //     if(img2X>=0 and img2X < width and img2Y>=0 and img2Y <height){
  //       Pixel(x,y) = 0.5*temp.Pixel(x,y) + 0.5*image2.Pixel(img2X,img2Y);
  //     }
  //   }
  // }

  //draw the feature pairs
  for (int num = 0; num < 150; num++){
    int x1 = feature1.at(num).centerX;
    int y1 = feature1.at(num).centerY;
    int x2 = feature2.at(num).centerX;
    int y2 = feature2.at(num).centerY;
    double hpointZ = HMatrix[3][1]*x1+HMatrix[3][2]*y1+HMatrix[3][3];
    double hpointX = (HMatrix[1][1]*x1+HMatrix[1][2]*y1+HMatrix[1][3])/hpointZ;
    double hpointY = (HMatrix[2][1]*x1+HMatrix[2][2]*y1+HMatrix[2][3])/hpointZ;
    double difference = sqrt(pow(hpointX-x2,2)+pow(hpointY-y2,2));
    std::cout << "difference is " << difference << std::endl;
    if(difference < 5){
      // draw them green
      for(int a = -5; a < 6; a++){
        otherImage->Pixel(x1 + a, y1 - 5) = *greenPixel;
        otherImage->Pixel(x1 + a, y1 + 5) = *greenPixel;
        otherImage->Pixel(x1 - 5, y1 + a) = *greenPixel;
        otherImage->Pixel(x1 + 5, y1 + a) = *greenPixel;
        otherImage->Pixel(x2 + a, y2 - 5) = *greenPixel;
        otherImage->Pixel(x2 + a, y2 + 5) = *greenPixel;
        otherImage->Pixel(x2 - 5, y2 + a) = *greenPixel;
        otherImage->Pixel(x2 + 5, y2 + a) = *greenPixel;
      }
      otherImage->line(x1, x2, y1, y2, 0.0, 1.0, 0.0);
    }else{
      for(int a = -5; a < 6; a++){
        otherImage->Pixel(x1 + a, y1 - 5) = *redPixel;
        otherImage->Pixel(x1 + a, y1 + 5) = *redPixel;
        otherImage->Pixel(x1 - 5, y1 + a) = *redPixel;
        otherImage->Pixel(x1 + 5, y1 + a) = *redPixel;
        otherImage->Pixel(x2 + a, y2 - 5) = *redPixel;
        otherImage->Pixel(x2 + a, y2 + 5) = *redPixel;
        otherImage->Pixel(x2 - 5, y2 + a) = *redPixel;
        otherImage->Pixel(x2 + 5, y2 + a) = *redPixel;
      }
      otherImage->line(x1, x2, y1, y2, 1.0, 0.0, 0.0);
    }
  }




}

////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  npixels = width = height = 0;

  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);
  
  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
  else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
  unsigned char rgbtBlue;
  unsigned char rgbtGreen;
  unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
  unsigned char rgbBlue;
  unsigned char rgbGreen;
  unsigned char rgbRed;
  unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
  // Read a unsigned short int from a file in little endian format 
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
  // Write a unsigned short int to a file in little endian format
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp); 
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
  // Read a unsigned int word from a file in little endian format 
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
  // Write a unsigned int to a file in little endian format 
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
  // Read a int word from a file in little endian format 
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
  // Write a int to a file in little endian format 
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  /* Read file header */
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);
  
  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);
  
  /* Read info header */
  BITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);
  
  // Check info header 
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and number of pixels
  width = bmih.biWidth;
  height = bmih.biHeight;
  npixels = width * height;

  // Allocate unsigned char buffer for reading pixels
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = bmih.biSizeImage;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read buffer 
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Allocate pixels for image
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double b = (double) *(p++) / 255;
      double g = (double) *(p++) / 255;
      double r = (double) *(p++) / 255;
      R2Pixel pixel(r, g, b, 1);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Compute number of bytes in row
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Write file header 
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header 
  BITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      double r = 255.0 * pixel.Red();
      double g = 255.0 * pixel.Green();
      double b = 255.0 * pixel.Blue();
      if (r >= 255) r = 255;
      if (g >= 255) g = 255;
      if (b >= 255) b = 255;
      fputc((unsigned char) b, fp);
      fputc((unsigned char) g, fp);
      fputc((unsigned char) r, fp);
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }
  
  // Close file
  fclose(fp);

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Read PPM file magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }
  
  // Read max value
  double max_value;
  if (fscanf(fp, "%lf", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }
  
  // Allocate image pixels
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check if raw or ascii file
  if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read raw image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
        double r = (double) getc(fp) / max_value;
        double g = (double) getc(fp) / max_value;
        double b = (double) getc(fp) / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }
  else {
    // Read asci image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
  // Read pixel values
  int red, green, blue;
  if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
    fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
    fclose(fp);
    return 0;
  }

  // Assign pixel values
  double r = (double) red / max_value;
  double g = (double) green / max_value;
  double b = (double) blue / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }

    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
        if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s", filename);
      return 0;
    }
    
    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%c%c%c", r, g, b);
      }
    }
    
    // Close file
    fclose(fp);
  }

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
  extern "C" { 
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
  };
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  npixels = width * height;
  int ncomponents = cinfo.output_components;

  // Allocate pixels for image
  pixels = new R2Pixel [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Allocate unsigned char buffer for reading image
  int rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines 
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double r, g, b, a;
      if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
        p++;
      }
      else if (ncomponents == 3) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 4) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = (double) *(p++) / 255;
      }
      else {
        fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
        return 0;
      }
      R2Pixel pixel(r, g, b, a);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}


  

int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width;  /* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;   /* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB;   /* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 95, TRUE);
  jpeg_start_compress(&cinfo, TRUE);
  
  // Allocate unsigned char buffer for reading image
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Fill buffer with pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      int r = (int) (255 * pixel.Red());
      int g = (int) (255 * pixel.Green());
      int b = (int) (255 * pixel.Blue());
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;
      *(p++) = r;
      *(p++) = g;
      *(p++) = b;
    }
  }



  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return number of bytes written
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}






