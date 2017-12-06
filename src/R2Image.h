// Include file for image class
#ifndef R2_IMAGE_INCLUDED
#define R2_IMAGE_INCLUDED

#include "vector"



// Constant definitions

typedef enum {
  R2_IMAGE_RED_CHANNEL,
  R2_IMAGE_GREEN_CHANNEL,
  R2_IMAGE_BLUE_CHANNEL,
  R2_IMAGE_ALPHA_CHANNEL,
  R2_IMAGE_NUM_CHANNELS
} R2ImageChannel;

typedef enum {
  R2_IMAGE_POINT_SAMPLING,
  R2_IMAGE_BILINEAR_SAMPLING,
  R2_IMAGE_GAUSSIAN_SAMPLING,
  R2_IMAGE_NUM_SAMPLING_METHODS
} R2ImageSamplingMethod;

typedef enum {
  R2_IMAGE_OVER_COMPOSITION,
  R2_IMAGE_IN_COMPOSITION,
  R2_IMAGE_OUT_COMPOSITION,
  R2_IMAGE_ATOP_COMPOSITION,
  R2_IMAGE_XOR_COMPOSITION,
} R2ImageCompositeOperation;



// Class definition

class R2Image {
 public:
  // Constructors/destructor
  R2Image(void);
  R2Image(const char *filename);
  R2Image(int width, int height);
  R2Image(int width, int height, const R2Pixel *pixels);
  R2Image(const R2Image& image);
  ~R2Image(void);

  // Image properties
  int NPixels(void) const;
  int Width(void) const;
  int Height(void) const;

  // Pixel access/update
  R2Pixel& Pixel(int x, int y);
  R2Pixel *Pixels(void);
  R2Pixel *Pixels(int row);
  R2Pixel *operator[](int row);
  const R2Pixel *operator[](int row) const;
  void SetPixel(int x, int y,  const R2Pixel& pixel);

  // additional functionality
  void line(int x0, int x1, int y0, int y1, float r, float g, float b);
  void square(int a0, int a1, int b0, int b1, int c0, int c1, int d0, int d1); 

  public:
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

  // Image processing
  R2Image& operator=(const R2Image& image);

  // Per-pixel operations
  void Brighten(double factor);
  void ChangeSaturation(double factor);

  // show how SVD works
  void svdTest();

  // Linear filtering operations
  void SobelX();
  void SobelY();
  void LoG();
  void Blur(double sigma);
  void Harris(double sigma);
  void Sharpen(void);
  void FeatureDetection();


  // further operations
  void blendOtherImageTranslated(R2Image * otherImage);
  void blendOtherImageHomography(R2Image * otherImage);

  void computeHomography(R2Point p1, R2Point p2, R2Point p3, R2Point p4, R2Point p5, R2Point p6, R2Point p7, R2Point p8);


  // Non Linear filtering operations
  void Bilateral(double sigma, double tolerance);
  void Median(int kernelSize, double sigma);
  void LensDistortion();
  void ScaleInvariantHarris();

  // video operations
  void FirstFrameProcessing();
  void FrameProcessing(R2Image * prevImage, R2Image * currentImage, std::vector<Feature> temp);

  // File reading/writing
  int Read(const char *filename);
  int ReadBMP(const char *filename);
  int ReadPPM(const char *filename);
  int ReadJPEG(const char *filename);
  int Write(const char *filename) const;
  int WriteBMP(const char *filename) const;
  int WritePPM(const char *filename, int ascii = 0) const;
  int WriteJPEG(const char *filename) const;

 private:
  // Utility functions
  void Resize(int width, int height);
  R2Pixel Sample(double u, double v,  int sampling_method);

 private:
  R2Pixel *pixels;
  int npixels;
  int width;
  int height;

  public:
    std::vector<Feature> prevStoredFeature;
    std::vector<Feature> currStoredFeature;

};

// Inline functions

inline int R2Image::
NPixels(void) const
{
  // Return total number of pixels
  return npixels;
}



inline int R2Image::
Width(void) const
{
  // Return width
  return width;
}



inline int R2Image::
Height(void) const
{
  // Return height
  return height;
}



inline R2Pixel& R2Image::
Pixel(int x, int y)
{
  // Return pixel value at (x,y)
  // (pixels start at lower-left and go in row-major order)
  return pixels[x*height + y];
}



inline R2Pixel *R2Image::
Pixels(void)
{
  // Return pointer to pixels for whole image 
  // (pixels start at lower-left and go in row-major order)
  return pixels;
}



inline R2Pixel *R2Image::
Pixels(int x)
{
  // Return pixels pointer for row at x
  // (pixels start at lower-left and go in row-major order)
  return &pixels[x*height];
}



inline R2Pixel *R2Image::
operator[](int x) 
{
  // Return pixels pointer for row at x
  return Pixels(x);
}



inline const R2Pixel *R2Image::
operator[](int x) const
{
  // Return pixels pointer for row at x
  // (pixels start at lower-left and go in row-major order)
  return &pixels[x*height];
}



inline void R2Image::
SetPixel(int x, int y, const R2Pixel& pixel)
{
  // Set pixel
  pixels[x*height + y] = pixel;
}



#endif
