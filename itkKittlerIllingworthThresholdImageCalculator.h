
#ifndef __itkKittlerIllingworthThresholdImageCalculator_h
#define __itkKittlerIllingworthThresholdImageCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkNumericTraits.h"

namespace itk
{

/** \class KittlerIllingworthThresholdImageCalculator
 * \brief Computes the KittlerIllingworth's threshold for an image.
 *
 * Kittler and J. Illingworth, "Minimum error thresholding," Pattern Recognition, vol. 19, pp. 41-47, 1986.
 * C. A. Glasbey, "An analysis of histogram-based thresholding algorithms," CVGIP: Graphical Models and Image Processing, vol. 55, pp. 532-537, 1993.
 * Ported to ImageJ plugin by G.Landini from Antti Niemisto's Matlab code (GPL)
 * Original Matlab code Copyright (C) 2004 Antti Niemisto
 * See http://www.cs.tut.fi/~ant/histthresh/ for an excellent slide presentation
 * and the original Matlab code.
 *
 * Ported to ITK from ImageJ version.
 *
 * \author Richard Beare.
 * This class is templated over the input image type.
 *
 * \warning This method assumes that the input image consists of scalar pixel
 * types.
 *
 * \ingroup Operators
 */
template <class TInputImage>
class ITK_EXPORT KittlerIllingworthThresholdImageCalculator : public Object 
{
public:
  /** Standard class typedefs. */
  typedef KittlerIllingworthThresholdImageCalculator Self;
  typedef Object                       Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(KittlerIllingworthThresholdImageCalculator, Object);

  /** Type definition for the input image. */
  typedef TInputImage  ImageType;

  /** Pointer type for the image. */
  typedef typename TInputImage::Pointer  ImagePointer;
  
  /** Const Pointer type for the image. */
  typedef typename TInputImage::ConstPointer ImageConstPointer;

  /** Type definition for the input image pixel type. */
  typedef typename TInputImage::PixelType PixelType;
  
  /** Type definition for the input image region type. */
  typedef typename TInputImage::RegionType RegionType;
  
  /** Set the input image. */
  itkSetConstObjectMacro(Image,ImageType);

  /** Compute the KittlerIllingworth's threshold for the input image. */
  void Compute(void);

  /** Return the KittlerIllingworth's threshold value. */
  itkGetConstMacro(Threshold,PixelType);
  
  /** Set/Get the number of histogram bins. Default is 128. */
  itkSetClampMacro( NumberOfHistogramBins, unsigned long, 1, 
                    NumericTraits<unsigned long>::max() );
  itkGetConstMacro( NumberOfHistogramBins, unsigned long );

  /** Set the region over which the values will be computed */
  void SetRegion( const RegionType & region );

protected:
  KittlerIllingworthThresholdImageCalculator();
  virtual ~KittlerIllingworthThresholdImageCalculator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  int Mean(std::vector<double> data);
  double A(std::vector<double> y, int j);
  double B(std::vector<double> y, int j);
  double C(std::vector<double> y, int j);

private:
  KittlerIllingworthThresholdImageCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  PixelType            m_Threshold;
  unsigned long        m_NumberOfHistogramBins;
  ImageConstPointer    m_Image;
  RegionType           m_Region;
  bool                 m_RegionSetByUser;

};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkKittlerIllingworthThresholdImageCalculator.txx"
#endif

#endif
