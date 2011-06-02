
#ifndef __itkIntermodesThresholdImageCalculator_h
#define __itkIntermodesThresholdImageCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkNumericTraits.h"

namespace itk
{

/** \class IntermodesThresholdImageCalculator
 * \brief Computes the Intermodes's threshold for an image.
 *
 * J. M. S. Prewitt and M. L. Mendelsohn, "The analysis of cell images," in
 * Annals of the New York Academy of Sciences, vol. 128, pp. 1035-1053, 1966.
 *  *
 * Assumes a bimodal histogram. The histogram needs is smoothed (using a
 * running average of size 3, iteratively) until there are only two local maxima.
 * j and k
 * Threshold t is (j+k)/2.
 * Images with histograms having extremely unequal peaks or a broad and
 * ﬂat valley are unsuitable for this method.
 *
 * Ported to ImageJ plugin by G.Landini from Antti Niemisto's Matlab code (GPL)
 * Original Matlab code Copyright (C) 2004 Antti Niemisto
 * See http://www.cs.tut.fi/~ant/histthresh/ for an excellent slide presentation
 * and the original Matlab code.
 *
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
class ITK_EXPORT IntermodesThresholdImageCalculator : public Object 
{
public:
  /** Standard class typedefs. */
  typedef IntermodesThresholdImageCalculator Self;
  typedef Object                       Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(IntermodesThresholdImageCalculator, Object);

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

  /** Compute the Intermodes's threshold for the input image. */
  void Compute(void);

  /** Return the Intermodes's threshold value. */
  itkGetConstMacro(Threshold,PixelType);
  
  /** Set/Get the number of histogram bins. Default is 128. */
  itkSetClampMacro( NumberOfHistogramBins, unsigned long, 1, 
                    NumericTraits<unsigned long>::max() );
  itkGetConstMacro( NumberOfHistogramBins, unsigned long );

  
  itkSetMacro( MaxSmoothingIterations, unsigned long);
  itkGetConstMacro( MaxSmoothingIterations, unsigned long );



  /** Set the region over which the values will be computed */
  void SetRegion( const RegionType & region );

protected:
  IntermodesThresholdImageCalculator();
  virtual ~IntermodesThresholdImageCalculator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  bool bimodalTest(std::vector<double> h);

private:
  IntermodesThresholdImageCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  PixelType            m_Threshold;
  unsigned             m_MaxSmoothingIterations;
  unsigned long        m_NumberOfHistogramBins;
  ImageConstPointer    m_Image;
  RegionType           m_Region;
  bool                 m_RegionSetByUser;

};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkIntermodesThresholdImageCalculator.txx"
#endif

#endif
