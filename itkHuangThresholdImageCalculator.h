
#ifndef __itkHuangThresholdImageCalculator_h
#define __itkHuangThresholdImageCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkNumericTraits.h"

namespace itk
{

/** \class HuangThresholdImageCalculator
 * \brief Computes the Huang's threshold for an image.
 * 
 * This calculator computes the Huang's fuzzy threshold which separates an image
 * into foreground and background components. Uses Shannon's entropy
 * function (one can also use Yager's entropy function)  
 * Huang L.-K. and Wang M.-J.J. (1995) "Image Thresholding by Minimizing  
 * the Measures of Fuzziness" Pattern Recognition, 28(1): 41-51
 * Reimplemented (to handle 16-bit efficiently) by Johannes Schindelin Jan 31, 2011

 * Ported from the ImageJ implementation. 
 *
 * This class is templated over the input image type.
 * \author Richard Beare
 * \warning This method assumes that the input image consists of scalar pixel
 * types.
 *
 * \ingroup Operators
 */
template <class TInputImage>
class ITK_EXPORT HuangThresholdImageCalculator : public Object 
{
public:
  /** Standard class typedefs. */
  typedef HuangThresholdImageCalculator Self;
  typedef Object                       Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(HuangThresholdImageCalculator, Object);

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

  /** Compute the Huang's threshold for the input image. */
  void Compute(void);

  /** Return the Huang's threshold value. */
  itkGetConstMacro(Threshold,PixelType);
  
  /** Set/Get the number of histogram bins. Default is 128. */
  itkSetClampMacro( NumberOfHistogramBins, unsigned long, 1, 
                    NumericTraits<unsigned long>::max() );
  itkGetConstMacro( NumberOfHistogramBins, unsigned long );


  itkSetClampMacro(LowThresh, double, 0.0, 1.0);
  itkGetConstMacro(LowThresh, double);

  itkSetClampMacro(HighThresh, double, 0.0, 1.0);
  itkGetConstMacro(HighThresh, double);

  /** Set the region over which the values will be computed */
  void SetRegion( const RegionType & region );

protected:
  HuangThresholdImageCalculator();
  virtual ~HuangThresholdImageCalculator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  HuangThresholdImageCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  PixelType            m_Threshold;
  double               m_LowThresh;
  double               m_HighThresh;
  unsigned long        m_NumberOfHistogramBins;
  ImageConstPointer    m_Image;
  RegionType           m_Region;
  bool                 m_RegionSetByUser;

};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHuangThresholdImageCalculator.txx"
#endif

#endif
