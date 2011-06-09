#ifndef __itkTriangleThresholdImageCalculator_h
#define __itkTriangleThresholdImageCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkNumericTraits.h"

namespace itk
{

/** \class TriangleThresholdImageCalculator
 * \brief Computes the Triangle's threshold for an image.
 * 
 * This calculator computes the Triangle's threshold which separates an image
 * into foreground and background components. The method relies on a
 * histogram of image intensities. A line is drawn between the peak
 * point in the hist and the furthest zero point (robustly estimated
 * as the 1% or 99% point). The threshold is the position of maximum
 * difference between the line and the original histogram.
 *
 * This class is templated over the input image type.
 *
 * \warning This method assumes that the input image consists of scalar pixel
 * types.
 *
 * \ingroup Operators
 */
template <class TInputImage>
class ITK_EXPORT TriangleThresholdImageCalculator : public Object 
{
public:
  /** Standard class typedefs. */
  typedef TriangleThresholdImageCalculator Self;
  typedef Object                       Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(TriangleThresholdImageCalculator, Object);

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

  /** Compute the Triangle's threshold for the input image. */
  void Compute(void);

  /** Return the Triangle's threshold value. */
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
  TriangleThresholdImageCalculator();
  virtual ~TriangleThresholdImageCalculator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  TriangleThresholdImageCalculator(const Self&); //purposely not implemented
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
#include "itkTriangleThresholdImageCalculator.txx"
#endif

#endif
