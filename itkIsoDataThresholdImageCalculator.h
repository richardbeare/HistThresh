
#ifndef __itkIsoDataThresholdImageCalculator_h
#define __itkIsoDataThresholdImageCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkNumericTraits.h"

namespace itk
{

/** \class IsoDataThresholdImageCalculator
 * \brief Computes the IsoData threshold for an image. Aka intermeans
 *
 * Iterative procedure based on the isodata algorithm [T.W. Ridler, S. Calvard, Picture 
 * thresholding using an iterative selection method, IEEE Trans. System, Man and 
 * Cybernetics, SMC-8 (1978) 630-632.] 
 * The procedure divides the image into objects and background by taking an initial threshold,
 * then the averages of the pixels at or below the threshold and pixels above are computed. 
 * The averages of those two values are computed, the threshold is incremented and the 
 * process is repeated until the threshold is larger than the composite average. That is,
 * threshold = (average background + average objects)/2
 *
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
class ITK_EXPORT IsoDataThresholdImageCalculator : public Object 
{
public:
  /** Standard class typedefs. */
  typedef IsoDataThresholdImageCalculator Self;
  typedef Object                       Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(IsoDataThresholdImageCalculator, Object);

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

  /** Compute the IsoData's threshold for the input image. */
  void Compute(void);

  /** Return the IsoData's threshold value. */
  itkGetConstMacro(Threshold,PixelType);
  
  /** Set/Get the number of histogram bins. Default is 128. */
  itkSetClampMacro( NumberOfHistogramBins, unsigned long, 1, 
                    NumericTraits<unsigned long>::max() );
  itkGetConstMacro( NumberOfHistogramBins, unsigned long );


  /** Set the region over which the values will be computed */
  void SetRegion( const RegionType & region );

protected:
  IsoDataThresholdImageCalculator();
  virtual ~IsoDataThresholdImageCalculator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  IsoDataThresholdImageCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  PixelType            m_Threshold;
  unsigned long        m_NumberOfHistogramBins;
  ImageConstPointer    m_Image;
  RegionType           m_Region;
  bool                 m_RegionSetByUser;

};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkIsoDataThresholdImageCalculator.txx"
#endif

#endif
