
#ifndef __itkShanbhagThresholdImageCalculator_h
#define __itkShanbhagThresholdImageCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkNumericTraits.h"

namespace itk
{

/** \class ShanbhagThresholdImageCalculator
 * \brief Computes the Shanbhag's threshold for an image.
 *
 * Shanhbag A.G. (1994) "Utilization of Information Measure as a Means of
 * Image Thresholding" Graphical Models and Image Processing, 56(5): 414-419
 * Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines
 * Ported from the ImageJ implementation. http://pacific.mpi-cbg.de/wiki/index.php/Auto_Threshold
 *
 * This class is templated over the input image type.
 * \author Richard Beare
 * \warning This method assumes that the input image consists of scalar pixel
 * types.
 *
 * \ingroup Operators
 */
template <class TInputImage>
class ITK_EXPORT ShanbhagThresholdImageCalculator : public Object 
{
public:
  /** Standard class typedefs. */
  typedef ShanbhagThresholdImageCalculator Self;
  typedef Object                       Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ShanbhagThresholdImageCalculator, Object);

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

  /** Compute the Shanbhag's threshold for the input image. */
  void Compute(void);

  /** Return the Shanbhag's threshold value. */
  itkGetConstMacro(Threshold,PixelType);
  
  /** Set/Get the number of histogram bins. Default is 128. */
  itkSetClampMacro( NumberOfHistogramBins, unsigned long, 1, 
                    NumericTraits<unsigned long>::max() );
  itkGetConstMacro( NumberOfHistogramBins, unsigned long );


  /** Set the region over which the values will be computed */
  void SetRegion( const RegionType & region );

protected:
  ShanbhagThresholdImageCalculator();
  virtual ~ShanbhagThresholdImageCalculator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  ShanbhagThresholdImageCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  PixelType            m_Threshold;
  unsigned long        m_NumberOfHistogramBins;
  ImageConstPointer    m_Image;
  RegionType           m_Region;
  bool                 m_RegionSetByUser;

};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkShanbhagThresholdImageCalculator.txx"
#endif

#endif
