
#ifndef __itkLiThresholdImageCalculator_h
#define __itkLiThresholdImageCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkNumericTraits.h"

namespace itk
{

/** \class LiThresholdImageCalculator
 * \brief Computes the Li's threshold for an image.
 * 
 * Ported from the ImageJ implementation. 
 * Implements Li's Minimum Cross Entropy thresholding method
 * This implementation is based on the iterative version (Ref. 2) of the algorithm.
 * 1) Li C.H. and Lee C.K. (1993) "Minimum Cross Entropy Thresholding" 
 *    Pattern Recognition, 26(4): 617-625
 * 2) Li C.H. and Tam P.K.S. (1998) "An Iterative Algorithm for Minimum 
 *    Cross Entropy Thresholding"Pattern Recognition Letters, 18(8): 771-776
 * 3) Sezgin M. and Sankur B. (2004) "Survey over Image Thresholding 
 *    Techniques and Quantitative Performance Evaluation" Journal of 
 *    Electronic Imaging, 13(1): 146-165 
 *    http://citeseer.ist.psu.edu/sezgin04survey.html
 * Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines
 *
 *
 * This class is templated over the input image type.
 * \author Richard Beare
 * \warning This method assumes that the input image consists of scalar pixel
 * types.
 *
 * \ingroup Operators
 */
template <class TInputImage>
class ITK_EXPORT LiThresholdImageCalculator : public Object 
{
public:
  /** Standard class typedefs. */
  typedef LiThresholdImageCalculator Self;
  typedef Object                       Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(LiThresholdImageCalculator, Object);

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

  /** Compute the Li's threshold for the input image. */
  void Compute(void);

  /** Return the Li's threshold value. */
  itkGetConstMacro(Threshold,PixelType);
  
  /** Set/Get the number of histogram bins. Default is 128. */
  itkSetClampMacro( NumberOfHistogramBins, unsigned long, 1, 
                    NumericTraits<unsigned long>::max() );
  itkGetConstMacro( NumberOfHistogramBins, unsigned long );


  /** Set the region over which the values will be computed */
  void SetRegion( const RegionType & region );

protected:
  LiThresholdImageCalculator();
  virtual ~LiThresholdImageCalculator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  LiThresholdImageCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  PixelType            m_Threshold;
  unsigned long        m_NumberOfHistogramBins;
  ImageConstPointer    m_Image;
  RegionType           m_Region;
  bool                 m_RegionSetByUser;

};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLiThresholdImageCalculator.txx"
#endif

#endif
