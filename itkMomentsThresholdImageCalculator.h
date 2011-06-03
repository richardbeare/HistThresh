
#ifndef __itkMomentsThresholdImageCalculator_h
#define __itkMomentsThresholdImageCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkNumericTraits.h"

namespace itk
{

/** \class MomentsThresholdImageCalculator
 * \brief Computes the Moments's threshold for an image.
 * 
 *  W. Tsai, "Moment-preserving thresholding: a new approach," Computer Vision,
 * Graphics, and Image Processing, vol. 29, pp. 377-393, 1985.
 * Ported to ImageJ plugin by G.Landini from the the open source project FOURIER 0.8
 * by  M. Emre Celebi , Department of Computer Science,  Louisiana State University in Shreveport
 * Shreveport, LA 71115, USA
 *  http://sourceforge.net/projects/fourier-ipal
 *  http://www.lsus.edu/faculty/~ecelebi/fourier.htm
 *
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
class ITK_EXPORT MomentsThresholdImageCalculator : public Object 
{
public:
  /** Standard class typedefs. */
  typedef MomentsThresholdImageCalculator Self;
  typedef Object                       Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MomentsThresholdImageCalculator, Object);

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

  /** Compute the Moments's threshold for the input image. */
  void Compute(void);

  /** Return the Moments's threshold value. */
  itkGetConstMacro(Threshold,PixelType);
  
  /** Set/Get the number of histogram bins. Default is 128. */
  itkSetClampMacro( NumberOfHistogramBins, unsigned long, 1, 
                    NumericTraits<unsigned long>::max() );
  itkGetConstMacro( NumberOfHistogramBins, unsigned long );

  /** Set the region over which the values will be computed */
  void SetRegion( const RegionType & region );

protected:
  MomentsThresholdImageCalculator();
  virtual ~MomentsThresholdImageCalculator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  MomentsThresholdImageCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  PixelType            m_Threshold;
  unsigned long        m_NumberOfHistogramBins;
  ImageConstPointer    m_Image;
  RegionType           m_Region;
  bool                 m_RegionSetByUser;

};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMomentsThresholdImageCalculator.txx"
#endif

#endif
