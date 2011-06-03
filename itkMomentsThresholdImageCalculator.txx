
#ifndef __itkMomentsThresholdImageCalculator_txx
#define __itkMomentsThresholdImageCalculator_txx

#include "itkMomentsThresholdImageCalculator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "vnl/vnl_math.h"

namespace itk
{ 
    
/**
 * Constructor
 */
template<class TInputImage>
MomentsThresholdImageCalculator<TInputImage>
::MomentsThresholdImageCalculator()
{
  m_Image = NULL;
  m_Threshold = NumericTraits<PixelType>::Zero;
  m_NumberOfHistogramBins = 128;
  m_RegionSetByUser = false;
}


/*
 * Compute the Moments's threshold
 */
template<class TInputImage>
void
MomentsThresholdImageCalculator<TInputImage>
::Compute(void)
{

  if ( !m_Image ) { return; }
  if( !m_RegionSetByUser )
    {
    m_Region = m_Image->GetRequestedRegion();
    }

  double totalPixels = (double) m_Region.GetNumberOfPixels();
  if ( totalPixels == 0 ) { return; }


  // compute image max and min
  typedef MinimumMaximumImageCalculator<TInputImage> RangeCalculator;
  typename RangeCalculator::Pointer rangeCalculator = RangeCalculator::New();
  rangeCalculator->SetImage( m_Image );
  rangeCalculator->Compute();

  PixelType imageMin = rangeCalculator->GetMinimum();
  PixelType imageMax = rangeCalculator->GetMaximum();

  if ( imageMin >= imageMax )
    {
    m_Threshold = imageMin;
    return;
    }

  // create a histogram
  std::vector<double> relativeFrequency;
  relativeFrequency.resize( m_NumberOfHistogramBins );
  std::fill(relativeFrequency.begin(), relativeFrequency.end(), 0.0);


  double binMultiplier = (double) m_NumberOfHistogramBins /
    (double) ( imageMax - imageMin );

  typedef ImageRegionConstIteratorWithIndex<TInputImage> Iterator;
  Iterator iter( m_Image, m_Region );

  while ( !iter.IsAtEnd() )
    {
    unsigned int binNumber;
    PixelType value = iter.Get();

    if ( value == imageMin ) 
      {
      binNumber = 0;
      }
    else
      {
      binNumber = (unsigned int) vcl_ceil((value - imageMin) * binMultiplier ) - 1;
      if ( binNumber == m_NumberOfHistogramBins ) // in case of rounding errors
        {
        binNumber -= 1;
        }
      }

    relativeFrequency[binNumber] += 1.0;
    ++iter;

    }
  double total =0;
  double m0=1.0, m1=0.0, m2 =0.0, m3 =0.0, sum =0.0, p0=0.0;
  double cd, c0, c1, z0, z1;	/* auxiliary variables */
  int threshold = -1;
  
  std::vector<double> histo(relativeFrequency.size());
  
  for (unsigned i=0; i<relativeFrequency.size(); i++)
    total+=relativeFrequency[i];
  
  for (unsigned i=0; i<relativeFrequency.size(); i++)
    histo[i]=(double)(relativeFrequency[i]/total); //normalised histogram
  
  /* Calculate the first, second, and third order moments */
  for ( unsigned i = 0; i < relativeFrequency.size(); i++ )
    {
    m1 += i * histo[i];
    m2 += i * i * histo[i];
    m3 += i * i * i * histo[i];
    }
  // 
  // First 4 moments of the gray-level image should match the first 4 moments
  // of the target binary image. This leads to 4 equalities whose solutions 
  // are given in the Appendix of Ref. 1 
  //
  cd = m0 * m2 - m1 * m1;
  c0 = ( -m2 * m2 + m1 * m3 ) / cd;
  c1 = ( m0 * -m3 + m2 * m1 ) / cd;
  z0 = 0.5 * ( -c1 - vcl_sqrt ( c1 * c1 - 4.0 * c0 ) );
  z1 = 0.5 * ( -c1 + vcl_sqrt ( c1 * c1 - 4.0 * c0 ) );
  p0 = ( z1 - m1 ) / ( z1 - z0 );  /* Fraction of the object pixels in the target binary image */

  // The threshold is the gray-level closest  
  // to the p0-tile of the normalized histogram 
  sum=0;
  for (unsigned i=0; i<relativeFrequency.size(); i++)
    {
    sum+=histo[i];
    if (sum>p0) 
      {
      threshold = i;
      break;
      }
    }

    m_Threshold = static_cast<PixelType>( imageMin + 
					  ( threshold ) / binMultiplier );


}

template<class TInputImage>
void
MomentsThresholdImageCalculator<TInputImage>
::SetRegion( const RegionType & region )
{
  m_Region = region;
  m_RegionSetByUser = true;
}

  
template<class TInputImage>
void
MomentsThresholdImageCalculator<TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Threshold: " << m_Threshold << std::endl;
  os << indent << "NumberOfHistogramBins: " << m_NumberOfHistogramBins << std::endl;
  os << indent << "Image: " << m_Image.GetPointer() << std::endl;
}

} // end namespace itk

#endif
