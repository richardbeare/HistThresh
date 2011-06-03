
#ifndef __itkYenThresholdImageCalculator_txx
#define __itkYenThresholdImageCalculator_txx

#include "itkYenThresholdImageCalculator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "vnl/vnl_math.h"

namespace itk
{ 
    
/**
 * Constructor
 */
template<class TInputImage>
YenThresholdImageCalculator<TInputImage>
::YenThresholdImageCalculator()
{
  m_Image = NULL;
  m_Threshold = NumericTraits<PixelType>::Zero;
  m_NumberOfHistogramBins = 128;
  m_RegionSetByUser = false;
}


/*
 * Compute the Yen's threshold
 */
template<class TInputImage>
void
YenThresholdImageCalculator<TInputImage>
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
  int threshold;
  int ih, it;
  double crit;
  double max_crit;
  std::vector<double> norm_histo(relativeFrequency.size()); /* normalized histogram */
  std::vector<double> P1(relativeFrequency.size()); /* cumulative normalized histogram */
  std::vector<double> P1_sq(relativeFrequency.size());
  std::vector<double> P2_sq(relativeFrequency.size());
  
  int total =0;
  for (ih = 0; (unsigned)ih < relativeFrequency.size(); ih++ )
    total+=relativeFrequency[ih];
  
  for (ih = 0; (unsigned)ih < relativeFrequency.size(); ih++ )
    norm_histo[ih] = (double)relativeFrequency[ih]/total;
  
  P1[0]=norm_histo[0];
  for (ih = 1; (unsigned)ih < relativeFrequency.size(); ih++ )
    P1[ih]= P1[ih-1] + norm_histo[ih];
  
  P1_sq[0]=norm_histo[0]*norm_histo[0];
  for (ih = 1; (unsigned)ih < relativeFrequency.size(); ih++ )
    P1_sq[ih]= P1_sq[ih-1] + norm_histo[ih] * norm_histo[ih];
  
  P2_sq[relativeFrequency.size() - 1] = 0.0;
  for ( ih = (unsigned)relativeFrequency.size()-2; ih >= 0; ih-- )
    P2_sq[ih] = P2_sq[ih + 1] + norm_histo[ih + 1] * norm_histo[ih + 1];
  
  /* Find the threshold that maximizes the criterion */
  threshold = -1;
  max_crit = itk::NumericTraits<double>::NonpositiveMin();
  for ( it = 0; (unsigned)it < relativeFrequency.size(); it++ ) 
    {
    crit = -1.0 * (( P1_sq[it] * P2_sq[it] )> 0.0? vcl_log( P1_sq[it] * P2_sq[it]):0.0) +  2 * ( ( P1[it] * ( 1.0 - P1[it] ) )>0.0? vcl_log(  P1[it] * ( 1.0 - P1[it] ) ): 0.0);
    if ( crit > max_crit ) 
      {
      max_crit = crit;
      threshold = it;
      }
    }

  m_Threshold = static_cast<PixelType>( imageMin + 
                                        ( threshold ) / binMultiplier );


}

template<class TInputImage>
void
YenThresholdImageCalculator<TInputImage>
::SetRegion( const RegionType & region )
{
  m_Region = region;
  m_RegionSetByUser = true;
}

  
template<class TInputImage>
void
YenThresholdImageCalculator<TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Threshold: " << m_Threshold << std::endl;
  os << indent << "NumberOfHistogramBins: " << m_NumberOfHistogramBins << std::endl;
  os << indent << "Image: " << m_Image.GetPointer() << std::endl;
}

} // end namespace itk

#endif
