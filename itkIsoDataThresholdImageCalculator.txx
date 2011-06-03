
#ifndef __itkIsoDataThresholdImageCalculator_txx
#define __itkIsoDataThresholdImageCalculator_txx

#include "itkIsoDataThresholdImageCalculator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "vnl/vnl_math.h"

namespace itk
{ 
    
/**
 * Constructor
 */
template<class TInputImage>
IsoDataThresholdImageCalculator<TInputImage>
::IsoDataThresholdImageCalculator()
{
  m_Image = NULL;
  m_Threshold = NumericTraits<PixelType>::Zero;
  m_NumberOfHistogramBins = 128;
  m_RegionSetByUser = false;
  m_LowThresh=0.01;
  m_HighThresh=0.99;
}


/*
 * Compute the IsoData's threshold
 */
template<class TInputImage>
void
IsoDataThresholdImageCalculator<TInputImage>
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
  std::vector<double> cumSum;
  std::vector<double> triangle;
  relativeFrequency.resize( m_NumberOfHistogramBins );
  cumSum.resize( m_NumberOfHistogramBins );
  triangle.resize( m_NumberOfHistogramBins );

  std::fill(relativeFrequency.begin(), relativeFrequency.end(), 0.0);
  std::fill(cumSum.begin(), cumSum.end(), 0.0);
  std::fill(triangle.begin(), triangle.end(), 0.0);

  // for ( j = 0; j < m_NumberOfHistogramBins; j++ )
  //   {
  //   relativeFrequency[j] = 0.0;
  //   cumSum[j]=0.0;
  //   }

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

  int i, l, toth, totl, h, g=0;
  for (i = 1; (unsigned)i < relativeFrequency.size(); i++)
    {
    if (relativeFrequency[i] > 0)
      {
      g = i + 1;
      break;
      }
    }
  while (true)
    {
    l = 0;
    totl = 0;
    for (i = 0; i < g; i++) 
      {
      totl = totl + relativeFrequency[i];
      l = l + (relativeFrequency[i] * i);
      }
    h = 0;
    toth = 0;
    for (i = g + 1; (unsigned)i < relativeFrequency.size(); i++)
      {
      toth += relativeFrequency[i];
      h += (relativeFrequency[i]*i);
      }
    if (totl > 0 && toth > 0)
      {
      l /= totl;
      h /= toth;
      if (g == (int) round((l + h) / 2.0))
	break;
      }
    g++;
    if ((unsigned)g >relativeFrequency.size()-2)
      {
      itkWarningMacro(<<"IsoData Threshold not found.");
      return;
      }
    }


  m_Threshold = static_cast<PixelType>( imageMin + 
                                        ( g ) / binMultiplier );


}

template<class TInputImage>
void
IsoDataThresholdImageCalculator<TInputImage>
::SetRegion( const RegionType & region )
{
  m_Region = region;
  m_RegionSetByUser = true;
}

  
template<class TInputImage>
void
IsoDataThresholdImageCalculator<TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Threshold: " << m_Threshold << std::endl;
  os << indent << "NumberOfHistogramBins: " << m_NumberOfHistogramBins << std::endl;
  os << indent << "Image: " << m_Image.GetPointer() << std::endl;
}

} // end namespace itk

#endif
