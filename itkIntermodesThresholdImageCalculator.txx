#ifndef __itkIntermodesThresholdImageCalculator_txx
#define __itkIntermodesThresholdImageCalculator_txx

#include "itkIntermodesThresholdImageCalculator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "vnl/vnl_math.h"

namespace itk
{ 
    
/**
 * Constructor
 */
template<class TInputImage>
IntermodesThresholdImageCalculator<TInputImage>
::IntermodesThresholdImageCalculator()
{
  m_Image = NULL;
  m_Threshold = NumericTraits<PixelType>::Zero;
  m_NumberOfHistogramBins = 128;
  m_RegionSetByUser = false;
  m_MaxSmoothingIterations = 10000;
  m_UseInterMode = true;
}
template<class TInputImage>
bool
IntermodesThresholdImageCalculator<TInputImage>
::bimodalTest(std::vector<double> h)
{
  int modes = 0;
  
  const unsigned len = h.size();
  for (unsigned k = 1; k < len - 1 ; k++)
    {
    if (h[k-1] < h[k] && h[k+1] < h[k]) 
      {
      modes++;
      if (modes>2)
	return false;
      }

    }

  return (modes == 2);
}

/*
 * Compute the Intermodes's threshold
 */
template<class TInputImage>
void
IntermodesThresholdImageCalculator<TInputImage>
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

  // smooth the histogram
  std::vector<double> smoothedHist;
  smoothedHist.resize(m_NumberOfHistogramBins);
  std::copy(relativeFrequency.begin(), relativeFrequency.end(), smoothedHist.begin());

  unsigned SmIter = 0;

  while (!bimodalTest(smoothedHist))
    {
    // smooth with a 3 point running mean
    double previous = 0, current = 0, next = smoothedHist[0];
    for (unsigned i = 0; i < smoothedHist.size() - 1; i++) 
      {
      previous = current;
      current = next;
      next = smoothedHist[i + 1];
      smoothedHist[i] = (previous + current + next) / 3;
      }
    smoothedHist[smoothedHist.size() - 1] = (current + next) / 3;
    SmIter++;
    if (SmIter > m_MaxSmoothingIterations )
      {
      m_Threshold = -1;
      itkWarningMacro( << "Exceeded maximum iterations for histogram smoothing." );
      return;
      }
    }
  if (m_UseInterMode)
    {
    // The threshold is the mean between the two peaks.
    unsigned tt=0;
    for (unsigned i=1; i<smoothedHist.size() - 1; i++) 
      {
      if (smoothedHist[i-1] < smoothedHist[i] && smoothedHist[i+1] < smoothedHist[i])
	{
	tt += i;
	}
      }
    
    m_Threshold = static_cast<PixelType>( imageMin + 
					  ( tt/2.0 ) / binMultiplier );
    }
  else
    {
    // The threshold is the mean between peaks
    unsigned firstpeak=0, lastpeak=0;
    for (unsigned i=1; i<smoothedHist.size() - 1; i++) 
      {
      if (smoothedHist[i-1] < smoothedHist[i] && smoothedHist[i+1] < smoothedHist[i])
	{
	firstpeak = i;
	break;
	}
      }
    double MinVal = smoothedHist[firstpeak];
    unsigned MinPos=firstpeak;

    for (unsigned i=firstpeak + 1; i<smoothedHist.size() - 1; i++) 
      {
      if (smoothedHist[i] < MinVal)
	{
	MinVal = smoothedHist[i];
	MinPos = i;
	}
      if (smoothedHist[i-1] < smoothedHist[i] && smoothedHist[i+1] < smoothedHist[i])
	{
	lastpeak = i;
	break;
	}
      }
    m_Threshold = static_cast<PixelType>( imageMin + 
					  ( MinPos) / binMultiplier );

    }

}

template<class TInputImage>
void
IntermodesThresholdImageCalculator<TInputImage>
::SetRegion( const RegionType & region )
{
  m_Region = region;
  m_RegionSetByUser = true;
}

  
template<class TInputImage>
void
IntermodesThresholdImageCalculator<TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Threshold: " << m_Threshold << std::endl;
  os << indent << "NumberOfHistogramBins: " << m_NumberOfHistogramBins << std::endl;
  os << indent << "Image: " << m_Image.GetPointer() << std::endl;
}

} // end namespace itk

#endif
