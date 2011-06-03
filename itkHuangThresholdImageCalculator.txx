
#ifndef __itkHuangThresholdImageCalculator_txx
#define __itkHuangThresholdImageCalculator_txx

#include "itkHuangThresholdImageCalculator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "vnl/vnl_math.h"

namespace itk
{ 
    
/**
 * Constructor
 */
template<class TInputImage>
HuangThresholdImageCalculator<TInputImage>
::HuangThresholdImageCalculator()
{
  m_Image = NULL;
  m_Threshold = NumericTraits<PixelType>::Zero;
  m_NumberOfHistogramBins = 128;
  m_RegionSetByUser = false;
}


/*
 * Compute the Huang's threshold
 */
template<class TInputImage>
void
HuangThresholdImageCalculator<TInputImage>
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

  // find first and last non-empty bin - could replace with stl
  int first, last;
  for (first = 0; (unsigned)first < relativeFrequency.size() && relativeFrequency[first] == 0; first++); // do nothing
  for (last = relativeFrequency.size() - 1; last > first && relativeFrequency[last] == 0; last--); // do nothing


  if (first == last)
    {
    itkWarningMacro(<< "No data in histogram");
    return;
    }
  // calculate the cumulative density and the weighted cumulative density
  std::vector<double> S(last+1);
  std::vector<double> W(last+1);

  S[0] = relativeFrequency[0];

  for (int i = std::max(1, first); i <= last; i++) 
    {
    S[i] = S[i - 1] + relativeFrequency[i];
    W[i] = W[i - 1] + i * relativeFrequency[i];
    }
  
  // precalculate the summands of the entropy given the absolute difference x - mu (integral)
  double C = last - first;
  std::vector<double> Smu(last + 1 - first);

  for (int i = 1; (unsigned)i < Smu.size(); i++) 
    {
    double mu = 1 / (1 + vcl_abs(i) / C);
    Smu[i] = -mu * vcl_log(mu) - (1 - mu) * vcl_log(1 - mu);
    }
  
  // calculate the threshold
  int bestThreshold = 0;
  double bestEntropy = itk::NumericTraits<double>::max();
  for (int threshold = first; threshold <= last; threshold++) 
    {
    double entropy = 0;
    int mu = (int)round(W[threshold] / S[threshold]);
    for (int i = first; i <= threshold; i++)
      entropy += Smu[vcl_abs(i - mu)] * relativeFrequency[i];
    mu = (int)round((W[last] - W[threshold]) / (S[last] - S[threshold]));
    for (int i = threshold + 1; i <= last; i++)
      entropy += Smu[vcl_abs(i - mu)] * relativeFrequency[i];
  
    if (bestEntropy > entropy) 
      {
      bestEntropy = entropy;
      bestThreshold = threshold;
      }
    }


  m_Threshold = static_cast<PixelType>( imageMin + 
                                        ( bestThreshold ) / binMultiplier );


}

template<class TInputImage>
void
HuangThresholdImageCalculator<TInputImage>
::SetRegion( const RegionType & region )
{
  m_Region = region;
  m_RegionSetByUser = true;
}

  
template<class TInputImage>
void
HuangThresholdImageCalculator<TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Threshold: " << m_Threshold << std::endl;
  os << indent << "NumberOfHistogramBins: " << m_NumberOfHistogramBins << std::endl;
  os << indent << "Image: " << m_Image.GetPointer() << std::endl;
}

} // end namespace itk

#endif
