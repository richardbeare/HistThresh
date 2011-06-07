#ifndef __itkTriangleThresholdImageCalculator_txx
#define __itkTriangleThresholdImageCalculator_txx

#include "itkTriangleThresholdImageCalculator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "vnl/vnl_math.h"

namespace itk
{ 
    
/**
 * Constructor
 */
template<class TInputImage>
TriangleThresholdImageCalculator<TInputImage>
::TriangleThresholdImageCalculator()
{
  m_Image = NULL;
  m_Threshold = NumericTraits<PixelType>::Zero;
  m_NumberOfHistogramBins = 128;
  m_RegionSetByUser = false;
  m_LowThresh=0.01;
  m_HighThresh=0.99;
}


/*
 * Compute the Triangle's threshold
 */
template<class TInputImage>
void
TriangleThresholdImageCalculator<TInputImage>
::Compute(void)
{

  unsigned int j;

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

  if (this->GetDebug())
    {
    for (unsigned i = 0;i < m_NumberOfHistogramBins; i++)
      {
      double count = relativeFrequency[i];
      double bin = ( imageMin + ( i + 1 ) / binMultiplier );
      std::cout << bin << "," << count << std::endl;
      }
    }
  // Triangle method needs the maximum and minimum indexes
  // Minimum indexes for this purpose are poorly defined - can't just
  // take a index with zero entries.
  double Mx = itk::NumericTraits<double>::min();
  unsigned long MxIdx=0;

  for ( j = 0; j < m_NumberOfHistogramBins; j++ )
    {
    //std::cout << relativeFrequency[j] << std::endl;
    if (relativeFrequency[j] > Mx)
      {
      MxIdx=j;
      Mx=relativeFrequency[j];
      }
    }


  cumSum[0]=relativeFrequency[0];
  for ( j = 1; j < m_NumberOfHistogramBins; j++ )
    {
    cumSum[j] = relativeFrequency[j] + cumSum[j-1];
    }


  double total = cumSum[m_NumberOfHistogramBins - 1];
  // find 1% and 99% levels
  double onePC = total * m_LowThresh;
  unsigned onePCIdx=0;
  for (j=0; j < m_NumberOfHistogramBins; j++ )
    {
    if (cumSum[j] > onePC)
      {
      onePCIdx = j;
      break;
      }
    }

  double nnPC = total * m_HighThresh;
  unsigned nnPCIdx=m_NumberOfHistogramBins;
  for (j=0; j < m_NumberOfHistogramBins; j++ )
    {
    if (cumSum[j] > nnPC)
      {
      nnPCIdx = j;
      break;
      }
    }


  
  // figure out which way we are looking - we want to construct our
  // line between the max index and the further of 1% and 99%
  unsigned ThreshIdx=0;
  if (fabs((float)MxIdx - (float)onePCIdx) > fabs((float)MxIdx - (float)nnPCIdx))
    {
    // line to 1 %
    double slope = Mx/(MxIdx - onePCIdx);
    for (unsigned k=onePCIdx; k<MxIdx; k++)
      {
      float line = (slope*(k-onePCIdx));
      triangle[k]= line - relativeFrequency[k];
      // std::cout << relativeFrequency[k] << "," << line << "," << triangle[k] << std::endl;
      }

    ThreshIdx = onePCIdx + std::distance(&(triangle[onePCIdx]), std::max_element(&(triangle[onePCIdx]), &(triangle[MxIdx]))) ;
    }
  else
    {
    // line to 99 %
    double slope = -Mx/(nnPCIdx - MxIdx);
    for (unsigned k=MxIdx; k < nnPCIdx; k++)
      {
      float line = (slope*(k-MxIdx) + Mx);
      triangle[k]= line - relativeFrequency[k];
//      std::cout << relativeFrequency[k] << "," << line << "," << triangle[k] << std::endl;
      }
    ThreshIdx = MxIdx + std::distance(&(triangle[MxIdx]), std::max_element(&(triangle[MxIdx]), &(triangle[nnPCIdx]))) ;
    }

  m_Threshold = static_cast<PixelType>( imageMin + 
                                        ( ThreshIdx + 1 ) / binMultiplier );


  // for (unsigned k = 0; k < m_NumberOfHistogramBins ; k++)
  //   {
  //   std::cout << relativeFrequency[k] << std::endl;
  //   }

}

template<class TInputImage>
void
TriangleThresholdImageCalculator<TInputImage>
::SetRegion( const RegionType & region )
{
  m_Region = region;
  m_RegionSetByUser = true;
}

  
template<class TInputImage>
void
TriangleThresholdImageCalculator<TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Threshold: " << m_Threshold << std::endl;
  os << indent << "NumberOfHistogramBins: " << m_NumberOfHistogramBins << std::endl;
  os << indent << "Image: " << m_Image.GetPointer() << std::endl;
}

} // end namespace itk

#endif
