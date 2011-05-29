/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTriangleThresholdImageCalculator.txx,v $
  Language:  C++
  Date:      $Date: 2009-01-26 21:45:54 $
  Version:   $Revision: 1.9 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
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

  // Triangle method needs the maximum and minimum indexes
  // Minimum indexes for this purpose are poorly defined - can't just
  // take a index with zero entries.
  double Mx = itk::NumericTraits<double>::min();
  unsigned long MxIdx=0;
  //std::cout << "-------" << std::endl;
  for ( j = 0; j < m_NumberOfHistogramBins; j++ )
    {
    //std::cout << relativeFrequency[j] << std::endl;
    if (relativeFrequency[j] > Mx)
      {
      MxIdx=j;
      Mx=relativeFrequency[j];
      }
    }
  //std::cout << "-------" << std::endl;

  cumSum[0]=relativeFrequency[0];
  for ( j = 1; j < m_NumberOfHistogramBins; j++ )
    {
    cumSum[j] = relativeFrequency[j] + cumSum[j-1];
    }

  std::cout << MxIdx << " " << Mx << std::endl;
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

  std::cout << onePCIdx << " " << nnPCIdx << std::endl;
  
  // figure out which way we are looking - we want to construct our
  // line between the max index and the further of 1% and 99%
  unsigned ThreshIdx=0;
  if (fabs(MxIdx - onePCIdx) > fabs(MxIdx - nnPCIdx))
    {
    // line to 1 %
    double slope = Mx/(MxIdx - onePCIdx);
    for (unsigned k=0; k<MxIdx; k++)
      {
      triangle[k]=slope*(k-MxIdx) + Mx - relativeFrequency[k];
      }

    ThreshIdx = onePCIdx + std::distance(&(triangle[onePCIdx]), std::max_element(&(triangle[onePCIdx]), &(triangle[MxIdx]))) ;
    }
  else
    {
    // line to 99 %
    double slope = -Mx/(nnPCIdx - MxIdx);
    std::cout << "Here " << slope << std::endl;
    std::cout << "------------" << std::endl;
    for (unsigned k=MxIdx; k < nnPCIdx; k++)
      {
      std::cout << slope*(k-MxIdx) + Mx - relativeFrequency[k] << std::endl;
      triangle[k]=(slope*(k-MxIdx) + Mx) - relativeFrequency[k];
      }
    std::cout << "------------" << std::endl;

    ThreshIdx = MxIdx + std::distance(&(triangle[MxIdx]), std::max_element(&(triangle[MxIdx]), &(triangle[nnPCIdx]))) ;
    }

  std::cout << ThreshIdx << " " << imageMin << " " << binMultiplier << std::endl;
  m_Threshold = static_cast<PixelType>( imageMin + 
                                        ( ThreshIdx + 1 ) / binMultiplier );

#if 0 
  // normalize the frequencies
  double totalMean = 0.0;
  for ( j = 0; j < m_NumberOfHistogramBins; j++ )
    {
    relativeFrequency[j] /= totalPixels;
    totalMean += (j+1) * relativeFrequency[j];
    }


  // compute Otsu's threshold by maximizing the between-class
  // variance
  double freqLeft = relativeFrequency[0];
  double meanLeft = 1.0;
  double meanRight = ( totalMean - freqLeft ) / ( 1.0 - freqLeft );

  double maxVarBetween = freqLeft * ( 1.0 - freqLeft ) *
    vnl_math_sqr( meanLeft - meanRight );
  int maxBinNumber = 0;

  double freqLeftOld = freqLeft;
  double meanLeftOld = meanLeft;

  for ( j = 1; j < m_NumberOfHistogramBins; j++ )
    {
    freqLeft += relativeFrequency[j];
    meanLeft = ( meanLeftOld * freqLeftOld + 
                 (j+1) * relativeFrequency[j] ) / freqLeft;
    if (freqLeft == 1.0)
      {
      meanRight = 0.0;
      }
    else
      {
      meanRight = ( totalMean - meanLeft * freqLeft ) / 
        ( 1.0 - freqLeft );
      }
    double varBetween = freqLeft * ( 1.0 - freqLeft ) *
      vnl_math_sqr( meanLeft - meanRight );
   
    if ( varBetween > maxVarBetween )
      {
      maxVarBetween = varBetween;
      maxBinNumber = j;
      }

    // cache old values
    freqLeftOld = freqLeft;
    meanLeftOld = meanLeft; 

    } 

  m_Threshold = static_cast<PixelType>( imageMin + 
                                        ( maxBinNumber + 1 ) / binMultiplier );

#endif

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
