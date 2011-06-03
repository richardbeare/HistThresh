
#ifndef __itkLiThresholdImageCalculator_txx
#define __itkLiThresholdImageCalculator_txx

#include "itkLiThresholdImageCalculator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "vnl/vnl_math.h"

namespace itk
{ 
    
/**
 * Constructor
 */
template<class TInputImage>
LiThresholdImageCalculator<TInputImage>
::LiThresholdImageCalculator()
{
  m_Image = NULL;
  m_Threshold = NumericTraits<PixelType>::Zero;
  m_NumberOfHistogramBins = 128;
  m_RegionSetByUser = false;
}


/*
 * Compute the Li's threshold
 */
template<class TInputImage>
void
LiThresholdImageCalculator<TInputImage>
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

  int threshold;
  int ih;
  int num_pixels;
  int sum_back; /* sum of the background pixels at a given threshold */
  int sum_obj;  /* sum of the object pixels at a given threshold */
  int num_back; /* number of background pixels at a given threshold */
  int num_obj;  /* number of object pixels at a given threshold */
  double old_thresh;
  double new_thresh;
  double mean_back; /* mean of the background pixels at a given threshold */
  double mean_obj;  /* mean of the object pixels at a given threshold */
  double mean;  /* mean gray-level in the image */
  double tolerance; /* threshold tolerance */
  double temp;

  tolerance=0.5;
  num_pixels = 0;
  for (ih = 0; (unsigned)ih < relativeFrequency.size(); ih++ )
    num_pixels += relativeFrequency[ih];

  /* Calculate the mean gray-level */
  mean = 0.0;
  for ( ih = 0; (unsigned)ih < relativeFrequency.size(); ih++ ) //0 + 1?
    mean += ih * relativeFrequency[ih];
  mean /= num_pixels;
  /* Initial estimate */
  new_thresh = mean;

  do{
  old_thresh = new_thresh;
  threshold = (int) (old_thresh + 0.5);	/* range */
  /* Calculate the means of background and object pixels */
  /* Background */
  sum_back = 0;
  num_back = 0;
  for ( ih = 0; ih <= threshold; ih++ ) 
    {
    sum_back += ih * relativeFrequency[ih];
    num_back += relativeFrequency[ih];
    }
  mean_back = ( num_back == 0 ? 0.0 : ( sum_back / ( double ) num_back ) );
  /* Object */
  sum_obj = 0;
  num_obj = 0;
  for ( ih = threshold + 1; (unsigned)ih < relativeFrequency.size(); ih++ ) 
    {
    sum_obj += ih * relativeFrequency[ih];
    num_obj += relativeFrequency[ih];
    }
  mean_obj = ( num_obj == 0 ? 0.0 : ( sum_obj / ( double ) num_obj ) );
  
  /* Calculate the new threshold: Equation (7) in Ref. 2 */
  //new_thresh = simple_round ( ( mean_back - mean_obj ) / ( Math.log ( mean_back ) - Math.log ( mean_obj ) ) );
  //simple_round ( double x ) {
  // return ( int ) ( IS_NEG ( x ) ? x - .5 : x + .5 );
  //}
  //
  //#define IS_NEG( x ) ( ( x ) < -DBL_EPSILON ) 
  //DBL_EPSILON = 2.220446049250313E-16
  temp = ( mean_back - mean_obj ) / ( vcl_log ( mean_back ) - vcl_log ( mean_obj ) );
  
  if (temp < -2.220446049250313E-16)
    new_thresh = (int) (temp - 0.5);
  else
    new_thresh = (int) (temp + 0.5);
  /*  Stop the iterations when the difference between the
			new and old threshold values is less than the tolerance */
  }
  while ( vcl_abs ( new_thresh - old_thresh ) > tolerance );

  m_Threshold = static_cast<PixelType>( imageMin + 
                                        ( threshold ) / binMultiplier );


}

template<class TInputImage>
void
LiThresholdImageCalculator<TInputImage>
::SetRegion( const RegionType & region )
{
  m_Region = region;
  m_RegionSetByUser = true;
}

  
template<class TInputImage>
void
LiThresholdImageCalculator<TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Threshold: " << m_Threshold << std::endl;
  os << indent << "NumberOfHistogramBins: " << m_NumberOfHistogramBins << std::endl;
  os << indent << "Image: " << m_Image.GetPointer() << std::endl;
}

} // end namespace itk

#endif
