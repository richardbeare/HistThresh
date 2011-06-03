#ifndef __itkKittlerIllingworthThresholdImageCalculator_txx
#define __itkKittlerIllingworthThresholdImageCalculator_txx

#include "itkKittlerIllingworthThresholdImageCalculator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "vnl/vnl_math.h"

namespace itk
{ 
    
/**
 * Constructor
 */
template<class TInputImage>
KittlerIllingworthThresholdImageCalculator<TInputImage>
::KittlerIllingworthThresholdImageCalculator()
{
  m_Image = NULL;
  m_Threshold = NumericTraits<PixelType>::Zero;
  m_NumberOfHistogramBins = 128;
  m_RegionSetByUser = false;
  }
template<class TInputImage>
int
KittlerIllingworthThresholdImageCalculator<TInputImage>
::Mean(std::vector<double> data)
{
  int threshold = -1;
  double tot=0, sum=0;
  for (unsigned i=0; i<data.size(); i++)
    {
    tot+= data[i];
    sum+=(i*data[i]);
    }
  threshold =(int) vcl_floor(sum/tot);
  return threshold;
  }

 
template<class TInputImage>
double 
KittlerIllingworthThresholdImageCalculator<TInputImage>
::A(std::vector<double> y, int j) 
{
  double x = 0;
  for (int i=0;i<=j;i++)
    x+=y[i];
  return x;
}

template<class TInputImage>
double 
KittlerIllingworthThresholdImageCalculator<TInputImage>
::B(std::vector<double> y, int j) 
{
  double x = 0;
  for (int i=0;i<=j;i++)
    x+=i*y[i];
  return x;
}

template<class TInputImage>
double 
KittlerIllingworthThresholdImageCalculator<TInputImage>
::C(std::vector<double> y, int j) 
{
  double x = 0;
  for (int i=0;i<=j;i++)
    x+=i*i*y[i];
  return x;
}


/*
 * Compute the KittlerIllingworth's threshold
 */
template<class TInputImage>
void
KittlerIllingworthThresholdImageCalculator<TInputImage>
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

  int threshold = Mean(relativeFrequency);
  int Tprev =-2;
  double mu, nu, p, q, sigma2, tau2, w0, w1, w2, sqterm, temp;
  //int counter=1;
  while (threshold!=Tprev)
    {
    //Calculate some statistics.
    mu = B(relativeFrequency, threshold)/A(relativeFrequency, threshold);
    nu = (B(relativeFrequency, relativeFrequency.size() - 1)-B(relativeFrequency, threshold))/(A(relativeFrequency, relativeFrequency.size() - 1)-A(relativeFrequency, threshold));
    p = A(relativeFrequency, threshold)/A(relativeFrequency, relativeFrequency.size() - 1);
    q = (A(relativeFrequency, relativeFrequency.size() - 1)-A(relativeFrequency, threshold)) / A(relativeFrequency, relativeFrequency.size() - 1);
    sigma2 = C(relativeFrequency, threshold)/A(relativeFrequency, threshold)-(mu*mu);
    tau2 = (C(relativeFrequency, relativeFrequency.size() - 1)-C(relativeFrequency, threshold)) / (A(relativeFrequency, relativeFrequency.size() - 1)-A(relativeFrequency, threshold)) - (nu*nu);
    
    //The terms of the quadratic equation to be solved.
    w0 = 1.0/sigma2-1.0/tau2;
    w1 = mu/sigma2-nu/tau2;
    w2 = (mu*mu)/sigma2 - (nu*nu)/tau2 + vcl_log10((sigma2*(q*q))/(tau2*(p*p)));
    
    //If the next threshold would be imaginary, return with the current one.
    sqterm = (w1*w1)-w0*w2;
    if (sqterm < 0) 
      {
      itkWarningMacro( << "MinError(I): not converging. Try \'Ignore black/white\' options");
      m_Threshold=threshold;
      return;
      }
  
    //The updated threshold is the integer part of the solution of the quadratic equation.
    Tprev = threshold;
    temp = (w1+vcl_sqrt(sqterm))/w0;
  
    if (vnl_math_isnan(temp)) 
      {
      itkWarningMacro (<< "MinError(I): NaN, not converging. Try \'Ignore black/white\' options");
      threshold = Tprev;
      }
    else
      {
      threshold =(int) vcl_floor(temp);
      }
  }
  m_Threshold = static_cast<PixelType>( imageMin + 
					( threshold) / binMultiplier );


}

template<class TInputImage>
void
KittlerIllingworthThresholdImageCalculator<TInputImage>
::SetRegion( const RegionType & region )
{
  m_Region = region;
  m_RegionSetByUser = true;
}

  
template<class TInputImage>
void
KittlerIllingworthThresholdImageCalculator<TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Threshold: " << m_Threshold << std::endl;
  os << indent << "NumberOfHistogramBins: " << m_NumberOfHistogramBins << std::endl;
  os << indent << "Image: " << m_Image.GetPointer() << std::endl;
}

} // end namespace itk

#endif
